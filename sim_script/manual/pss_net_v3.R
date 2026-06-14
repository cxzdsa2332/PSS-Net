rm(list = ls())

################################################################################
# pss_net_v3.R  —  PSS-Net: 8-species microbial community
#                   Sparse additive nonparametric ODE  (CLAUDE.md framework)
#
# Model:    dx_j/dt = μ_j + Σ_i f_{ji}(x_i) + u_j
# Basis:    ψ(x) = (x, x², …, x^M)  [satisfies f_{ji}(0)=0, CLAUDE.md eq.2]
# At SS (centered):  Ψ_c · θ_j = −u_{c,j}
# Intercept:         μ̂_j = −ū_j − Ψ̄ · θ̂_j          (CLAUDE.md eq.10)
# Jacobian:          J_{ji} = ψ'(x_i^wt)ᵀ θ̂_{ji}     (CLAUDE.md eq.11)
# Edge detection:    ‖θ̂_{ji,·}‖₂ ≥ τ                  (CLAUDE.md §7.3)
#
# Simulation: linear GLV  (true f_{ji}(x) = A_{ji}·x, so M=2 recovers it)
# Regression: ADSIHT (primary) + Group Lasso (secondary)
################################################################################

library(deSolve)
library(ADSIHT)
library(grpreg)   # group Lasso/MCP/SCAD, more stable CV than gglasso

set.seed(42)

# ── Tuning ────────────────────────────────────────────────────────────────────
M        <- 2    # polynomial degree: ψ(x) = (x, x², …, x^M)
n_sp     <- 8
sp_names <- paste0("Sp", seq_len(n_sp))

# ── True parameters ───────────────────────────────────────────────────────────
r_true   <- c(0.8, 1.2, 0.6, 1.0, 0.7, 1.1, 0.5, 0.9)
gam_true <- c(1.5, 1.8, 1.2, 1.6, 1.4, 1.5, 1.0, 1.6)

A_true <- matrix(0, n_sp, n_sp)
A_true[1,3] <-  0.40; A_true[1,5] <- -0.30
A_true[2,1] <-  0.30; A_true[2,4] <- -0.40; A_true[2,7] <-  0.20
A_true[3,2] <- -0.30; A_true[3,6] <-  0.30
A_true[4,1] <-  0.40; A_true[4,3] <- -0.20; A_true[4,8] <-  0.30
A_true[5,2] <-  0.30; A_true[5,6] <- -0.20
A_true[6,4] <-  0.20; A_true[6,5] <-  0.30
A_true[7,3] <- -0.30; A_true[7,8] <-  0.40
A_true[8,1] <-  0.20; A_true[8,6] <- -0.20

n_edges <- sum(A_true != 0)

cat("================================================================\n")
cat("  sindy_ss_v3 (PSS-Net): 8-species microbial community\n")
cat(sprintf("  %d species | %d true edges | sparsity %.1f%% | M=%d\n",
            n_sp, n_edges, 100 * n_edges / (n_sp*(n_sp-1)), M))
cat("================================================================\n\n")

# ── ODE (linear GLV) ──────────────────────────────────────────────────────────
ode_func <- function(t, state, parms) {
  x <- pmax(state, 0)
  list(r_true + as.numeric(A_true %*% x) - gam_true * x + parms$u)
}

# ── Step 1: Simulate perturbed steady states ──────────────────────────────────
cat("Step 1: Simulating perturbed steady states...\n")

N_cond <- 300
U_mat  <- matrix(runif(N_cond * n_sp, -0.4, 0.8), N_cond, n_sp)
U_mat[1, ] <- 0

ss_mat <- matrix(NA, N_cond, n_sp)
for (k in seq_len(N_cond)) {
  tryCatch({
    out <- ode(runif(n_sp, 0.5, 2.0), c(0, 5000), ode_func,
               list(u = U_mat[k, ]), method = "lsoda",
               rtol = 1e-9, atol = 1e-11)
    ss_mat[k, ] <- out[nrow(out), 2:(n_sp + 1)]
  }, error = function(e) NULL)
}

ok    <- apply(ss_mat, 1, function(r) all(is.finite(r) & r > 0))
X_obs <- ss_mat[ok, ] + matrix(rnorm(sum(ok)*n_sp, 0, 0.03), sum(ok), n_sp)
X_obs <- pmax(X_obs, 1e-6)
U_obs <- U_mat[ok, ]
N     <- nrow(X_obs)

# x_wt: noise-free WT steady state via independent ODE integration (u=0).
# Do NOT use X_obs[which(rowSums(...)==0), ] — that row contains σ=0.03 noise
# which biases polynomial basis evaluation and Jacobian estimation.
out_wt <- ode(rep(1, n_sp), c(0, 1e4), ode_func, list(u = rep(0, n_sp)),
              method = "lsoda", rtol = 1e-12, atol = 1e-14)
x_wt   <- as.numeric(out_wt[nrow(out_wt), 2:(n_sp + 1)])
cat(sprintf("  Valid conditions: %d / %d\n", N, N_cond))
cat(sprintf("  x_wt (noise-free): %s\n\n",
            paste(sprintf("%.3f", x_wt), collapse = ", ")))

# ── Step 2: Polynomial basis  ψ(x) = (x, x², …, x^M) ────────────────────────
# Monomial basis satisfies ψ_m(0)=0 (CLAUDE.md eq.2) and maps the linear GLV
# true model f_{ji}(x)=A_{ji}·x onto a single nonzero coefficient (m=1 only),
# preserving the group-sparse structure that ADSIHT/grLasso exploit.
# Column standardization below is sufficient to control condition number at M=2.
cat(sprintf("Step 2: Building polynomial basis (M=%d)...\n", M))

Psi <- matrix(0, N, n_sp * M)
for (i in seq_len(n_sp))
  for (m in seq_len(M))
    Psi[, (i - 1) * M + m] <- X_obs[, i]^m

group    <- rep(seq_len(n_sp), each = M)
Psi_bar  <- colMeans(Psi)
Psi_c    <- sweep(Psi, 2, Psi_bar)
# Column standardization: equalizes scale across polynomial degrees
Psi_sd   <- apply(Psi_c, 2, sd)
Psi_sd   <- pmax(Psi_sd, 1e-10)   # guard against constant columns
Psi_cs   <- sweep(Psi_c, 2, Psi_sd, "/")
U_bar    <- colMeans(U_obs)
U_c      <- sweep(U_obs, 2, U_bar)

cat(sprintf("  Design matrix: %d × %d  (%d groups × %d bases)\n\n",
            N, ncol(Psi), n_sp, M))

# ── Step 3: Sparse regression per target species ──────────────────────────────
cat("Step 3: Sparse regression...\n")

ALPHA_ads <- matrix(0, n_sp, n_sp * M)
ALPHA_gl  <- matrix(0, n_sp, n_sp * M)
MU_ads    <- numeric(n_sp)
MU_gl     <- numeric(n_sp)

for (j in seq_len(n_sp)) {
  rhs_c <- -U_c[, j]

  # ADSIHT on standardized matrix; rescale coefficients back to original scale
  fit_a         <- ADSIHT(Psi_cs, matrix(rhs_c), group, ic.type = "dsic")
  best_a        <- which.min(fit_a$ic)
  ALPHA_ads[j,] <- fit_a$beta[, best_a] / Psi_sd   # undo standardization
  MU_ads[j]     <- -U_bar[j] - sum(Psi_bar * ALPHA_ads[j, ])

  # grpreg on standardized matrix; coef() returns lambda.min coefficients
  cv            <- cv.grpreg(Psi_cs, rhs_c, group = group,
                             penalty = "grLasso", nfolds = 5)
  ALPHA_gl[j, ] <- coef(cv)[-1] / Psi_sd
  MU_gl[j]      <- -U_bar[j] - sum(Psi_bar * ALPHA_gl[j, ])
}
cat("  Done.\n\n")

# ── Step 4: Jacobian  J_{ji} = ψ'(x_i^wt)ᵀ θ̂_{ji}  ─────────────────────────
# ψ'_m(x) = m · x^{m-1}  (analytic derivative of monomial basis)
get_jacobian <- function(alpha_mat) {
  j_mat <- matrix(0, n_sp, n_sp)
  for (j in seq_len(n_sp)) {
    for (i in seq_len(n_sp)) {
      cols        <- (i - 1) * M + seq_len(M)
      dpsi_i      <- seq_len(M) * x_wt[i]^pmax(seq_len(M) - 1, 0)
      j_mat[j, i] <- sum(alpha_mat[j, cols] * dpsi_i)
    }
  }
  j_mat
}

J_ads <- get_jacobian(ALPHA_ads)
J_gl  <- get_jacobian(ALPHA_gl)

get_A   <- function(J) { A <- J; diag(A) <- 0; A }
get_gam <- function(J) -diag(J)

A_ads <- get_A(J_ads);  g_ads <- get_gam(J_ads)
A_gl  <- get_A(J_gl);   g_gl  <- get_gam(J_gl)

# ── Step 5: Binary adjacency + metrics ────────────────────────────────────────
# Edge (j←i) ∈ Ê  iff  ‖θ̂_{ji,·}‖₂ ≥ τ  (CLAUDE.md §7.3, group L2 norm).
#
# Threshold choice:
#   ADSIHT  — IHT produces exact group zeros; τ_ads = 1e-10 recovers the
#              exact sparsity pattern without shrinkage bias.
#   grLasso — continuous shrinkage never gives exact zeros; use data-adaptive
#              τ_gl = max(‖θ̂_{ji}‖₂) * rel_thr to accommodate scale variation.
tau_ads <- 1e-10   # float-safe exact-zero detection for ADSIHT
rel_thr <- 0.01    # grLasso: edge present if group norm ≥ 1% of maximum norm

group_norm_mat <- function(ALPHA) {
  mat <- matrix(0, n_sp, n_sp)
  for (j in seq_len(n_sp))
    for (i in seq_len(n_sp)) {
      cols      <- (i - 1) * M + seq_len(M)
      mat[j, i] <- sqrt(sum(ALPHA[j, cols]^2))
    }
  mat
}

norm_ads <- group_norm_mat(ALPHA_ads)
norm_gl  <- group_norm_mat(ALPHA_gl)

# Adaptive threshold for grLasso: relative to largest group norm per row
tau_gl_vec <- apply(norm_gl, 1, max) * rel_thr

get_adj_ads <- function(norm_mat, tau) {
  adj <- norm_mat >= tau
  diag(adj) <- FALSE
  adj * 1L
}

get_adj_gl <- function(norm_mat, tau_vec) {
  adj <- matrix(FALSE, n_sp, n_sp)
  for (j in seq_len(n_sp))
    adj[j, ] <- norm_mat[j, ] >= tau_vec[j]
  diag(adj) <- FALSE
  adj * 1L
}

adj_true <- (A_true != 0) * 1L
adj_ads  <- get_adj_ads(norm_ads, tau_ads)
adj_gl   <- get_adj_gl(norm_gl,  tau_gl_vec)

# True coefficient matrix: θ_true[j,(i-1)M+1] = A_true[j,i]; higher = 0
alpha_true <- matrix(0, n_sp, n_sp * M)
for (j in seq_len(n_sp))
  for (i in seq_len(n_sp))
    alpha_true[j, (i - 1) * M + 1] <- A_true[j, i]

mets <- function(adj_h, ALPHA_h) {
  TP <- sum(adj_true == 1 & adj_h == 1)
  FP <- sum(adj_true == 0 & adj_h == 1)
  FN <- sum(adj_true == 1 & adj_h == 0)
  TN <- sum(adj_true == 0 & adj_h == 0)
  pr  <- TP / (TP + FP + 1e-9)
  re  <- TP / (TP + FN + 1e-9)
  f1  <- 2 * pr * re / (pr + re + 1e-9)
  denom <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  mcc <- if (denom < 1e-9) 0 else (TP * TN - FP * FN) / denom
  coef_l2 <- mean(sqrt(rowSums((ALPHA_h - alpha_true)^2)))
  list(TP = TP, FP = FP, FN = FN, TN = TN,
       Pr = pr, Re = re, F1 = f1, mcc = mcc,
       coef_l2 = coef_l2,
       jac_rmse = sqrt(mean((get_A(get_jacobian(ALPHA_h)) - A_true)^2)))
}

m_ads <- mets(adj_ads, ALPHA_ads)
m_gl  <- mets(adj_gl,  ALPHA_gl)

# ── Print adjacency detail ─────────────────────────────────────────────────────
print_adj <- function(label, adj_h) {
  cat(sprintf("\n%s  adjacency (row=target j, col=source i):\n", label))
  mat <- matrix(sprintf("%d", adj_h), n_sp, n_sp,
                dimnames = list(sp_names, sp_names))
  print(noquote(mat))
  tp_idx <- which(adj_true == 1 & adj_h == 1, arr.ind = TRUE)
  fp_idx <- which(adj_true == 0 & adj_h == 1, arr.ind = TRUE)
  fn_idx <- which(adj_true == 1 & adj_h == 0, arr.ind = TRUE)
  fmt_pairs <- function(idx)
    if (nrow(idx) == 0) "none"
    else paste(sprintf("%s<-%s", sp_names[idx[,1]], sp_names[idx[,2]]),
               collapse = ", ")
  cat(sprintf("  TP: %s\n", fmt_pairs(tp_idx)))
  cat(sprintf("  FP: %s\n", fmt_pairs(fp_idx)))
  cat(sprintf("  FN: %s\n", fmt_pairs(fn_idx)))
}

mu_true <- r_true   # f_{ji}(0)=0 ⟹ μ_j = r_j  (CLAUDE.md eq.2)

# ── Parameter recovery table ──────────────────────────────────────────────────
cat("\n── r (intrinsic growth) recovery ──────────────────────────────────\n")
cat(sprintf("  %-6s  %7s  %9s  %9s\n", "Sp", "r_true", "r_ADSIHT", "r_grLasso"))
for (j in seq_len(n_sp))
  cat(sprintf("  %-6s  %7.4f  %9.4f  %9.4f\n",
              sp_names[j], r_true[j], MU_ads[j], MU_gl[j]))
cat(sprintf("  RMSE   ADSIHT=%.4f  GroupLasso=%.4f\n",
            sqrt(mean((MU_ads - mu_true)^2)),
            sqrt(mean((MU_gl  - mu_true)^2))))

cat("\n── gamma (death/self-regulation) recovery ──────────────────────────\n")
cat(sprintf("  %-6s  %9s  %9s  %9s\n",
            "Sp", "gam_true", "g_ADSIHT", "g_grLasso"))
for (j in seq_len(n_sp))
  cat(sprintf("  %-6s  %9.4f  %9.4f  %9.4f\n",
              sp_names[j], gam_true[j], g_ads[j], g_gl[j]))
cat(sprintf("  RMSE   ADSIHT=%.4f  GroupLasso=%.4f\n",
            sqrt(mean((g_ads - gam_true)^2)),
            sqrt(mean((g_gl  - gam_true)^2))))

print_adj("ADSIHT",      adj_ads)
print_adj("Group-Lasso", adj_gl)

cat("\n================================================================\n")
cat(sprintf("  Polynomial basis M=%d  |  N=%d conditions\n", M, N))
cat(sprintf("  True edges: %d / %d off-diagonal\n",
            n_edges, n_sp * (n_sp - 1)))
fmt_row <- function(label, m)
  cat(sprintf(
    "  %-12s  TP=%2d FP=%2d FN=%2d TN=%2d  Pr=%.3f Re=%.3f F1=%.3f MCC=%.3f  CoefL2=%.4f JacRMSE=%.4f\n",
    label, m$TP, m$FP, m$FN, m$TN,
    m$Pr, m$Re, m$F1, m$mcc, m$coef_l2, m$jac_rmse
  ))

fmt_row("ADSIHT",      m_ads)
fmt_row("Group-Lasso", m_gl)
cat("================================================================\n")

# ── Step 6: Dynamics comparison — true ODE vs ADSIHT-reconstructed ODE ────────
library(ggplot2)

# ψ(x[i]) = (x[i]^1, …, x[i]^M): monomial basis evaluated at current state
ode_recon <- function(t, state, parms) {
  x     <- pmax(state, 0)
  u     <- parms$u
  alpha <- parms$alpha
  mu    <- parms$mu
  dx    <- numeric(n_sp)
  for (j in seq_len(n_sp)) {
    fhat_sum <- 0
    for (i in seq_len(n_sp)) {
      cols     <- (i - 1) * M + seq_len(M)
      fhat_sum <- fhat_sum + sum(alpha[j, cols] * x[i]^seq_len(M))
    }
    dx[j] <- mu[j] + fhat_sum + u[j]
  }
  list(dx)
}

u_tests <- list(
  "u: Sp1=-0.5, Sp5=+0.6" = { u <- rep(0, n_sp); u[1] <- -0.5; u[5] <-  0.6; u },
  "u: Sp3=+0.5, Sp7=-0.4" = { u <- rep(0, n_sp); u[3] <-  0.5; u[7] <- -0.4; u },
  "u: all +0.3"            = { rep(0.3, n_sp) }
)

t_span <- seq(0, 200, by = 0.5)

wide_to_long <- function(df_wide, model_label, perturb_label) {
  do.call(rbind, lapply(sp_names, function(sp) {
    data.frame(
      time    = df_wide$time,
      species = sp,
      x       = df_wide[[sp]],
      model   = model_label,
      perturb = perturb_label,
      stringsAsFactors = FALSE
    )
  }))
}

traj_df <- do.call(rbind, lapply(names(u_tests), function(uname) {
  u_vec <- u_tests[[uname]]

  out_true <- ode(x_wt, t_span, ode_func, list(u = u_vec),
                  method = "lsoda", rtol = 1e-8, atol = 1e-10)
  df_true           <- as.data.frame(out_true)
  colnames(df_true) <- c("time", sp_names)

  out_recon <- tryCatch(
    ode(x_wt, t_span, ode_recon,
        list(u = u_vec, alpha = ALPHA_ads, mu = MU_ads),
        method = "lsoda", rtol = 1e-8, atol = 1e-10),
    error = function(e) NULL
  )
  if (is.null(out_recon)) return(NULL)
  df_recon           <- as.data.frame(out_recon)
  colnames(df_recon) <- c("time", sp_names)

  rbind(wide_to_long(df_true,  "True",   uname),
        wide_to_long(df_recon, "ADSIHT", uname))
}))

traj_df$species <- factor(traj_df$species, levels = sp_names)
traj_df$model   <- factor(traj_df$model,   levels = c("True", "ADSIHT"))
traj_df$perturb <- factor(traj_df$perturb, levels = names(u_tests))

wt_df <- data.frame(species = factor(sp_names, levels = sp_names),
                    xwt     = x_wt)

p_dyn <- ggplot(traj_df, aes(x = time, y = x,
                             color = model, linetype = model)) +
  geom_hline(data = wt_df, aes(yintercept = xwt),
             color = "grey80", linewidth = 0.35, linetype = "dotted") +
  geom_line(linewidth = 0.75) +
  facet_grid(species ~ perturb, scales = "free_y") +
  scale_color_manual(values = c(True = "grey20", ADSIHT = "tomato3")) +
  scale_linetype_manual(values = c(True = "solid", ADSIHT = "dashed")) +
  labs(
    title    = sprintf(
      "Dynamics comparison: true ODE vs ADSIHT-reconstructed ODE  [M=%d]", M),
    subtitle = "Dotted line = WT steady state  |  starting from x_wt",
    x = "Time", y = "Abundance", color = NULL, linetype = NULL
  ) +
  theme_bw(base_size = 8.5) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    strip.text.x     = element_text(size = 7,   face = "bold"),
    strip.text.y     = element_text(size = 7.5, face = "bold"),
    legend.position  = "bottom",
    plot.title       = element_text(size = 10, face = "bold"),
    plot.subtitle    = element_text(size = 7.5, color = "grey45")
  )

print(p_dyn)

# ── Step 7: Effect decomposition — true vs ADSIHT, k selected species ─────────
k_targets   <- c(1, 2, 4)
decomp_pert <- names(u_tests)[[1]]

sub_idx   <- traj_df$perturb == decomp_pert
traj_sub  <- traj_df[sub_idx, ]
time_vals <- sort(unique(traj_sub$time))
n_t       <- length(time_vals)

xi_true_mat  <- matrix(NA, n_t, n_sp, dimnames = list(NULL, sp_names))
xi_recon_mat <- matrix(NA, n_t, n_sp, dimnames = list(NULL, sp_names))
for (sp in sp_names) {
  idx_t <- traj_sub$model == "True"   & traj_sub$species == sp
  idx_r <- traj_sub$model == "ADSIHT" & traj_sub$species == sp
  xi_true_mat[, sp]  <- traj_sub$x[idx_t][order(traj_sub$time[idx_t])]
  xi_recon_mat[, sp] <- traj_sub$x[idx_r][order(traj_sub$time[idx_r])]
}

decomp_rows <- vector("list", length(k_targets) * n_sp)
row_idx <- 1L
for (j in k_targets) {
  for (i in seq_len(n_sp)) {
    cols     <- (i - 1) * M + seq_len(M)
    xi_t     <- xi_true_mat[,  sp_names[i]]
    xi_r     <- xi_recon_mat[, sp_names[i]]
    f_true   <- A_true[j, i] * xi_t
    theta_ji <- ALPHA_ads[j, cols]
    psi_mat  <- sapply(seq_len(M), function(m) xi_r^m)
    f_hat    <- as.numeric(psi_mat %*% theta_ji)
    decomp_rows[[row_idx]] <- data.frame(
      time   = rep(time_vals, 2L),
      source = sp_names[i],
      target = sp_names[j],
      f      = c(f_true, f_hat),
      model  = rep(c("True", "ADSIHT"), each = n_t)
    )
    row_idx <- row_idx + 1L
  }
}
decomp_df        <- do.call(rbind, decomp_rows)
decomp_df$source <- factor(decomp_df$source, levels = sp_names)
decomp_df$target <- factor(decomp_df$target, levels = sp_names[k_targets])
decomp_df$model  <- factor(decomp_df$model,  levels = c("True", "ADSIHT"))

p_decomp <- ggplot(decomp_df,
                   aes(x = time, y = f,
                       color = model, linetype = model)) +
  geom_hline(yintercept = 0, color = "grey75", linewidth = 0.3) +
  geom_line(linewidth = 0.7) +
  facet_grid(target ~ source, scales = "free_y",
             labeller = labeller(
               target = function(x) paste0("Target: ", x),
               source = function(x) paste0("Source: ", x)
             )) +
  scale_color_manual(values = c(True = "grey25", ADSIHT = "tomato3")) +
  scale_linetype_manual(values = c(True = "solid", ADSIHT = "dashed")) +
  labs(
    title    = sprintf(
      "Effect decomposition: f_{ji}(x_i(t))  |  targets: %s  |  M=%d",
      paste(sp_names[k_targets], collapse = ", "), M),
    subtitle = sprintf("Perturbation: %s", decomp_pert),
    x = "Time", y = expression(f[ji](x[i](t))),
    color = NULL, linetype = NULL
  ) +
  theme_bw(base_size = 8) +
  theme(
    strip.background = element_rect(fill = "grey93"),
    strip.text.x     = element_text(size = 6.5, face = "bold"),
    strip.text.y     = element_text(size = 6.5, face = "bold"),
    legend.position  = "bottom",
    plot.title       = element_text(size = 9,  face = "bold"),
    plot.subtitle    = element_text(size = 7.5, color = "grey45")
  )

print(p_decomp)

# ── Step 8: Network comparison — true vs ADSIHT (igraph) ──────────────────────
# Two side-by-side directed graphs sharing the same circular layout.
# Edge colour encodes sign of interaction (A_true for true, A_ads for inferred).
# Edge classification: TP (correctly recovered), FP (spurious), FN (missed).
library(igraph)

# ── Build graphs ──────────────────────────────────────────────────────────────
# igraph edge list: from = source i (col), to = target j (row)
make_graph_from_adj <- function(adj_mat, weight_mat) {
  el <- which(adj_mat != 0, arr.ind = TRUE)   # [row=j, col=i]
  if (nrow(el) == 0) return(make_empty_graph(n_sp, directed = TRUE))
  g <- graph_from_edgelist(cbind(el[, 2], el[, 1]), directed = TRUE)
  E(g)$weight <- weight_mat[el]
  V(g)$name   <- sp_names
  g
}

g_true <- make_graph_from_adj(adj_true, A_true)
g_ads  <- make_graph_from_adj(adj_ads,  A_ads)

# Shared circular layout (consistent node positions)
lay <- layout_in_circle(g_true)

# ── Edge colour: positive (promotion) = tomato3, negative (inhibition) = steelblue3
edge_col <- function(w) ifelse(w > 0, "tomato3", "steelblue3")

# ── Edge classification for ADSIHT graph ──────────────────────────────────────
# Each ADSIHT edge (i→j) is labelled TP, FP, or FN (missed = in true not ads)
# We draw TP/FP edges in g_ads, and FN edges as dashed overlays on g_true layout
tp_mat <- adj_true == 1 & adj_ads == 1
fp_mat <- adj_true == 0 & adj_ads == 1
fn_mat <- adj_true == 1 & adj_ads == 0

# ── Plot ──────────────────────────────────────────────────────────────────────
par(mfrow = c(1, 2), mar = c(1, 1, 2.5, 1), bg = "white")

# Panel 1: True network
e_col_true <- edge_col(E(g_true)$weight)
plot(g_true,
     layout         = lay,
     vertex.color   = "grey88",
     vertex.frame.color = "grey55",
     vertex.size    = 26,
     vertex.label   = sp_names,
     vertex.label.cex   = 0.85,
     vertex.label.color = "grey10",
     edge.color     = e_col_true,
     edge.width     = 2,
     edge.arrow.size = 0.55,
     edge.curved    = 0.25,
     main           = "True network")
legend("bottomleft", bty = "n", cex = 0.75,
       legend = c("Promotion (+)", "Inhibition (−)"),
       col    = c("tomato3", "steelblue3"),
       lwd    = 2)

# Panel 2: ADSIHT inferred network
# ends() returns integer vertex IDs; index tp_mat with integers directly
ads_el    <- ends(g_ads, E(g_ads), names = FALSE)  # n_edge × 2 integer matrix
ads_class <- ifelse(
  mapply(function(src, tgt) tp_mat[tgt, src], ads_el[, 1], ads_el[, 2]),
  "TP", "FP"
)
ads_w    <- E(g_ads)$weight
ads_ecol <- ifelse(ads_class == "TP",
                   ifelse(ads_w > 0, "tomato3", "steelblue3"),
                   "orange2")   # FP always orange

plot(g_ads,
     layout         = lay,
     vertex.color   = "grey88",
     vertex.frame.color = "grey55",
     vertex.size    = 26,
     vertex.label   = sp_names,
     vertex.label.cex   = 0.85,
     vertex.label.color = "grey10",
     edge.color     = ads_ecol,
     edge.width     = ifelse(ads_class == "FP", 1.2, 2),
     edge.lty       = ifelse(ads_class == "FP", 2,   1),
     edge.arrow.size = 0.55,
     edge.curved    = 0.25,
     main           = sprintf(
       "ADSIHT inferred  [TP=%d FP=%d FN=%d  MCC=%.2f]",
       m_ads$TP, m_ads$FP, m_ads$FN, m_ads$mcc))

# Overlay FN edges as dotted grey arrows
# Use g_true as base (same n_sp nodes) and hide all edges then re-add FN only
fn_el <- which(fn_mat, arr.ind = TRUE)   # [row=j, col=i]
if (nrow(fn_el) > 0) {
  g_fn <- make_empty_graph(n_sp, directed = TRUE)
  g_fn <- add_edges(g_fn, as.vector(t(cbind(fn_el[, 2], fn_el[, 1]))))
  plot(g_fn, layout = lay, add = TRUE,
       vertex.color       = NA,
       vertex.frame.color = NA,
       vertex.label       = NA,
       vertex.size        = 26,
       edge.color         = "grey60",
       edge.width         = 1.2,
       edge.lty           = 3,
       edge.arrow.size    = 0.45,
       edge.curved        = 0.25)
}

legend("bottomleft", bty = "n", cex = 0.75,
       legend = c("TP promotion", "TP inhibition", "FP", "FN (missed)"),
       col    = c("tomato3", "steelblue3", "orange2", "grey60"),
       lwd    = c(2, 2, 1.2, 1.2),
       lty    = c(1, 1, 2, 3))

par(mfrow = c(1, 1))

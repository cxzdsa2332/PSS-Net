rm(list = ls())

################################################################################
# pss_net_v0.R  —  PSS-Net: 8-species microbial community
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
x_wt  <- X_obs[1, ]   # wild-type (u=0) steady state
cat(sprintf("  Valid conditions: %d / %d\n\n", N, N_cond))

# ── Step 2: Polynomial basis  ψ(x) = (x, x², …, x^M) ─────────────────────────
# Each source species i contributes M columns; group i = columns (i-1)M+1 : iM
# All basis functions satisfy ψ_m(0)=0  (CLAUDE.md eq.2)
cat(sprintf("Step 2: Building polynomial basis (M=%d)...\n", M))

poly_basis <- function(x, m) x^m   # scalar or vector

Psi <- matrix(0, N, n_sp * M)
for (i in seq_len(n_sp))
  for (m in seq_len(M))
    Psi[, (i-1)*M + m] <- poly_basis(X_obs[, i], m)

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
  fit_a             <- ADSIHT(Psi_cs, matrix(rhs_c), group, ic.type = "dsic")
  best_a            <- which.min(fit_a$ic)
  ALPHA_ads[j,]     <- fit_a$beta[, best_a] / Psi_sd   # undo standardization
  
  MU_ads[j] <- -U_bar[j] - sum(Psi_bar * ALPHA_ads[j, ])
  
  # grpreg on standardized matrix; coef() returns lambda.min coefficients
  cv            <- cv.grpreg(Psi_cs, rhs_c, group = group,
                             penalty = "grLasso", nfolds = 5)
  ALPHA_gl[j, ] <- coef(cv)[-1] / Psi_sd
  MU_gl[j]      <- -U_bar[j] - sum(Psi_bar * ALPHA_gl[j, ])
}
cat("  Done.\n\n")

# ── Step 4: Jacobian  J_{ji} = ψ'(x_i^wt)ᵀ θ̂_{ji}  ─────────────────────────
# ψ'_m(x) = m · x^{m-1}  (analytic derivative of monomial basis)
dpoly <- function(x, m) m * x^(m-1)

get_jacobian <- function(ALPHA) {
  J <- matrix(0, n_sp, n_sp)
  for (j in seq_len(n_sp))
    for (i in seq_len(n_sp)) {
      cols   <- (i-1)*M + seq_len(M)
      dpsi_i <- sapply(seq_len(M), function(m) dpoly(x_wt[i], m))
      J[j,i] <- sum(ALPHA[j, cols] * dpsi_i)
    }
  J
}

J_ads <- get_jacobian(ALPHA_ads)
J_gl  <- get_jacobian(ALPHA_gl)

get_A   <- function(J) { A <- J; diag(A) <- 0; A }
get_gam <- function(J) -diag(J)

A_ads   <- get_A(J_ads);   g_ads   <- get_gam(J_ads)
A_gl    <- get_A(J_gl);    g_gl    <- get_gam(J_gl)

# ── Step 5: Binary adjacency + metrics ────────────────────────────────────────
# Edge (j←i) exists iff ANY basis coefficient |θ_{ji,m}| > eps (float-safe)
eps_coef <- 1e-10

get_adj <- function(ALPHA) {
  adj <- matrix(FALSE, n_sp, n_sp)
  for (j in seq_len(n_sp))
    for (i in seq_len(n_sp)) {
      cols     <- (i - 1) * M + seq_len(M)
      adj[j, i] <- any(abs(ALPHA[j, cols]) > eps_coef)
    }
  diag(adj) <- FALSE
  adj * 1L
}

adj_true <- (A_true != 0) * 1L
adj_ads  <- get_adj(ALPHA_ads)
adj_gl   <- get_adj(ALPHA_gl)

# True coefficient matrix: θ_true[j, (i-1)M+1] = A_true[j,i]; higher = 0
# (linear GLV: f_{ji}(x) = A_{ji}·x → only first basis coefficient nonzero)
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
  # MCC: Matthews Correlation Coefficient
  denom <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  mcc <- if (denom < 1e-9) 0 else (TP * TN - FP * FN) / denom
  # Coefficient L2: ‖θ̂ - θ_true‖₂ (full pM vector, per-node then averaged)
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
  else paste(
    sprintf("%s<-%s", sp_names[idx[, 1]], sp_names[idx[, 2]]),
    collapse = ", "
  )
  cat(sprintf("  TP: %s\n", fmt_pairs(tp_idx)))
  cat(sprintf("  FP: %s\n", fmt_pairs(fp_idx)))
  cat(sprintf("  FN: %s\n", fmt_pairs(fn_idx)))
}

# mu_true = r_true because f_{ji}(0) = 0 for all i (CLAUDE.md eq.2)
mu_true <- r_true

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
# Reconstructed ODE:  dx_j/dt = μ̂_j + Σ_i ψ(x_i)ᵀ θ̂_{ji}  + u_j
# Test: start from WT steady state, apply 3 held-out perturbations,
# integrate both true and reconstructed systems forward to t=200,
# compare species trajectories.
library(ggplot2)

# Reconstructed ODE using ADSIHT coefficients
ode_recon <- function(t, state, parms) {
  x     <- pmax(state, 0)
  u     <- parms$u
  alpha <- parms$alpha   # n_sp × (n_sp*M) coefficient matrix
  mu    <- parms$mu      # length n_sp intercepts
  dx <- numeric(n_sp)
  for (j in seq_len(n_sp)) {
    fhat_sum <- 0
    for (i in seq_len(n_sp)) {
      cols     <- (i - 1) * M + seq_len(M)
      psi_xi   <- x[i]^seq_len(M)
      fhat_sum <- fhat_sum + sum(alpha[j, cols] * psi_xi)
    }
    dx[j] <- mu[j] + fhat_sum + u[j]
  }
  list(dx)
}

# 3 test perturbations (held-out, not in training set)
u_tests <- list(
  "u: Sp1=-0.5, Sp5=+0.6" = { u <- rep(0, n_sp); u[1] <- -0.5; u[5] <-  0.6; u },
  "u: Sp3=+0.5, Sp7=-0.4" = { u <- rep(0, n_sp); u[3] <-  0.5; u[7] <- -0.4; u },
  "u: all +0.3"            = { rep(0.3, n_sp) }
)

t_span <- seq(0, 200, by = 0.5)

# wide-to-long helper (base R, no dplyr/tidyr)
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
  df_true        <- as.data.frame(out_true)
  colnames(df_true) <- c("time", sp_names)
  
  out_recon <- tryCatch(
    ode(x_wt, t_span, ode_recon,
        list(u = u_vec, alpha = ALPHA_ads, mu = MU_ads),
        method = "lsoda", rtol = 1e-8, atol = 1e-10),
    error = function(e) NULL
  )
  if (is.null(out_recon)) return(NULL)
  df_recon        <- as.data.frame(out_recon)
  colnames(df_recon) <- c("time", sp_names)
  
  rbind(wide_to_long(df_true,  "True",   uname),
        wide_to_long(df_recon, "ADSIHT", uname))
}))

traj_df$species <- factor(traj_df$species, levels = sp_names)
traj_df$model   <- factor(traj_df$model,   levels = c("True", "ADSIHT"))
traj_df$perturb <- factor(traj_df$perturb, levels = names(u_tests))

# Mark WT steady state as horizontal reference
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
    subtitle = "Dotted line = WT steady state  |  starting point = x_wt for all runs",
    x = "Time", y = "Abundance", color = NULL, linetype = NULL
  ) +
  theme_bw(base_size = 8.5) +
  theme(
    strip.background  = element_rect(fill = "grey92"),
    strip.text.x      = element_text(size = 7, face = "bold"),
    strip.text.y      = element_text(size = 7.5, face = "bold"),
    legend.position   = "bottom",
    plot.title        = element_text(size = 10, face = "bold"),
    plot.subtitle     = element_text(size = 7.5, color = "grey45")
  )

print(p_dyn)


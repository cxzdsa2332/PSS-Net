rm(list = ls())

################################################################################
# pss_net_v1.R  —  PSS-Net: 8-species community with nonlinear interactions
#
# Model:    dx_j/dt = μ_j + Σ_i f_{ji}(x_i) + u_j
#           f_{ji}(x) = A_{ji}·x  +  B_{ji}·x²   (mixed linear + quadratic)
# Basis:    ψ(x) = (x, x²)  [M=2, exactly represents true f_{ji}]
# At SS:    Ψ_c · θ_j = −u_{c,j}
#
# Difference from v3: true ODE contains quadratic interaction terms B_{ji}·x²
# for a subset of edges, testing whether ADSIHT/grLasso can recover nonlinear
# edge functions under the monomial basis with M=2.
################################################################################

library(deSolve)
library(ADSIHT)
library(grpreg)

set.seed(42)

# ── Tuning ────────────────────────────────────────────────────────────────────
M        <- 2
n_sp     <- 8
sp_names <- paste0("Sp", seq_len(n_sp))

# ── True parameters ───────────────────────────────────────────────────────────
r_true   <- c(0.8, 1.2, 0.6, 1.0, 0.7, 1.1, 0.5, 0.9)
gam_true <- c(1.5, 1.8, 1.2, 1.6, 1.4, 1.5, 1.0, 1.6)

# Linear interaction coefficients (same network topology as v3)
A_true <- matrix(0, n_sp, n_sp)
A_true[1,3] <-  0.40; A_true[1,5] <- -0.30
A_true[2,1] <-  0.30; A_true[2,4] <- -0.40; A_true[2,7] <-  0.20
A_true[3,2] <- -0.30; A_true[3,6] <-  0.30
A_true[4,1] <-  0.40; A_true[4,3] <- -0.20; A_true[4,8] <-  0.30
A_true[5,2] <-  0.30; A_true[5,6] <- -0.20
A_true[6,4] <-  0.20; A_true[6,5] <-  0.30
A_true[7,3] <- -0.30; A_true[7,8] <-  0.40
A_true[8,1] <-  0.20; A_true[8,6] <- -0.20

# Quadratic interaction coefficients (nonzero for 5 of the 18 edges)
# f_{ji}(x) = A_{ji}·x + B_{ji}·x²;  B≠0 adds curvature to edge function
B_true <- matrix(0, n_sp, n_sp)
B_true[2,1] <-  0.15   # mutualism with saturation
B_true[3,6] <- -0.20   # competition strengthens at high abundance
B_true[4,8] <- -0.15
B_true[6,5] <-  0.18
B_true[7,8] <-  0.20

n_edges    <- sum(A_true != 0)
n_nl_edges <- sum(B_true != 0)   # nonlinear edges

cat("================================================================\n")
cat("  sindy_ss_v1 (PSS-Net + nonlinear): 8-species community\n")
cat(sprintf("  %d species | %d true edges (%d nonlinear) | M=%d\n",
            n_sp, n_edges, n_nl_edges, M))
cat("  Nonlinear edges (f=Ax+Bx²):\n")
for (j in seq_len(n_sp))
  for (i in seq_len(n_sp))
    if (B_true[j, i] != 0)
      cat(sprintf("    %s <- %s :  A=%.2f  B=%.2f\n",
                  sp_names[j], sp_names[i], A_true[j,i], B_true[j,i]))
cat("================================================================\n\n")

# ── ODE (linear + quadratic GLV) ─────────────────────────────────────────────
ode_func <- function(t, state, parms) {
  x <- pmax(state, 0)
  dx <- r_true + as.numeric(A_true %*% x) + as.numeric(B_true %*% x^2) -
        gam_true * x + parms$u
  list(dx)
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

out_wt <- ode(rep(1, n_sp), c(0, 1e4), ode_func, list(u = rep(0, n_sp)),
              method = "lsoda", rtol = 1e-12, atol = 1e-14)
x_wt   <- as.numeric(out_wt[nrow(out_wt), 2:(n_sp + 1)])
cat(sprintf("  Valid conditions: %d / %d\n", N, N_cond))
cat(sprintf("  x_wt: %s\n\n", paste(sprintf("%.3f", x_wt), collapse = ", ")))

# ── Step 2: Polynomial basis  ψ(x) = (x, x²) ────────────────────────────────
cat(sprintf("Step 2: Building polynomial basis (M=%d)...\n", M))

Psi <- matrix(0, N, n_sp * M)
for (i in seq_len(n_sp))
  for (m in seq_len(M))
    Psi[, (i - 1) * M + m] <- X_obs[, i]^m

group    <- rep(seq_len(n_sp), each = M)
Psi_bar  <- colMeans(Psi)
Psi_c    <- sweep(Psi, 2, Psi_bar)
Psi_sd   <- pmax(apply(Psi_c, 2, sd), 1e-10)
Psi_cs   <- sweep(Psi_c, 2, Psi_sd, "/")
U_bar    <- colMeans(U_obs)
U_c      <- sweep(U_obs, 2, U_bar)

cat(sprintf("  Design matrix: %d × %d\n\n", N, ncol(Psi)))

# ── Step 3: Sparse regression ─────────────────────────────────────────────────
cat("Step 3: Sparse regression...\n")

ALPHA_ads <- matrix(0, n_sp, n_sp * M)
ALPHA_gl  <- matrix(0, n_sp, n_sp * M)
MU_ads    <- numeric(n_sp)
MU_gl     <- numeric(n_sp)

for (j in seq_len(n_sp)) {
  rhs_c <- -U_c[, j]

  fit_a         <- ADSIHT(Psi_cs, matrix(rhs_c), group, ic.type = "dsic")
  best_a        <- which.min(fit_a$ic)
  ALPHA_ads[j,] <- fit_a$beta[, best_a] / Psi_sd
  MU_ads[j]     <- -U_bar[j] - sum(Psi_bar * ALPHA_ads[j, ])

  cv            <- cv.grpreg(Psi_cs, rhs_c, group = group,
                             penalty = "grLasso", nfolds = 5)
  ALPHA_gl[j, ] <- coef(cv)[-1] / Psi_sd
  MU_gl[j]      <- -U_bar[j] - sum(Psi_bar * ALPHA_gl[j, ])
}
cat("  Done.\n\n")

# ── Step 4: Jacobian  J_{ji} = ψ'(x_wt[i])ᵀ θ̂_{ji} ─────────────────────────
get_jacobian <- function(alpha_mat) {
  j_mat <- matrix(0, n_sp, n_sp)
  for (j in seq_len(n_sp))
    for (i in seq_len(n_sp)) {
      cols        <- (i - 1) * M + seq_len(M)
      dpsi_i      <- seq_len(M) * x_wt[i]^pmax(seq_len(M) - 1, 0)
      j_mat[j, i] <- sum(alpha_mat[j, cols] * dpsi_i)
    }
  j_mat
}

J_ads <- get_jacobian(ALPHA_ads)
J_gl  <- get_jacobian(ALPHA_gl)

get_A   <- function(J) { A <- J; diag(A) <- 0; A }
get_gam <- function(J) -diag(J)

A_ads <- get_A(J_ads);  g_ads <- get_gam(J_ads)
A_gl  <- get_A(J_gl);   g_gl  <- get_gam(J_gl)

# ── Step 5: Edge detection + metrics ─────────────────────────────────────────
tau_ads <- 1e-10
rel_thr <- 0.01

group_norm_mat <- function(ALPHA) {
  mat <- matrix(0, n_sp, n_sp)
  for (j in seq_len(n_sp))
    for (i in seq_len(n_sp)) {
      cols      <- (i - 1) * M + seq_len(M)
      mat[j, i] <- sqrt(sum(ALPHA[j, cols]^2))
    }
  mat
}

norm_ads   <- group_norm_mat(ALPHA_ads)
norm_gl    <- group_norm_mat(ALPHA_gl)
tau_gl_vec <- apply(norm_gl, 1, max) * rel_thr

get_adj_ads <- function(nm, tau) { adj <- nm >= tau; diag(adj) <- FALSE; adj * 1L }
get_adj_gl  <- function(nm, tv)  {
  adj <- matrix(FALSE, n_sp, n_sp)
  for (j in seq_len(n_sp)) adj[j, ] <- nm[j, ] >= tv[j]
  diag(adj) <- FALSE; adj * 1L
}

adj_true <- (A_true != 0) * 1L
adj_ads  <- get_adj_ads(norm_ads, tau_ads)
adj_gl   <- get_adj_gl(norm_gl,  tau_gl_vec)

# True Jacobian: J_{ji} = A_{ji} + 2·B_{ji}·x_wt[i]  (at x_wt)
J_true    <- A_true + 2 * B_true * matrix(x_wt, n_sp, n_sp, byrow = TRUE)
A_jac_true <- J_true; diag(A_jac_true) <- 0

# True alpha: theta_{ji,1} = A_{ji},  theta_{ji,2} = B_{ji}
alpha_true <- matrix(0, n_sp, n_sp * M)
for (j in seq_len(n_sp))
  for (i in seq_len(n_sp)) {
    alpha_true[j, (i-1)*M + 1] <- A_true[j, i]
    alpha_true[j, (i-1)*M + 2] <- B_true[j, i]
  }

mets <- function(adj_h, ALPHA_h) {
  TP <- sum(adj_true == 1 & adj_h == 1)
  FP <- sum(adj_true == 0 & adj_h == 1)
  FN <- sum(adj_true == 1 & adj_h == 0)
  TN <- sum(adj_true == 0 & adj_h == 0)
  pr    <- TP / (TP + FP + 1e-9)
  re    <- TP / (TP + FN + 1e-9)
  f1    <- 2 * pr * re / (pr + re + 1e-9)
  denom <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  mcc   <- if (denom < 1e-9) 0 else (TP * TN - FP * FN) / denom
  coef_l2  <- mean(sqrt(rowSums((ALPHA_h - alpha_true)^2)))
  jac_rmse <- sqrt(mean((get_A(get_jacobian(ALPHA_h)) - A_jac_true)^2))
  list(TP=TP, FP=FP, FN=FN, TN=TN, Pr=pr, Re=re, F1=f1, mcc=mcc,
       coef_l2=coef_l2, jac_rmse=jac_rmse)
}

m_ads <- mets(adj_ads, ALPHA_ads)
m_gl  <- mets(adj_gl,  ALPHA_gl)

# ── Print adjacency detail ────────────────────────────────────────────────────
print_adj <- function(label, adj_h) {
  cat(sprintf("\n%s adjacency (row=target j, col=source i):\n", label))
  mat <- matrix(sprintf("%d", adj_h), n_sp, n_sp,
                dimnames = list(sp_names, sp_names))
  print(noquote(mat))
  tp_idx <- which(adj_true == 1 & adj_h == 1, arr.ind = TRUE)
  fp_idx <- which(adj_true == 0 & adj_h == 1, arr.ind = TRUE)
  fn_idx <- which(adj_true == 1 & adj_h == 0, arr.ind = TRUE)
  fmt <- function(idx)
    if (nrow(idx) == 0) "none"
    else paste(sprintf("%s<-%s", sp_names[idx[,1]], sp_names[idx[,2]]), collapse=", ")
  cat(sprintf("  TP: %s\n  FP: %s\n  FN: %s\n", fmt(tp_idx), fmt(fp_idx), fmt(fn_idx)))
}

# ── Parameter recovery: r and gamma ──────────────────────────────────────────
mu_true <- r_true

cat("── r recovery ───────────────────────────────────────────────────\n")
cat(sprintf("  %-6s  %7s  %9s  %9s\n", "Sp","r_true","r_ADSIHT","r_grLasso"))
for (j in seq_len(n_sp))
  cat(sprintf("  %-6s  %7.4f  %9.4f  %9.4f\n",
              sp_names[j], r_true[j], MU_ads[j], MU_gl[j]))
cat(sprintf("  RMSE   ADSIHT=%.4f  GroupLasso=%.4f\n",
            sqrt(mean((MU_ads-mu_true)^2)), sqrt(mean((MU_gl-mu_true)^2))))

cat("\n── gamma recovery ───────────────────────────────────────────────\n")
cat(sprintf("  %-6s  %9s  %9s  %9s\n","Sp","gam_true","g_ADSIHT","g_grLasso"))
for (j in seq_len(n_sp))
  cat(sprintf("  %-6s  %9.4f  %9.4f  %9.4f\n",
              sp_names[j], gam_true[j], g_ads[j], g_gl[j]))
cat(sprintf("  RMSE   ADSIHT=%.4f  GroupLasso=%.4f\n",
            sqrt(mean((g_ads-gam_true)^2)), sqrt(mean((g_gl-gam_true)^2))))

# ── Nonlinear coefficient recovery ────────────────────────────────────────────
cat("\n── Nonlinear coefficient recovery (B_{ji}, m=2 basis) ──────────\n")
cat(sprintf("  %-18s  %7s  %9s  %9s\n","Edge","B_true","B_ADSIHT","B_grLasso"))
for (j in seq_len(n_sp))
  for (i in seq_len(n_sp))
    if (A_true[j,i] != 0) {
      col2 <- (i-1)*M + 2
      cat(sprintf("  %s <- %-6s  %7.4f  %9.4f  %9.4f\n",
                  sp_names[j], sp_names[i],
                  B_true[j,i], ALPHA_ads[j,col2], ALPHA_gl[j,col2]))
    }

print_adj("ADSIHT",      adj_ads)
print_adj("Group-Lasso", adj_gl)

cat("\n================================================================\n")
cat(sprintf("  v1 (nonlinear)  M=%d | N=%d | %d edges (%d nonlinear)\n",
            M, N, n_edges, n_nl_edges))
fmt_row <- function(label, m)
  cat(sprintf(
    "  %-12s  TP=%2d FP=%2d FN=%2d TN=%2d  Pr=%.3f Re=%.3f F1=%.3f MCC=%.3f  CoefL2=%.4f JacRMSE=%.4f\n",
    label, m$TP, m$FP, m$FN, m$TN, m$Pr, m$Re, m$F1, m$mcc, m$coef_l2, m$jac_rmse))
fmt_row("ADSIHT",      m_ads)
fmt_row("Group-Lasso", m_gl)
cat("================================================================\n")

# ── Effect decomposition: true f_{ji}(x) vs inferred, at x_wt ────────────────
library(ggplot2)
x_range  <- seq(0.2, max(x_wt) * 1.6, length.out = 200)
nl_pairs <- which(B_true != 0, arr.ind = TRUE)  # nonlinear edges only

curve_rows <- vector("list", nrow(nl_pairs) * 2)
idx <- 1L
for (k in seq_len(nrow(nl_pairs))) {
  j <- nl_pairs[k, 1]; i <- nl_pairs[k, 2]
  cols   <- (i-1)*M + seq_len(M)
  f_true <- A_true[j,i]*x_range + B_true[j,i]*x_range^2
  f_ads  <- ALPHA_ads[j, cols[1]]*x_range + ALPHA_ads[j, cols[2]]*x_range^2
  label  <- sprintf("%s <- %s  (A=%.2f, B=%.2f)", sp_names[j], sp_names[i],
                    A_true[j,i], B_true[j,i])
  curve_rows[[idx]]   <- data.frame(x=x_range, f=f_true, model="True",   edge=label)
  curve_rows[[idx+1]] <- data.frame(x=x_range, f=f_ads,  model="ADSIHT", edge=label)
  idx <- idx + 2L
}
curve_df       <- do.call(rbind, curve_rows)
curve_df$model <- factor(curve_df$model, levels = c("True","ADSIHT"))

xwt_df <- data.frame(
  edge = sapply(seq_len(nrow(nl_pairs)), function(k)
    sprintf("%s <- %s  (A=%.2f, B=%.2f)",
            sp_names[nl_pairs[k,1]], sp_names[nl_pairs[k,2]],
            A_true[nl_pairs[k,1], nl_pairs[k,2]],
            B_true[nl_pairs[k,1], nl_pairs[k,2]])),
  xwt = x_wt[nl_pairs[, 2]]
)

p_nl <- ggplot(curve_df, aes(x=x, y=f, color=model, linetype=model)) +
  geom_hline(yintercept=0, color="grey80", linewidth=0.3) +
  geom_vline(data=xwt_df, aes(xintercept=xwt),
             color="grey60", linewidth=0.4, linetype="dotted") +
  geom_line(linewidth=0.9) +
  facet_wrap(~ edge, scales="free_y", ncol=2) +
  scale_color_manual(values=c(True="grey20", ADSIHT="tomato3")) +
  scale_linetype_manual(values=c(True="solid", ADSIHT="dashed")) +
  labs(
    title    = sprintf("Nonlinear edge function recovery  [M=%d, %d NL edges]", M, n_nl_edges),
    subtitle = "Dotted line = x_wt of source species",
    x = "Source abundance x_i", y = expression(f[ji](x[i])),
    color=NULL, linetype=NULL
  ) +
  theme_bw(base_size=9) +
  theme(strip.background=element_rect(fill="grey92"),
        strip.text=element_text(size=7.5, face="bold"),
        legend.position="bottom",
        plot.title=element_text(size=10, face="bold"),
        plot.subtitle=element_text(size=8, color="grey45"))

print(p_nl)

# ── Step 6: Dynamics comparison — true ODE vs ADSIHT-reconstructed ODE ────────
# Reconstructed ODE uses inferred θ̂ with monomial basis evaluation
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
  "u: all +0.3"            = rep(0.3, n_sp)
)

t_span <- seq(0, 200, by = 0.5)

wide_to_long <- function(df_wide, model_label, perturb_label) {
  do.call(rbind, lapply(sp_names, function(sp) {
    data.frame(time=df_wide$time, species=sp, x=df_wide[[sp]],
               model=model_label, perturb=perturb_label,
               stringsAsFactors=FALSE)
  }))
}

traj_df <- do.call(rbind, lapply(names(u_tests), function(uname) {
  u_vec <- u_tests[[uname]]

  out_true <- ode(x_wt, t_span, ode_func, list(u=u_vec),
                  method="lsoda", rtol=1e-8, atol=1e-10)
  df_true           <- as.data.frame(out_true)
  colnames(df_true) <- c("time", sp_names)

  out_recon <- tryCatch(
    ode(x_wt, t_span, ode_recon,
        list(u=u_vec, alpha=ALPHA_ads, mu=MU_ads),
        method="lsoda", rtol=1e-8, atol=1e-10),
    error = function(e) NULL)
  if (is.null(out_recon)) return(NULL)
  df_recon           <- as.data.frame(out_recon)
  colnames(df_recon) <- c("time", sp_names)

  rbind(wide_to_long(df_true,  "True",   uname),
        wide_to_long(df_recon, "ADSIHT", uname))
}))

traj_df$species <- factor(traj_df$species, levels=sp_names)
traj_df$model   <- factor(traj_df$model,   levels=c("True","ADSIHT"))
traj_df$perturb <- factor(traj_df$perturb, levels=names(u_tests))

wt_df <- data.frame(species=factor(sp_names, levels=sp_names), xwt=x_wt)

p_dyn <- ggplot(traj_df, aes(x=time, y=x, color=model, linetype=model)) +
  geom_hline(data=wt_df, aes(yintercept=xwt),
             color="grey80", linewidth=0.35, linetype="dotted") +
  geom_line(linewidth=0.75) +
  facet_grid(species ~ perturb, scales="free_y") +
  scale_color_manual(values=c(True="grey20", ADSIHT="tomato3")) +
  scale_linetype_manual(values=c(True="solid", ADSIHT="dashed")) +
  labs(
    title    = sprintf("Dynamics: true ODE vs ADSIHT  [M=%d, %d NL edges]", M, n_nl_edges),
    subtitle = "Dotted = WT steady state  |  starting from x_wt",
    x="Time", y="Abundance", color=NULL, linetype=NULL
  ) +
  theme_bw(base_size=8.5) +
  theme(strip.background=element_rect(fill="grey92"),
        strip.text.x=element_text(size=7, face="bold"),
        strip.text.y=element_text(size=7.5, face="bold"),
        legend.position="bottom",
        plot.title=element_text(size=10, face="bold"),
        plot.subtitle=element_text(size=7.5, color="grey45"))

print(p_dyn)

# ── Step 7: Network comparison — true vs ADSIHT (igraph) ──────────────────────
library(igraph)

make_graph_from_adj <- function(adj_mat, weight_mat) {
  el <- which(adj_mat != 0, arr.ind=TRUE)   # [row=j, col=i]
  if (nrow(el) == 0) return(make_empty_graph(n_sp, directed=TRUE))
  g <- graph_from_edgelist(cbind(el[,2], el[,1]), directed=TRUE)
  E(g)$weight <- weight_mat[el]
  V(g)$name   <- sp_names
  g
}

# Use Jacobian-derived A for edge sign (captures nonlinear slope at x_wt)
g_true <- make_graph_from_adj(adj_true,  A_jac_true)
g_ads  <- make_graph_from_adj(adj_ads,   A_ads)

lay <- layout_in_circle(g_true)

edge_col <- function(w) ifelse(w > 0, "tomato3", "steelblue3")

tp_mat <- adj_true == 1 & adj_ads == 1
fp_mat <- adj_true == 0 & adj_ads == 1
fn_mat <- adj_true == 1 & adj_ads == 0

par(mfrow=c(1,2), mar=c(1,1,2.5,1), bg="white")

# Panel 1: True network
plot(g_true,
     layout=lay, vertex.color="grey88", vertex.frame.color="grey55",
     vertex.size=26, vertex.label=sp_names, vertex.label.cex=0.85,
     vertex.label.color="grey10",
     edge.color=edge_col(E(g_true)$weight),
     edge.width=2, edge.arrow.size=0.55, edge.curved=0.25,
     main=sprintf("True network  (%d NL edges marked *)", n_nl_edges))

# Mark nonlinear edge sources with a star annotation
nl_nodes <- unique(nl_pairs[,2])
for (ni in nl_nodes) {
  coords <- lay[ni, ]
  text(coords[1]*1.22, coords[2]*1.22, "*", cex=1.3, col="darkorange")
}
legend("bottomleft", bty="n", cex=0.72,
       legend=c("Promotion (+)","Inhibition (−)","* NL source"),
       col=c("tomato3","steelblue3","darkorange"), lwd=c(2,2,NA), pch=c(NA,NA,8))

# Panel 2: ADSIHT inferred
ads_el    <- ends(g_ads, E(g_ads), names=FALSE)
ads_class <- ifelse(
  mapply(function(src,tgt) tp_mat[tgt,src], ads_el[,1], ads_el[,2]),
  "TP", "FP")
ads_w    <- E(g_ads)$weight
ads_ecol <- ifelse(ads_class=="TP", edge_col(ads_w), "orange2")

plot(g_ads,
     layout=lay, vertex.color="grey88", vertex.frame.color="grey55",
     vertex.size=26, vertex.label=sp_names, vertex.label.cex=0.85,
     vertex.label.color="grey10",
     edge.color=ads_ecol,
     edge.width=ifelse(ads_class=="FP", 1.2, 2),
     edge.lty=ifelse(ads_class=="FP", 2, 1),
     edge.arrow.size=0.55, edge.curved=0.25,
     main=sprintf("ADSIHT  [TP=%d FP=%d FN=%d  MCC=%.2f]",
                  m_ads$TP, m_ads$FP, m_ads$FN, m_ads$mcc))

fn_el <- which(fn_mat, arr.ind=TRUE)
if (nrow(fn_el) > 0) {
  g_fn <- make_empty_graph(n_sp, directed=TRUE)
  g_fn <- add_edges(g_fn, as.vector(t(cbind(fn_el[,2], fn_el[,1]))))
  plot(g_fn, layout=lay, add=TRUE,
       vertex.color=NA, vertex.frame.color=NA, vertex.label=NA,
       vertex.size=26, edge.color="grey60", edge.width=1.2,
       edge.lty=3, edge.arrow.size=0.45, edge.curved=0.25)
}

legend("bottomleft", bty="n", cex=0.72,
       legend=c("TP promotion","TP inhibition","FP","FN (missed)"),
       col=c("tomato3","steelblue3","orange2","grey60"),
       lwd=c(2,2,1.2,1.2), lty=c(1,1,2,3))

par(mfrow=c(1,1))

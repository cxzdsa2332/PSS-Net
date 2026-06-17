rm(list = ls())

################################################################################
# pss_net_glv_ss.R  —  PSS-Net: 真乘性广义 Lotka-Volterra (gLV) 稳态回归验证
#
# 目的：验证 §1.1 的稳态约化——真乘性 gLV
#         dx_j/dt = x_j ( r_j + Σ_{i≠j} A_ji x_i - γ_j x_j + u_j )
#       在可行正稳态处除以 x_j 得加性线性关系 (2'):
#         r_j + Σ_{i≠j} A_ji x_i* - γ_j x_j* + u_j = 0
#       因此对真 gLV 积分得到的稳态，PSS 加性稳态回归（单项式 M=2 + ADSIHT）
#       应当（a）与"线性求解 (diag γ - A) x* = r+u"得到的稳态完全一致；
#            （b）成功恢复互作网络 A 与自调节 γ。
#
# 对照：同一 (A, r, γ, u) 下，加性线性 ODE 与乘性 gLV 稳态值是否相同。
#
# Input:   none（自生成）
# Output:  results/sim_results/glv_ss_verification.csv  — 每 seed 的恢复指标与稳态一致性
################################################################################

suppressMessages({
  library(ADSIHT)
  library(deSolve)
})

set.seed(1)
M_ord <- 2L

## ---------------------------------------------------------------- 真值系统 ----
# 对角占优保证可行正稳态：γ_j > Σ_i|A_ji|
make_system <- function(p = 8, n_in = 2, seed = 1) {
  set.seed(seed)
  A <- matrix(0, p, p)                         # A[j,i]: i->j，对角为 0
  for (j in seq_len(p)) {
    src <- sample(setdiff(seq_len(p), j), n_in)
    A[j, src] <- runif(n_in, 0.2, 0.5) * sample(c(-1, 1), n_in, replace = TRUE)
  }
  gamma <- rowSums(abs(A)) + runif(p, 1.0, 1.5)
  r <- runif(p, 0.8, 1.6)
  list(p = p, A = A, gamma = gamma, r = r, adj = (A != 0) * 1)
}

# 线性求解稳态： (diag(γ) - A) x* = r + u   <=>  (2')
ss_linear <- function(sys, u) as.numeric(solve(diag(sys$gamma) - sys$A, sys$r + u))

# 真乘性 gLV 积分至稳态： dx_j = x_j ( r_j + Σ A_ji x_i - γ_j x_j + u_j )
ss_glv <- function(sys, u, t_max = 3000) {
  deriv <- function(t, x, parms) {
    rate <- sys$r + as.numeric(sys$A %*% x) - sys$gamma * x + u
    list(x * rate)
  }
  x0 <- pmax(ss_linear(sys, u), 0.05)          # 以线性稳态为初值
  out <- tryCatch(
    ode(y = x0, times = c(0, t_max), func = deriv, parms = NULL,
        method = "lsoda", rtol = 1e-10, atol = 1e-12),
    error = function(e) NULL)
  if (is.null(out)) return(rep(NA, sys$p))
  as.numeric(out[2, -1])
}

## ------------------------------------------------------ 基函数 / 推断 ----
psi_row <- function(xvec) as.vector(sapply(xvec, function(x) x^(seq_len(M_ord))))

infer_network <- function(sys, U, X) {
  p <- sys$p
  Psi <- t(apply(X, 1, psi_row))
  group <- rep(seq_len(p), each = M_ord)
  Psi_c <- sweep(Psi, 2, colMeans(Psi))
  sdv <- apply(Psi_c, 2, sd); sdv[sdv < 1e-10] <- 1e-10
  Psi_cs <- sweep(Psi_c, 2, sdv, FUN = "/")
  adj_est <- matrix(0, p, p)
  theta_hat <- matrix(0, p * M_ord, p)
  for (j in seq_len(p)) {
    rhs <- -(U[, j] - mean(U[, j]))
    fit <- tryCatch(ADSIHT(Psi_cs, matrix(rhs), group, ic.type = "dsic"),
                    error = function(e) NULL)
    if (is.null(fit)) next
    th <- fit$beta[, which.min(fit$ic)] / sdv
    theta_hat[, j] <- th
    gnorm <- sapply(seq_len(p), function(i)
      sqrt(sum(th[((i - 1) * M_ord + 1):(i * M_ord)]^2)))
    adj_est[j, gnorm >= 1e-8] <- 1
  }
  diag(adj_est) <- 0
  list(adj_est = adj_est, theta_hat = theta_hat)
}

edge_metrics <- function(est, true) {
  off <- which(row(true) != col(true))
  e <- est[off]; t <- true[off]
  TP <- sum(e == 1 & t == 1); FP <- sum(e == 1 & t == 0)
  TN <- sum(e == 0 & t == 0); FN <- sum(e == 0 & t == 1)
  pr <- ifelse(TP + FP == 0, 0, TP / (TP + FP))
  re <- ifelse(TP + FN == 0, 0, TP / (TP + FN))
  den <- sqrt(as.numeric(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  mcc <- ifelse(den == 0, 0, (TP * TN - FP * FN) / den)
  c(Pr = pr, Re = re, MCC = mcc)
}

# 系数真值（稳态关系 2'）：i≠j 为互作 A_ji；i=j 为自调节 -γ_j；均为一阶
coef_l2 <- function(sys, theta_hat) {
  p <- sys$p; err <- 0
  for (j in seq_len(p)) {
    th_true <- numeric(p * M_ord)
    for (i in seq_len(p)) {
      th_true[(i - 1) * M_ord + 1] <- if (i == j) -sys$gamma[j] else sys$A[j, i]
    }
    err <- err + sqrt(sum((theta_hat[, j] - th_true)^2))
  }
  err / p
}

## ----------------------------------------------------------- 主循环 ----
R <- 10
N <- 200
sigma <- 0.04
rows <- list()
for (seed in seq_len(R)) {
  sys <- make_system(p = 8, n_in = 2, seed = 300 + seed)
  U <- matrix(runif(N * sys$p, -0.3, 0.5), N, sys$p)
  U[1, ] <- 0
  # 两种真值稳态
  X_glv <- t(apply(U, 1, function(u) ss_glv(sys, u)))
  X_lin <- t(apply(U, 1, function(u) ss_linear(sys, u)))
  ok <- apply(X_glv, 1, function(r) all(is.finite(r)) && all(r > 0))
  # 稳态一致性：gLV 积分 vs 线性求解
  ss_maxdiff <- max(abs(X_glv[ok, ] - X_lin[ok, ]))
  # 含噪观测 + 推断（用真 gLV 稳态）
  X_obs <- X_glv + matrix(rnorm(length(X_glv), sd = sigma), nrow(X_glv))
  res <- infer_network(sys, U[ok, , drop = FALSE], X_obs[ok, , drop = FALSE])
  m <- edge_metrics(res$adj_est, sys$adj)
  rows[[seed]] <- data.frame(
    seed = seed, N_eff = sum(ok), ss_maxdiff = ss_maxdiff,
    Pr = m["Pr"], Re = m["Re"], MCC = m["MCC"],
    CoefL2 = coef_l2(sys, res$theta_hat))
  cat(sprintf("seed %2d | N_eff=%3d | ss_maxdiff=%.2e | MCC=%.3f\n",
              seed, sum(ok), ss_maxdiff, m["MCC"]))
}
df <- do.call(rbind, rows); rownames(df) <- NULL

dir.create("results/sim_results", showWarnings = FALSE, recursive = TRUE)
write.csv(df, "results/sim_results/glv_ss_verification.csv", row.names = FALSE)

cat("\n===== 真乘性 gLV 稳态回归验证（", R, "seeds，N=", N, "）=====\n", sep = "")
cat(sprintf("稳态一致性 max|x_glv - x_linear|  : mean %.2e  (max %.2e)\n",
            mean(df$ss_maxdiff), max(df$ss_maxdiff)))
cat(sprintf("网络恢复  MCC                     : %.3f ± %.3f\n",
            mean(df$MCC), sd(df$MCC)))
cat(sprintf("          Precision / Recall      : %.3f / %.3f\n",
            mean(df$Pr), mean(df$Re)))
cat(sprintf("系数误差  CoefL2 (vs 真 A)        : %.3f ± %.3f\n",
            mean(df$CoefL2), sd(df$CoefL2)))
cat("\nSaved: results/sim_results/glv_ss_verification.csv\n")

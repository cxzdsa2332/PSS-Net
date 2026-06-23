rm(list = ls())

################################################################################
# pss_net_design_nl_seq.R  —  经过设计的扰动与随机扰动比较（强非线性版）
#
# 在 pss_net_design_nl.R 基础上加大二次项强度 B_ji（更强非线性），
# 对比三种扰动设计在固定预算 N 下的网络恢复质量：
#   random / maximin(u 空间填充) / dopt(AlgDesign exact D-optimal)
#
# 回归:  逐节点 ADSIHT（中心化+标准化），单项式基 M=2（含 x, x^2）
# Input:   none
# Output:  results/sim_results/design_nl_seq_comparison.csv
# Plot:    analysis_script/plot_design_curves.R
################################################################################

suppressMessages({
  library(ADSIHT)
  library(AlgDesign)
  library(deSolve)
})

set.seed(1)
M_ord <- 2L

## ---------------------------------------------------------------- 真值系统 ----
# f_ji(x)=A_ji*x + Bq_ji*x^2 (i!=j); 自调节 f_jj=-gamma_j*x_j
make_system_nl <- function(p = 8, n_in = 2, seed = 1) {
  set.seed(seed)
  A <- matrix(0, p, p)
  Bq <- matrix(0, p, p)
  for (j in seq_len(p)) {
    src <- sample(setdiff(seq_len(p), j), n_in)
    A[j, src] <- runif(n_in, 0.3, 0.7) * sample(c(-1, 1), n_in, replace = TRUE)
    # 一半边附加【强】二次项（强非线性互作）
    nl <- src[seq_len(max(1, floor(n_in / 2)))]
    Bq[j, nl] <- runif(length(nl), 0.40, 0.80) * sample(c(-1, 1), length(nl), TRUE)
  }
  gamma <- runif(p, 3.0, 4.0)                  # 强自调节，保证非线性稳态稳定为正
  mu <- runif(p, 1.0, 2.0)
  sys <- list(p = p, A = A, Bq = Bq, gamma = gamma, mu = mu,
              adj = ((A != 0) | (Bq != 0)) * 1)
  sys$x_wt <- steady_one(sys, rep(0, p), x0 = mu / gamma)
  sys
}

# 单条件稳态：从 x0 积分 ODE 至 t_max 取末值
steady_one <- function(sys, u, x0 = NULL, t_max = 2000) {
  p <- sys$p
  if (is.null(x0)) x0 <- sys$x_wt
  deriv <- function(t, x, parms) {
    inter <- numeric(p)
    for (j in seq_len(p)) {
      others <- setdiff(seq_len(p), j)
      inter[j] <- sum(sys$A[j, others] * x[others] +
                      sys$Bq[j, others] * x[others]^2)
    }
    list(sys$mu - sys$gamma * x + inter + u)
  }
  out <- tryCatch(
    ode(y = x0, times = c(0, t_max), func = deriv, parms = NULL,
        method = "lsoda", rtol = 1e-9, atol = 1e-11),
    error = function(e) NULL)
  if (is.null(out)) return(rep(NA, p))
  as.numeric(out[2, -1])
}

steady_many <- function(sys, U) {
  t(apply(U, 1, function(u) steady_one(sys, u)))
}

# 线性化 Jacobian 矩阵（at wt）：B_jac[j,i]=A+2Bq*x_wt (i!=j), 对角 -gamma
jac_lin <- function(sys) {
  p <- sys$p
  Bj <- matrix(0, p, p)
  for (j in seq_len(p)) for (i in seq_len(p)) {
    if (i == j) Bj[j, i] <- -sys$gamma[j]
    else Bj[j, i] <- sys$A[j, i] + 2 * sys$Bq[j, i] * sys$x_wt[i]
  }
  Bj
}

## ------------------------------------------------------- 基函数 ----
psi_row <- function(xvec) as.vector(sapply(xvec, function(x) x^(seq_len(M_ord))))
design_row <- function(xvec) psi_row(xvec)

make_pool <- function(p, n_pool = 2500) {
  matrix(runif(n_pool * p, -0.4, 0.8), n_pool, p)
}

## ------------------------------------------------- 三策略：返回 N_max 有序索引 ----
order_random <- function(sys, pool, N_max) seq_len(N_max)

order_maximin <- function(sys, pool, N_max) {
  idx <- 1L
  mind <- as.vector(sqrt(rowSums((sweep(pool, 2, pool[1, ]))^2)))
  for (k in 2:N_max) {
    nxt <- which.max(mind); idx <- c(idx, nxt)
    mind <- pmin(mind, sqrt(rowSums((sweep(pool, 2, pool[nxt, ]))^2)))
    mind[idx] <- -Inf
  }
  idx
}

prepare_design_features <- function(Phi) {
  Phi <- sweep(Phi, 2, colMeans(Phi))
  s <- apply(Phi, 2, sd)
  s[!is.finite(s) | s < 1e-10] <- 1
  Phi <- sweep(Phi, 2, s, "/")
  q <- qr(Phi, tol = 1e-9)
  if (q$rank < ncol(Phi)) Phi <- Phi[, sort(q$pivot[seq_len(q$rank)]), drop = FALSE]
  Phi
}

select_dopt <- function(sys, pool, N, seed) {
  Bj <- jac_lin(sys)
  X_lin <- t(sys$x_wt - solve(Bj, t(pool)))
  Phi <- prepare_design_features(t(apply(X_lin, 1, design_row)))
  if (N < ncol(Phi) + 1L) return(integer(0))  # 加截距后的 exact D-opt 模型秩
  set.seed(seed)
  fit <- AlgDesign::optFederov(
    frml = ~ ., data = as.data.frame(Phi), nTrials = N, criterion = "D",
    augment = TRUE, rows = 1L, maxIteration = 100, nRepeats = 1
  )
  as.integer(fit$rows)
}

## ------------------------------------------------------ 逐节点 ADSIHT 推断 ----
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

## ---------------------------------------------------------------- 指标 ----
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

# Jacobian 真值 J_ji = A + 2 Bq x_wt
jac_rmse <- function(sys, theta_hat) {
  p <- sys$p; Jhat <- matrix(0, p, p)
  for (j in seq_len(p)) for (i in seq_len(p)) {
    th <- theta_hat[((i - 1) * M_ord + 1):(i * M_ord), j]
    Jhat[j, i] <- sum(seq_len(M_ord) * sys$x_wt[i]^(seq_len(M_ord) - 1) * th)
  }
  Jtrue <- sys$A + 2 * sys$Bq * matrix(sys$x_wt, p, p, byrow = TRUE)
  off <- which(row(sys$A) != col(sys$A))
  sqrt(mean((Jhat[off] - Jtrue[off])^2))
}

# 系数真值 theta_ji = (A_ji, Bq_ji)
coef_l2 <- function(sys, theta_hat) {
  p <- sys$p; err <- 0
  for (j in seq_len(p)) {
    th_true <- numeric(p * M_ord)
    for (i in seq_len(p)) {
      th_true[(i - 1) * M_ord + 1] <- sys$A[j, i]
      th_true[(i - 1) * M_ord + 2] <- sys$Bq[j, i]
    }
    err <- err + sqrt(sum((theta_hat[, j] - th_true)^2))
  }
  err / p
}

## ----------------------------------------------------------- 主实验循环 ----
N_grid <- c(12, 16, 20, 30, 40, 60)
N_max <- max(N_grid)
strategies <- c("random", "maximin", "dopt")
R <- 20
sigma <- 0.04

rows <- list()
for (seed in seq_len(R)) {
  sys <- make_system_nl(p = 8, n_in = 2, seed = 200 + seed)
  if (any(!is.finite(sys$x_wt)) || any(sys$x_wt <= 0)) { cat("seed", seed, "skip (bad wt)\n"); next }
  pool <- make_pool(sys$p, n_pool = 2500)
  pool[1, ] <- 0                                # 三种策略共享 WT 条件
  index_by_strategy <- list(
    random = setNames(lapply(N_grid, function(N) order_random(sys, pool, N)), N_grid),
    maximin = {
      ord <- order_maximin(sys, pool, N_max)
      setNames(lapply(N_grid, function(N) ord[seq_len(N)]), N_grid)
    },
    dopt = setNames(lapply(N_grid, function(N)
      select_dopt(sys, pool, N, 300000L + 100L * seed + N)), N_grid)
  )
  union_idx <- unique(unlist(index_by_strategy, use.names = FALSE))
  X_union <- steady_many(sys, pool[union_idx, , drop = FALSE])
  E_union <- matrix(rnorm(length(X_union), sd = sigma), nrow(X_union))
  X_union_obs <- X_union + E_union
  for (strat in strategies) {
    for (N in N_grid) {
      idx <- index_by_strategy[[strat]][[as.character(N)]]
      if (length(idx) == 0L) next
      X_full <- X_union_obs[match(idx, union_idx), , drop = FALSE]
      ok <- apply(X_full, 1, function(r) all(is.finite(r)) && all(r > 0))
      use <- which(ok)
      if (length(use) < 8) next
      res <- infer_network(sys, pool[idx[use], , drop = FALSE], X_full[use, , drop = FALSE])
      m <- edge_metrics(res$adj_est, sys$adj)
      rows[[length(rows) + 1]] <- data.frame(
        seed = seed, N = N, n_eff = length(use), strategy = strat,
        Pr = m["Pr"], Re = m["Re"], MCC = m["MCC"],
        CoefL2 = coef_l2(sys, res$theta_hat),
        JacRMSE = jac_rmse(sys, res$theta_hat))
    }
  }
  cat("seed", seed, "done\n")
}
df <- do.call(rbind, rows); rownames(df) <- NULL

dir.create("results/sim_results", showWarnings = FALSE, recursive = TRUE)
write.csv(df, "results/sim_results/design_nl_seq_comparison.csv", row.names = FALSE)

agg <- aggregate(cbind(MCC, Pr, Re, CoefL2, JacRMSE) ~ strategy + N, df, mean)
agg <- agg[order(agg$N, agg$strategy), ]
cat("\n===== strong nonlinear, mean over seeds =====\n")
print(agg, row.names = FALSE)
cat("\nSaved: results/sim_results/design_nl_seq_comparison.csv\n")

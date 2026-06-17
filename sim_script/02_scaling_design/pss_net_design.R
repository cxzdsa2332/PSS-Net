rm(list = ls())

################################################################################
# pss_net_design.R  —  PSS-Net 创新点 A：最优扰动设计（Optimal Perturbation Design）
#
# 比较三种扰动设计策略在固定预算 N 下的网络恢复质量：
#   (1) random   —  Uniform 随机扰动（现状基线）
#   (2) maximin  —  u 空间贪心空间填充（Latin-hypercube 风格）
#   (3) dopt     —  序贯 D-最优主动设计（本方法，eq.3 + Sherman-Morrison）
#
# 系统:  线性 GLV，稳态闭式  x*(u) = x_wt - B^{-1} u   (methods/optimal_perturbation_design.md eq.1)
# 回归:  逐节点 ADSIHT（中心化 + 标准化），单项式基 M=2
# 指标:  MCC / Precision / Recall / CoefL2 / JacRMSE，扫描 N，重复 R 个种子
#
# Input:   none（自生成模拟数据）
# Output:  results/sim_results/design_comparison.csv  — 每 (strategy,N,seed) 指标
# Plot:    analysis_script/plot_design_curves.R
################################################################################

suppressMessages({
  library(ADSIHT)
  library(MASS)
})

set.seed(1)

## ---------------------------------------------------------------- 真值系统 ----
# 生成稀疏线性 GLV：B = A_offdiag - diag(gamma)，Hurwitz 且 x_wt > 0
make_system <- function(p = 8, n_in = 2, seed = 1) {
  set.seed(seed)
  A <- matrix(0, p, p)                       # A[j,i]: i -> j 的互作强度
  for (j in seq_len(p)) {
    src <- sample(setdiff(seq_len(p), j), n_in)
    A[j, src] <- runif(n_in, 0.3, 0.8) * sample(c(-1, 1), n_in, replace = TRUE)
  }
  gamma <- runif(p, 2.0, 3.0)                 # 强自调节，保证稳定与正稳态
  B <- A - diag(gamma)
  # 保证 Hurwitz：实部全负；必要时加大自调节
  while (max(Re(eigen(B, only.values = TRUE)$values)) > -0.2) {
    gamma <- gamma + 0.5
    B <- A - diag(gamma)
  }
  mu <- runif(p, 0.8, 1.6)                    # 内禀增长，使 x_wt > 0
  x_wt <- as.vector(-solve(B, mu))
  list(p = p, A = A, gamma = gamma, B = B, mu = mu, x_wt = x_wt,
       adj = (A != 0) * 1)                    # 真实邻接（行 j <- 列 i）
}

# 稳态闭式：x*(u) = -B^{-1}(mu + u)
steady_state <- function(sys, U) {
  # U: N x p 扰动矩阵；返回 N x p 稳态
  -t(solve(sys$B, t(sweep(U, 2, -sys$mu, FUN = "+"))))  # = -B^{-1}(mu + U)
}

## ------------------------------------------------------- 设计矩阵 / 基函数 ----
M_ord <- 2L
psi_row <- function(xvec) as.vector(sapply(xvec, function(x) x^(seq_len(M_ord))))  # length p*M
# 增广行 [1, psi] 用于 D-最优（含截距）
aug_row <- function(xvec) c(1, psi_row(xvec))

## ------------------------------------------------- 三种扰动设计策略 ----
# 公共候选池
make_pool <- function(p, n_pool = 4000) {
  matrix(runif(n_pool * p, -0.4, 0.8), n_pool, p)
}

design_random <- function(sys, N, pool) {
  pool[seq_len(N), , drop = FALSE]
}

design_maximin <- function(sys, N, pool) {
  # 贪心最大化最小欧氏间距（u 空间空间填充）
  idx <- 1L                                   # 从第一个起
  mind <- as.vector(sqrt(rowSums((sweep(pool, 2, pool[1, ]))^2)))
  for (k in 2:N) {
    nxt <- which.max(mind)
    idx <- c(idx, nxt)
    d_new <- sqrt(rowSums((sweep(pool, 2, pool[nxt, ]))^2))
    mind <- pmin(mind, d_new)
    mind[idx] <- -Inf
  }
  pool[idx, , drop = FALSE]
}

design_dopt <- function(sys, N, pool, lambda = 1e-2) {
  # 序贯 D-最优：贪心最大化候选行预测方差 psi' Minv psi (eq.3)
  q <- 1 + sys$p * M_ord                       # 增广维度
  X_pool <- steady_state(sys, pool)            # 候选稳态（闭式，快）
  Phi_pool <- t(apply(X_pool, 1, aug_row))     # n_pool x q 增广特征
  # 种子：u=0（野生型）
  seed_idx <- 1L
  pool[1, ] <- 0                               # 第一个候选设为 WT
  X_pool[1, ] <- sys$x_wt
  Phi_pool[1, ] <- aug_row(sys$x_wt)
  Minv <- diag(1 / lambda, q)                  # (lambda I)^{-1}
  add_point <- function(idx) {
    phi <- Phi_pool[idx, ]
    # Sherman-Morrison 秩一更新 Minv
    Mv <- Minv %*% phi
    denom <- as.numeric(1 + phi %*% Mv)
    Minv <<- Minv - (Mv %*% t(Mv)) / denom
  }
  add_point(seed_idx)
  sel <- seed_idx
  avail <- rep(TRUE, nrow(pool)); avail[seed_idx] <- FALSE
  for (k in 2:N) {
    # 预测方差 score = rowSums((Phi %*% Minv) * Phi)
    PM <- Phi_pool %*% Minv
    score <- rowSums(PM * Phi_pool)
    score[!avail] <- -Inf
    nxt <- which.max(score)
    add_point(nxt)
    sel <- c(sel, nxt); avail[nxt] <- FALSE
  }
  pool[sel, , drop = FALSE]
}

## ------------------------------------------------------ 逐节点 ADSIHT 推断 ----
infer_network <- function(sys, U, X) {
  p <- sys$p
  N <- nrow(U)
  # 设计矩阵 Psi (N x pM)
  Psi <- t(apply(X, 1, psi_row))
  group <- rep(seq_len(p), each = M_ord)
  # 中心化
  Psi_bar <- colMeans(Psi)
  Psi_c <- sweep(Psi, 2, Psi_bar)
  # 标准化
  sdv <- apply(Psi_c, 2, sd); sdv[sdv < 1e-10] <- 1e-10
  Psi_cs <- sweep(Psi_c, 2, sdv, FUN = "/")

  adj_est <- matrix(0, p, p)
  theta_hat <- matrix(0, p * M_ord, p)        # 列 j: 节点 j 的 pM 系数
  for (j in seq_len(p)) {
    rhs <- -(U[, j] - mean(U[, j]))
    fit <- tryCatch(
      ADSIHT(Psi_cs, matrix(rhs), group, ic.type = "dsic"),
      error = function(e) NULL)
    if (is.null(fit)) next
    best <- which.min(fit$ic)
    th_s <- fit$beta[, best]
    th <- th_s / sdv                          # 反标准化
    theta_hat[, j] <- th
    # 边判定：组 L2 范数
    gnorm <- sapply(seq_len(p), function(i)
      sqrt(sum(th[((i - 1) * M_ord + 1):(i * M_ord)]^2)))
    adj_est[j, gnorm >= 1e-8] <- 1
  }
  diag(adj_est) <- 0
  list(adj_est = adj_est, theta_hat = theta_hat, Psi_bar = Psi_bar)
}

## ---------------------------------------------------------------- 指标 ----
edge_metrics <- function(est, true) {
  p <- nrow(true)
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

# Jacobian (线性 GLV: J_ji = theta_ji1，因 d/dx (theta1 x + theta2 x^2)|_{x_wt})
jac_rmse <- function(sys, theta_hat) {
  p <- sys$p
  Jhat <- matrix(0, p, p)
  for (j in seq_len(p)) for (i in seq_len(p)) {
    th <- theta_hat[((i - 1) * M_ord + 1):(i * M_ord), j]
    xw <- sys$x_wt[i]
    Jhat[j, i] <- sum(seq_len(M_ord) * xw^(seq_len(M_ord) - 1) * th)
  }
  off <- which(row(sys$A) != col(sys$A))
  sqrt(mean((Jhat[off] - sys$A[off])^2))
}

coef_l2 <- function(sys, theta_hat) {
  # 真值: theta_ji = (A_ji, 0)
  p <- sys$p
  err <- 0
  for (j in seq_len(p)) {
    th_true <- numeric(p * M_ord)
    for (i in seq_len(p)) th_true[(i - 1) * M_ord + 1] <- sys$A[j, i]
    err <- err + sqrt(sum((theta_hat[, j] - th_true)^2))
  }
  err / p
}

## ----------------------------------------------------------- 主实验循环 ----
N_grid <- c(12, 16, 20, 30, 40, 60)
strategies <- c("random", "maximin", "dopt")
R <- 20
sigma <- 0.04                                  # 稳态观测测量噪声（绝对尺度，对齐 sindy_ss σ≈0.03）

rows <- list()
for (seed in seq_len(R)) {
  sys <- make_system(p = 8, n_in = 2, seed = 100 + seed)
  pool <- make_pool(sys$p, n_pool = 4000)     # 共享候选池
  for (N in N_grid) {
    for (strat in strategies) {
      # 设计基于模型预测（无噪），观测含测量噪声——符合实际：设计时不知噪声实现
      U <- switch(strat,
                  random  = design_random(sys, N, pool),
                  maximin = design_maximin(sys, N, pool),
                  dopt    = design_dopt(sys, N, pool))
      X_true <- steady_state(sys, U)
      X <- X_true + matrix(rnorm(length(X_true), sd = sigma), nrow(X_true))
      res <- infer_network(sys, U, X)
      m <- edge_metrics(res$adj_est, sys$adj)
      rows[[length(rows) + 1]] <- data.frame(
        seed = seed, N = N, strategy = strat,
        Pr = m["Pr"], Re = m["Re"], MCC = m["MCC"],
        CoefL2 = coef_l2(sys, res$theta_hat),
        JacRMSE = jac_rmse(sys, res$theta_hat))
    }
  }
  cat("seed", seed, "done\n")
}
df <- do.call(rbind, rows)
rownames(df) <- NULL

dir.create("results/sim_results", showWarnings = FALSE, recursive = TRUE)
write.csv(df, "results/sim_results/design_comparison.csv", row.names = FALSE)

## ---------------------------------------------------------------- 汇总 ----
agg <- aggregate(cbind(MCC, Pr, Re, CoefL2, JacRMSE) ~ strategy + N, df, mean)
agg <- agg[order(agg$N, agg$strategy), ]
cat("\n===== mean over", R, "seeds =====\n")
print(agg, row.names = FALSE)
cat("\nSaved: results/sim_results/design_comparison.csv\n")

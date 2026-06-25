rm(list = ls())

################################################################################
# Fig2b_design_linear.R  —  Fig2b：线性稳态系统中的最优扰动设计
#
# 比较三种扰动设计策略在固定预算 N 下的网络恢复质量：
#   (1) random   —  Uniform 随机扰动（现状基线）
#   (2) maximin  —  u 空间贪心空间填充（Latin-hypercube 风格）
#   (3) dopt     —  AlgDesign::optFederov() 精确 D-optimal design（已有工具）
#
# 系统:  稀疏可加线性 ODE  dx_j/dt = r_j + sum_i a_ji x_i - gamma_j x_j + u_j，
#        稳态闭式  x*(u) = (diag(gamma) - A)^{-1}(r + u)   (§1.1 (2')，与 Fig1c/Fig2a 一致)
# 回归:  逐节点 ADSIHT（中心化 + 标准化），单项式基 M=2
# 指标:  MCC / Precision / Recall / CoefL2 / JacRMSE，扫描 N，重复 R 个种子
#
# Input:   none（自生成模拟数据）
# Output:  results/sim_results/Fig2b_design_linear.csv  — 每 (strategy,N,seed) 指标
# Plot:    sim_script/02_scaling_design/Fig2.R 中 Fig2b 对象
################################################################################

suppressMessages({
  library(ADSIHT)
  library(AlgDesign)
})

set.seed(1)

## ---------------------------------------------------------------- 真值系统 ----
# 生成稀疏可加线性 ODE：A 为 off-diagonal 互作，gamma 为自调节。对角占优
# (gamma_j > sum_i |a_ji|) 保证 M = diag(gamma) - A 可逆且稳态稳定。
make_system <- function(p = 8, n_in = 2, seed = 1) {
  set.seed(seed)
  A <- matrix(0, p, p)                       # A[j,i]: i -> j 的可加线性互作强度
  for (j in seq_len(p)) {
    src <- sample(setdiff(seq_len(p), j), n_in)
    # 较强互作区间：主动设计的增益依赖足够的耦合信号（弱耦合下三种设计持平）。
    A[j, src] <- runif(n_in, 0.3, 0.8) * sample(c(-1, 1), n_in, replace = TRUE)
  }
  gamma <- rowSums(abs(A)) + runif(p, 1.0, 1.5)  # 对角占优自调节
  r <- runif(p, 0.8, 1.5)                     # 内禀项，使 x_wt > 0
  M <- diag(gamma) - A                        # 稳态算子：M x* = r + u
  x_wt <- as.vector(solve(M, r))
  list(p = p, A = A, gamma = gamma, M = M, r = r, x_wt = x_wt,
       adj = (A != 0) * 1)                    # 真实邻接（行 j <- 列 i）
}

# 可加线性稳态：0 = r + A x - gamma x + u  =>  x*(u) = (diag(gamma) - A)^{-1}(r + u)
steady_state <- function(sys, U) {
  # U: N x p 扰动矩阵；返回 N x p 稳态
  t(solve(sys$M, t(sweep(U, 2, sys$r, FUN = "+"))))   # = M^{-1}(r + U)
}

## ------------------------------------------------------- 设计矩阵 / 基函数 ----
M_ord <- 2L
psi_row <- function(xvec) as.vector(sapply(xvec, function(x) x^(seq_len(M_ord))))  # length p*M
# PSS 特征作为模型项传入；AlgDesign 的显式公式 ~ . 同时估计截距。
design_row <- function(xvec) psi_row(xvec)

## ------------------------------------------------- 三种扰动设计策略 ----
# 扰动幅度范围放宽到 u ~ U[-0.4, 0.8]，给 D-最优更大的候选空间以体现设计增益
# （噪声、基函数、横轴刻度仍与 Fig1c/Fig2a 对齐）。
u_lo <- -0.4
u_hi <- 0.8

# 公共候选池
make_pool <- function(p, n_pool = 4000) {
  matrix(runif(n_pool * p, u_lo, u_hi), n_pool, p)
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

prepare_design_features <- function(Phi) {
  Phi <- sweep(Phi, 2, colMeans(Phi))
  s <- apply(Phi, 2, sd)
  s[!is.finite(s) | s < 1e-10] <- 1
  Phi <- sweep(Phi, 2, s, "/")
  q <- qr(Phi, tol = 1e-9)
  if (q$rank < ncol(Phi)) Phi <- Phi[, sort(q$pivot[seq_len(q$rank)]), drop = FALSE]
  Phi
}

design_dopt <- function(sys, N, pool, seed) {
  # 有限候选集 exact D-optimal design；WT（第 1 行）作为所有策略共享条件。
  X_pool <- steady_state(sys, pool)
  Phi <- prepare_design_features(t(apply(X_pool, 1, design_row)))
  if (N < ncol(Phi) + 1L) return(NULL)  # 加截距后的 exact D-opt 模型秩
  set.seed(seed)
  fit <- AlgDesign::optFederov(
    frml = ~ ., data = as.data.frame(Phi), nTrials = N, criterion = "D",
    augment = TRUE, rows = 1L, maxIteration = 100, nRepeats = 1
  )
  pool[fit$rows, , drop = FALSE]
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
    if (is.null(fit) || length(fit$ic) == 0L) next
    best <- which.min(fit$ic)
    if (length(best) == 0L) next               # 退化解：该节点不选任何边
    th_s <- fit$beta[, best]
    if (length(th_s) != length(sdv)) next
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

# Jacobian (可加线性: J_ji = theta_ji1，因 d/dx (theta1 x + theta2 x^2)|_{x_wt})
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
# 样本预算沿用 Fig1c/Fig2a 的重标度刻度 N = ceil(k * s log p)，横轴用 N/(s log p)。
p_design <- 8L
n_in <- 2L
slogp <- n_in * log(p_design)
# 设计增益主要出现在数据稀缺（小预算）区间，故横轴沿用 Fig1c 的 N/(s log p)
# 刻度，但向下覆盖到 p*M 以下的欠定区；退化拟合由 infer_network 的守卫处理。
k_grid <- c(1.5, 2, 2.5, 3, 4, 6, 8)
N_grid <- unique(ceiling(k_grid * slogp))
strategies <- c("random", "maximin", "dopt")
R <- 30
snr_level <- 30                                # 与 Fig1c 一致的相对噪声水平

rows <- list()
for (seed in seq_len(R)) {
  sys <- make_system(p = p_design, n_in = n_in, seed = 100 + seed)
  pool <- make_pool(sys$p, n_pool = 4000)     # 共享候选池
  pool[1, ] <- 0                              # 三种策略共享 WT 条件
  # 噪声相对 SNR=30，但用一次大样本随机参考确定该系统的绝对 sigma，使三种设计
  # 共享同一测量噪声（测量噪声不依赖设计，保证公平比较）。
  U_ref <- make_pool(sys$p, n_pool = 500)
  X_ref <- steady_state(sys, U_ref)
  signal_scale <- mean(apply(X_ref, 2, sd))
  sigma <- signal_scale / snr_level
  for (N in N_grid) {
    for (strat in strategies) {
      # 设计基于模型预测（无噪），观测含测量噪声——符合实际：设计时不知噪声实现
      U <- switch(strat,
                  random  = design_random(sys, N, pool),
                  maximin = design_maximin(sys, N, pool),
                  dopt    = design_dopt(sys, N, pool,
                                        seed = 100000L + 100L * seed + N))
      if (is.null(U)) next
      X_true <- steady_state(sys, U)
      X <- X_true + matrix(rnorm(length(X_true), sd = sigma), nrow(X_true))
      res <- infer_network(sys, U, X)
      m <- edge_metrics(res$adj_est, sys$adj)
      rows[[length(rows) + 1]] <- data.frame(
        seed = seed, p = sys$p, s = n_in, N = N, N_over_slogp = N / slogp,
        strategy = strat, snr = snr_level, sigma = sigma,
        signal_scale = signal_scale,
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
write.csv(df, "results/sim_results/Fig2b_design_linear.csv", row.names = FALSE)

## ---------------------------------------------------------------- 汇总 ----
agg <- aggregate(cbind(MCC, Pr, Re, CoefL2, JacRMSE) ~ strategy + N, df, mean)
agg <- agg[order(agg$N, agg$strategy), ]
cat("\n===== mean over", R, "seeds =====\n")
print(agg, row.names = FALSE)
cat("\nSaved: results/sim_results/Fig2b_design_linear.csv\n")

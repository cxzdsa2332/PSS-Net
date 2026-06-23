rm(list = ls())

################################################################################
# Fig2a_highdim_sample_complexity.R -- PSS-Net high-dimensional sample complexity
#
# Fig2a connects the 8-node demonstration scale used in Fig1 with larger sparse
# systems. It scans p in {8, 50, 100}; p = 8 uses s = 2 incoming edges per node,
# matching the small-network setting, while larger systems use s = 3.
# For each setting, sample budgets follow N ≈ k * s * log(p).
# 用逐节点 ADSIHT 恢复网络，考察 MCC 是否随重标度变量 N/(s·log p) 坍缩到同一曲线
# （§11 (SC): N≳s·log(p/s)+s·s0·log(eM) 的经验印证）。
#
# 稳态：用 gLV/加性共享的闭式线性稳态 x*(u)=solve(diag(γ)-A, r+u)（§1.1 (2')），
#       对角占优保证可解、正值；避免高维 ODE 积分开销。
# 基：单项式 M=2。回归：逐节点 ADSIHT（中心化+标准化）。
#
# Output:  results/sim_results/Fig2a_highdim_sample_complexity.csv
################################################################################

suppressMessages(library(ADSIHT))
set.seed(1)
M_ord <- 2L

## ---------------------------------------------------------------- 真值系统 ----
# 稀疏 A（每节点 s 入边，弱耦合），对角占优 γ 保证 (diag(γ)-A) 可逆且 x*>0
make_system <- function(p, s_in, seed) {
  set.seed(seed)
  A <- matrix(0, p, p)
  for (j in seq_len(p)) {
    src <- sample(setdiff(seq_len(p), j), s_in)
    A[j, src] <- runif(s_in, 0.10, 0.30) * sample(c(-1, 1), s_in, replace = TRUE)
  }
  gamma <- rowSums(abs(A)) + runif(p, 1.0, 1.5)
  r <- runif(p, 0.8, 1.5)
  list(p = p, s = s_in, A = A, gamma = gamma, r = r, adj = (A != 0) * 1)
}

ss_lin <- function(sys, U) t(apply(U, 1, function(u)
  as.numeric(solve(diag(sys$gamma) - sys$A, sys$r + u))))

## ------------------------------------------------- 逐节点 ADSIHT 推断 ----
psi_row <- function(xv) as.vector(sapply(xv, function(x) x^(seq_len(M_ord))))

infer_adj <- function(sys, U, X) {
  p <- sys$p
  Psi <- t(apply(X, 1, psi_row))
  group <- rep(seq_len(p), each = M_ord)
  Psi_c <- sweep(Psi, 2, colMeans(Psi))
  sdv <- apply(Psi_c, 2, sd); sdv[sdv < 1e-10] <- 1e-10
  Psi_cs <- sweep(Psi_c, 2, sdv, FUN = "/")
  adj <- matrix(0, p, p)
  for (j in seq_len(p)) {
    fit <- tryCatch(ADSIHT(Psi_cs, matrix(-(U[, j] - mean(U[, j]))), group,
                           ic.type = "dsic"), error = function(e) NULL)
    if (is.null(fit)) next
    th <- fit$beta[, which.min(fit$ic)]
    gn <- sapply(seq_len(p), function(i)
      sqrt(sum(th[((i - 1) * M_ord + 1):(i * M_ord)]^2)))
    adj[j, gn >= 1e-8] <- 1
  }
  diag(adj) <- 0
  adj
}

mcc_of <- function(est, true) {
  off <- which(row(true) != col(true))
  e <- est[off]; t <- true[off]
  TP <- sum(e == 1 & t == 1); FP <- sum(e == 1 & t == 0)
  TN <- sum(e == 0 & t == 0); FN <- sum(e == 0 & t == 1)
  pr <- ifelse(TP + FP == 0, 0, TP / (TP + FP))
  re <- ifelse(TP + FN == 0, 0, TP / (TP + FN))
  den <- sqrt(as.numeric(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  c(Pr = pr, Re = re, MCC = ifelse(den == 0, 0, (TP * TN - FP * FN) / den))
}

## ----------------------------------------------------------- 主循环 ----
p_grid <- c(8, 50, 100)
k_grid <- c(1, 2, 3, 5, 8, 10, 12, 14, 16, 20, 24, 30)  # N = ceiling(k * s * log p)
R <- 5
sigma <- 0.03

rows <- list()
for (p in p_grid) {
  s_in <- if (p <= 8L) 2L else 3L
  base <- s_in * log(p)
  N_grid <- unique(pmax(M_ord * 2L, ceiling(k_grid * base)))
  for (seed in seq_len(R)) {
    sys <- make_system(p, s_in, seed = 700 + seed)
    for (N in N_grid) {
      U <- matrix(runif(N * p, -0.3, 0.5), N, p); U[1, ] <- 0
      X <- ss_lin(sys, U)
      ok <- apply(X, 1, function(r) all(is.finite(r)) && all(r > 0))
      if (sum(ok) < M_ord * 2L) next
      Xo <- X + matrix(rnorm(length(X), sd = sigma), nrow(X))
      adj <- infer_adj(sys, U[ok, , drop = FALSE], Xo[ok, , drop = FALSE])
      m <- mcc_of(adj, sys$adj)
      rows[[length(rows) + 1]] <- data.frame(
        p = p, s = s_in, N = N, n_eff = sum(ok),
        N_over_slogp = N / base, seed = seed,
        Pr = m["Pr"], Re = m["Re"], MCC = m["MCC"])
    }
    cat(sprintf("p=%3d seed=%d done\n", p, seed))
  }
}
df <- do.call(rbind, rows); rownames(df) <- NULL

dir.create("results/sim_results", showWarnings = FALSE, recursive = TRUE)
out <- "results/sim_results/Fig2a_highdim_sample_complexity.csv"
write.csv(df, out, row.names = FALSE)

agg <- aggregate(MCC ~ p + N + N_over_slogp, df, mean)
agg <- agg[order(agg$p, agg$N), ]
cat("\n===== 高维稀疏恢复（", R, " seeds, mean MCC）=====\n", sep = "")
print(agg, row.names = FALSE)
cat("Saved:", out, "\n")

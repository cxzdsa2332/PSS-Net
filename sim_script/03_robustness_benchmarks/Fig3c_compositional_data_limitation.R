rm(list = ls())

################################################################################
# Fig3c_compositional_data_limitation.R -- compositional-data limitation
#
# 致命缺口 #2：真实 16S 是组成型（只知相对丰度，总量未知），而 PSS 稳态关系
# (2') 写在【绝对】丰度上。本脚本量化组成型化对网络恢复的破坏，并比较常见处理：
#   abs       绝对丰度（oracle 上界）
#   rel       相对丰度（闭合 z=x/Σx，朴素，总量未知）
#   clr       中心对数比 CLR(z)=log z − mean(log z)
#   rel_x_T   相对丰度 × 实测总量（T 含噪，如 qPCR/OD600 校正 → 重建绝对）
#
# 真值：gLV/加性共享闭式稳态（对角占优，正值）。基：单项式 M=2，逐节点 ADSIHT。
#
# Output:  results/sim_results/Fig3c_compositional_data_limitation.csv
################################################################################

suppressMessages(library(ADSIHT))
set.seed(1)
M_ord <- 2L

make_system <- function(p = 10, s_in = 3, seed = 1) {
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
inputs <- c("abs", "rel", "clr", "rel_x_T")
R <- 10; N <- 200; p <- 10; sigma_rel <- 0.05; T_cv <- 0.15
rows <- list()
for (seed in seq_len(R)) {
  sys <- make_system(p = p, s_in = 3, seed = 500 + seed)
  U <- matrix(runif(N * p, -0.3, 0.5), N, p); U[1, ] <- 0
  X <- ss_lin(sys, U)
  ok <- apply(X, 1, function(r) all(is.finite(r)) && all(r > 0))
  Xok <- X[ok, , drop = FALSE]; Uok <- U[ok, , drop = FALSE]
  # 含噪绝对丰度（测量噪声，乘性对数正态更贴近测序计数）
  Xabs <- Xok * exp(matrix(rnorm(length(Xok), sd = sigma_rel), nrow(Xok)))
  Ttot <- rowSums(Xabs)
  Z <- Xabs / Ttot                                   # 相对丰度（闭合）
  CLR <- log(Z) - rowMeans(log(Z))                   # 中心对数比
  That <- Ttot * exp(rnorm(nrow(Xabs), sd = T_cv))   # 实测总量（含噪，如 qPCR）
  Xhat <- Z * That                                   # 相对 × 实测总量 → 重建绝对
  for (inp in inputs) {
    Xin <- switch(inp, abs = Xabs, rel = Z, clr = CLR, rel_x_T = Xhat)
    adj <- infer_adj(sys, Uok, Xin)
    m <- mcc_of(adj, sys$adj)
    rows[[length(rows) + 1]] <- data.frame(
      seed = seed, input = inp, n_eff = sum(ok),
      Pr = m["Pr"], Re = m["Re"], MCC = m["MCC"])
  }
  cat("seed", seed, "done\n")
}
df <- do.call(rbind, rows); rownames(df) <- NULL

dir.create("results/sim_results", showWarnings = FALSE, recursive = TRUE)
out <- "results/sim_results/Fig3c_compositional_data_limitation.csv"
write.csv(df, out, row.names = FALSE)

agg <- aggregate(cbind(Pr, Re, MCC) ~ input, df, mean)
agg <- agg[match(inputs, agg$input), ]
cat("\n===== 组成型数据网络恢复（p=", p, ", ", R, " seeds, mean）=====\n", sep = "")
print(agg, row.names = FALSE)
cat("\n判读：abs 为上界；rel/clr 衡量组成型化的破坏；rel_x_T 检验实测总量校正能否复原。\n")
cat("Saved:", out, "\n")

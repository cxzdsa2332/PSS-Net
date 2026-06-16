rm(list = ls())

################################################################################
# pss_net_joint_smalln.R  —  PSS-Net: 凸显联合估计优势的有利设定
#
# 联合(块对角 + 全局 DSIC)的增益只来自"跨节点共享稀疏度/噪声的全局模型选择"。
# 故构造对其最有利的数据：
#   (1) 均匀入度——每个目标恰好 s 条入边 → 全局共享稀疏度对所有节点都正确；
#   (2) 同质信号——所有 |A_ji| 相等 → 单一全局水平well-specified；
#   (3) 极小 N——逐节点 DSIC 样本太少、选择不稳；联合在 p·N 上池化更稳；
#   (4) 大 p、较高噪声——放大池化收益、加剧逐节点的不稳定。
#
# 分组用 v0.1.txt 方案：p*p 组、每组重复 M（CLAUDE.md 固定规则）。
# Output:  results/sim_results/joint_smalln.csv
################################################################################

suppressMessages(library(ADSIHT))
set.seed(1)
M_ord <- 2L

# 均匀入度 + 同质信号幅度
make_uniform <- function(p, s_in, amp = 0.25, seed = 1) {
  set.seed(seed)
  A <- matrix(0, p, p)
  for (j in seq_len(p)) {
    src <- sample(setdiff(seq_len(p), j), s_in)
    A[j, src] <- amp * sample(c(-1, 1), s_in, replace = TRUE)   # 同质幅度
  }
  gamma <- rowSums(abs(A)) + runif(p, 1.0, 1.5)
  r <- runif(p, 0.8, 1.5)
  list(p = p, s = s_in, A = A, gamma = gamma, r = r, adj = (A != 0) * 1)
}

ss_lin <- function(sys, U) t(apply(U, 1, function(u)
  as.numeric(solve(diag(sys$gamma) - sys$A, sys$r + u))))

psi_row <- function(xv) as.vector(sapply(xv, function(x) x^(seq_len(M_ord))))

build_design <- function(U, X) {
  p <- ncol(X)
  Psi <- t(apply(X, 1, psi_row))
  Psi_c <- sweep(Psi, 2, colMeans(Psi))
  sdv <- apply(Psi_c, 2, sd); sdv[sdv < 1e-10] <- 1e-10
  Psi_cs <- sweep(Psi_c, 2, sdv, FUN = "/")
  rhs <- sapply(seq_len(p), function(j) -(U[, j] - mean(U[, j])))
  list(Psi_cs = Psi_cs, rhs = rhs, p = p)
}

adj_from_theta <- function(th, p) {
  gn <- sapply(seq_len(p), function(i)
    sqrt(sum(th[((i - 1) * M_ord + 1):(i * M_ord)]^2)))
  as.integer(gn >= 1e-8)
}

infer_nodewise <- function(d) {
  p <- d$p; group <- rep(seq_len(p), each = M_ord)
  adj <- matrix(0, p, p)
  for (j in seq_len(p)) {
    fit <- tryCatch(ADSIHT(d$Psi_cs, matrix(d$rhs[, j]), group,
                           ic.type = "dsic"), error = function(e) NULL)
    if (is.null(fit)) next
    adj[j, ] <- adj_from_theta(fit$beta[, which.min(fit$ic)], p)
  }
  diag(adj) <- 0; adj
}

# 联合：X=I_p⊗Ψ，p*p 组、每组 M 次（CLAUDE.md 固定规则）
infer_joint <- function(d) {
  p <- d$p
  Xbig <- kronecker(diag(p), d$Psi_cs)
  Ybig <- as.vector(d$rhs)
  group <- rep(seq_len(p * p), each = M_ord)
  fit <- tryCatch(ADSIHT(Xbig, matrix(Ybig), group, ic.type = "dsic"),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  beta <- fit$beta[, which.min(fit$ic)]; pM <- p * M_ord
  adj <- matrix(0, p, p)
  for (j in seq_len(p))
    adj[j, ] <- adj_from_theta(beta[((j - 1) * pM + 1):(j * pM)], p)
  diag(adj) <- 0; adj
}

metrics <- function(est, true) {
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
p_grid <- c(50)
s_in <- 3L
k_grid <- c(1.5, 2.0, 2.5, 3.4)        # N≈18,24,30,40：部分恢复中间区
R <- 8
sigma <- 0.04                           # 中等噪声（落在 MCC~0.2–0.5 的可分辨区）
amp <- 0.30                             # 同质、稍强信号

rows <- list()
for (p in p_grid) {
  base <- s_in * log(p)
  N_grid <- unique(ceiling(k_grid * base))
  for (seed in seq_len(R)) {
    sys <- make_uniform(p, s_in, amp = amp, seed = 900 + seed)
    for (N in N_grid) {
      U <- matrix(runif(N * p, -0.3, 0.5), N, p); U[1, ] <- 0
      X <- ss_lin(sys, U)
      ok <- apply(X, 1, function(rr) all(is.finite(rr)) && all(rr > 0))
      Xo <- (X + matrix(rnorm(length(X), sd = sigma), nrow(X)))[ok, , drop = FALSE]
      d <- build_design(U[ok, , drop = FALSE], Xo)
      adj_nw <- infer_nodewise(d)
      adj_jt <- infer_joint(d)
      m_nw <- metrics(adj_nw, sys$adj)
      m_jt <- if (is.null(adj_jt)) c(Pr = NA, Re = NA, MCC = NA) else
        metrics(adj_jt, sys$adj)
      rows[[length(rows) + 1]] <- data.frame(
        p = p, N = N, seed = seed,
        MCC_nodewise = m_nw["MCC"], MCC_joint = m_jt["MCC"],
        Pr_nodewise = m_nw["Pr"], Pr_joint = m_jt["Pr"],
        Re_nodewise = m_nw["Re"], Re_joint = m_jt["Re"])
      cat(sprintf("p=%d N=%2d seed=%d | MCC nw=%.3f jt=%.3f\n",
                  p, N, seed, m_nw["MCC"], m_jt["MCC"]))
    }
  }
}
df <- do.call(rbind, rows); rownames(df) <- NULL

dir.create("results/sim_results", showWarnings = FALSE, recursive = TRUE)
write.csv(df, "results/sim_results/joint_smalln.csv", row.names = FALSE)

agg <- aggregate(cbind(MCC_nodewise, MCC_joint, Pr_nodewise, Pr_joint,
                       Re_nodewise, Re_joint) ~ p + N, df, mean)
agg <- agg[order(agg$p, agg$N), ]
cat("\n===== 有利设定（均匀入度+同质信号+极小N+高噪），", R, " seeds =====\n", sep = "")
print(agg, row.names = FALSE)
cat("\nSaved: results/sim_results/joint_smalln.csv\n")

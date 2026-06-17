rm(list = ls())

################################################################################
# pss_net_scalefree.R  —  PSS-Net: scale-free 网络下 联合 vs 逐节点 + hub 识别
#
# 重新设计真值网络为 scale-free（出度幂律、不固定 degree）：
#   - 用偏好连接(preferential attachment)生成有向边：目标 j 选源 i 的概率 ∝ 源的
#     当前出度 → 少数 hub 源积累高出度（影响很多目标），多数源低出度；
#   - 入度可变（Poisson），不固定。
# 这是"跨目标共享源支撑"最强的结构，检验联合(块对角, p*p 组)能否借力 + 找到 hub。
#
# 分组：CLAUDE.md 固定规则 group = rep(1:(p*p), each=M)。
# 指标：逐边 MCC；hub 识别——估计出度 vs 真出度 Spearman、top-k 命中。
# Output:  results/sim_results/scalefree_compare.csv
################################################################################

suppressMessages(library(ADSIHT))
set.seed(1)
M_ord <- 2L

# scale-free：偏好连接生成出度幂律；入度 ~ 1+Poisson
make_scalefree <- function(p, avg_in = 2, seed = 1) {
  set.seed(seed)
  A <- matrix(0, p, p)
  outdeg <- rep(1, p)                                  # +1 平滑，作偏好权重
  for (j in seq_len(p)) {
    k_j <- min(1 + rpois(1, avg_in - 1), p - 1)        # 可变入度
    chosen <- integer(0)
    for (t in seq_len(k_j)) {
      w <- outdeg; w[c(j, chosen)] <- 0                # 不自环、不重复
      i <- sample.int(p, 1, prob = w)                  # 概率 ∝ 出度（偏好连接）
      chosen <- c(chosen, i)
      A[j, i] <- runif(1, 0.15, 0.35) * sample(c(-1, 1), 1)
      outdeg[i] <- outdeg[i] + 1
    }
  }
  gamma <- rowSums(abs(A)) + runif(p, 1.0, 1.5)
  r <- runif(p, 0.8, 1.5)
  list(p = p, A = A, gamma = gamma, r = r, adj = (A != 0) * 1,
       outdeg_true = colSums(A != 0))
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

infer_joint <- function(d) {
  p <- d$p
  Xbig <- kronecker(diag(p), d$Psi_cs)
  Ybig <- as.vector(d$rhs)
  group <- rep(seq_len(p * p), each = M_ord)          # v0.1.txt: p*p 组
  fit <- tryCatch(ADSIHT(Xbig, matrix(Ybig), group, ic.type = "dsic"),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  beta <- fit$beta[, which.min(fit$ic)]; pM <- p * M_ord
  adj <- matrix(0, p, p)
  for (j in seq_len(p))
    adj[j, ] <- adj_from_theta(beta[((j - 1) * pM + 1):(j * pM)], p)
  diag(adj) <- 0; adj
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

# hub 识别：出度 Spearman + top-k 命中（k = 真出度最高的 10% 节点）
hub_of <- function(adj_est, sys) {
  out_est <- colSums(adj_est); out_true <- sys$outdeg_true
  k <- max(3, round(0.1 * sys$p))
  true_top <- order(out_true, decreasing = TRUE)[seq_len(k)]
  est_top <- order(out_est, decreasing = TRUE)[seq_len(k)]
  rho <- suppressWarnings(cor(out_est, out_true, method = "spearman"))
  c(rho = ifelse(is.na(rho), 0, rho),
    topk_hit = length(intersect(true_top, est_top)) / k)
}

## ----------------------------------------------------------- 主循环 ----
p <- 50L; avg_in <- 2
k_grid <- c(2, 3, 5)
R <- 5
sigma <- 0.03

rows <- list()
base <- avg_in * log(p)
N_grid <- unique(ceiling(k_grid * base))
for (seed in seq_len(R)) {
  sys <- make_scalefree(p, avg_in = avg_in, seed = 1000 + seed)
  max_out <- max(sys$outdeg_true); n_edge <- sum(sys$adj)
  for (N in N_grid) {
    U <- matrix(runif(N * p, -0.3, 0.5), N, p); U[1, ] <- 0
    X <- ss_lin(sys, U)
    ok <- apply(X, 1, function(rr) all(is.finite(rr)) && all(rr > 0))
    Xo <- (X + matrix(rnorm(length(X), sd = sigma), nrow(X)))[ok, , drop = FALSE]
    d <- build_design(U[ok, , drop = FALSE], Xo)
    adj_nw <- infer_nodewise(d); adj_jt <- infer_joint(d)
    m_nw <- mcc_of(adj_nw, sys$adj)
    m_jt <- if (is.null(adj_jt)) c(Pr = NA, Re = NA, MCC = NA) else mcc_of(adj_jt, sys$adj)
    h_nw <- hub_of(adj_nw, sys)
    h_jt <- if (is.null(adj_jt)) c(rho = NA, topk_hit = NA) else hub_of(adj_jt, sys)
    rows[[length(rows) + 1]] <- data.frame(
      p = p, N = N, seed = seed, n_edge = n_edge, max_outdeg = max_out,
      MCC_nodewise = m_nw["MCC"], MCC_joint = m_jt["MCC"],
      rho_nodewise = h_nw["rho"], rho_joint = h_jt["rho"],
      topk_nodewise = h_nw["topk_hit"], topk_joint = h_jt["topk_hit"])
    cat(sprintf("seed=%d N=%2d | maxOut=%2d | MCC nw=%.3f jt=%.3f | hubRho nw=%.2f jt=%.2f\n",
                seed, N, max_out, m_nw["MCC"], m_jt["MCC"], h_nw["rho"], h_jt["rho"]))
  }
}
df <- do.call(rbind, rows); rownames(df) <- NULL

dir.create("results/sim_results", showWarnings = FALSE, recursive = TRUE)
write.csv(df, "results/sim_results/scalefree_compare.csv", row.names = FALSE)

agg <- aggregate(cbind(MCC_nodewise, MCC_joint, rho_nodewise, rho_joint,
                       topk_nodewise, topk_joint) ~ N, df, mean)
cat("\n===== scale-free（p=50, avg_in=", avg_in, ", ", R, " seeds, 真出度最大≈",
    round(mean(df$max_outdeg)), "）=====\n", sep = "")
print(agg, row.names = FALSE)
cat("\nSaved: results/sim_results/scalefree_compare.csv\n")

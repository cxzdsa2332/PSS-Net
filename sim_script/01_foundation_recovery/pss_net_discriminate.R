rm(list = ls())

################################################################################
# pss_net_discriminate.R  —  PSS-Net: 生成机制可识别性判别（初步）
#
# 验证 note/identifiability_glv_vs_additive.md 的两条结论：
#   (1) 乘性 gLV 与加性线性模型：稳态不可区分（稳态值相同 ⇒ 判别量相同）；
#   (2) 线性 vs 非线性互作：可区分（M=1 vs M=2 的残差下降 / BIC）。
#
# 四种真值（同一稀疏支撑 A）：
#   glv_lin   乘性 gLV， 线性互作       dx=x(r+Ax-γx+u)
#   add_lin   加性，    线性互作       dx=r+Ax-γx+u        （与 glv_lin 稳态应相同）
#   add_quad  加性，    二次互作       f_ji=A x_i + B x_i^2
#   add_monod 加性，    Monod 饱和互作 f_ji=A x_i/(1+0.5 x_i)
#
# 判别量（逐节点 OLS，中心化无截距）：
#   relΔRSS = (RSS_{M=1} - RSS_{M=2}) / RSS_{M=1}   （加二次项的残差相对下降）
#   %BIC→M2 = BIC 选择 M=2 的节点比例
# 期望：线性真值 relΔRSS≈0、%BIC→M2 低；非线性真值二者明显升高；
#       glv_lin 与 add_lin 判别量几乎相同（不可区分）。
#
# Output:  results/sim_results/discriminate_gof.csv
################################################################################

suppressMessages({ library(deSolve); library(ADSIHT) })
set.seed(1)

## ---------------------------------------------------------------- 真值系统 ----
make_system <- function(p = 6, n_in = 2, seed = 1, b_scale = 1) {
  set.seed(seed)
  A <- matrix(0, p, p); B <- matrix(0, p, p)
  for (j in seq_len(p)) {
    src <- sample(setdiff(seq_len(p), j), n_in)
    A[j, src] <- runif(n_in, 0.25, 0.5) * sample(c(-1, 1), n_in, TRUE)
    B[j, src] <- b_scale * runif(n_in, 0.15, 0.30) *
      sample(c(-1, 1), n_in, TRUE)                                       # 二次项
  }
  gamma <- rowSums(abs(A)) + rowSums(abs(B)) + runif(p, 1.2, 1.8)
  r <- runif(p, 1.0, 1.8)
  list(p = p, A = A, B = B, gamma = gamma, r = r, adj = (A != 0) * 1)
}

ss_lin <- function(sys, u) as.numeric(solve(diag(sys$gamma) - sys$A, sys$r + u))

# 通用 ODE 积分至稳态；mult=TRUE 为乘性 gLV，否则加性；form 指定互作形式
integrate_ss <- function(sys, u, mult, form, t_max = 2500) {
  p <- sys$p
  interact <- function(x) {
    v <- numeric(p)
    for (j in seq_len(p)) {
      o <- setdiff(seq_len(p), j); xi <- x[o]
      fij <- switch(form,
        lin   = sys$A[j, o] * xi,
        quad  = sys$A[j, o] * xi + sys$B[j, o] * xi^2,
        monod = sys$A[j, o] * xi / (1 + 0.5 * xi))
      v[j] <- sum(fij)
    }
    v
  }
  deriv <- function(t, x, parms) {
    rate <- sys$r + interact(x) - sys$gamma * x + u
    list(if (mult) x * rate else rate)
  }
  x0 <- pmax(ss_lin(sys, u), 0.05)
  out <- tryCatch(ode(x0, c(0, t_max), deriv, NULL, method = "lsoda",
                      rtol = 1e-9, atol = 1e-11), error = function(e) NULL)
  if (is.null(out)) return(rep(NA, p))
  as.numeric(out[2, -1])
}

steady <- function(sys, U, truth) {
  switch(truth,
    glv_lin   = t(apply(U, 1, function(u) integrate_ss(sys, u, TRUE,  "lin"))),
    add_lin   = t(apply(U, 1, function(u) ss_lin(sys, u))),
    add_quad  = t(apply(U, 1, function(u) integrate_ss(sys, u, FALSE, "quad"))),
    add_monod = t(apply(U, 1, function(u) integrate_ss(sys, u, FALSE, "monod"))))
}

## ------------------------------------- 逐节点 GoF：M=1 vs M=2（OLS 残差/BIC）----
gof <- function(U, X) {
  p <- ncol(X); N <- nrow(X)
  cen <- function(M) sweep(M, 2, colMeans(M))
  X1 <- cen(X)                       # M=1 设计：x_i
  X2 <- cen(cbind(X, X^2))           # M=2 设计：x_i, x_i^2
  rel <- numeric(p); bic2 <- logical(p)
  for (j in seq_len(p)) {
    y <- -(U[, j] - mean(U[, j]))
    rss1 <- sum(lm.fit(X1, y)$residuals^2)
    rss2 <- sum(lm.fit(X2, y)$residuals^2)
    rel[j] <- (rss1 - rss2) / rss1
    b1 <- N * log(rss1 / N) + ncol(X1) * log(N)
    b2 <- N * log(rss2 / N) + ncol(X2) * log(N)
    bic2[j] <- b2 < b1
  }
  c(relRSS = mean(rel), bicM2 = mean(bic2))
}

## ------------------------------------------------- 网络恢复（ADSIHT, M=2）----
M_ord <- 2L
mcc_recover <- function(sys, U, X) {
  p <- sys$p
  Psi <- t(apply(X, 1, function(xv)
    as.vector(sapply(xv, function(x) x^(seq_len(M_ord))))))
  group <- rep(seq_len(p), each = M_ord)
  Psi_c <- sweep(Psi, 2, colMeans(Psi))
  sdv <- apply(Psi_c, 2, sd); sdv[sdv < 1e-10] <- 1e-10
  Psi_cs <- sweep(Psi_c, 2, sdv, "/")
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
  off <- which(row(sys$adj) != col(sys$adj))
  e <- adj[off]; t <- sys$adj[off]
  TP <- sum(e == 1 & t == 1); FP <- sum(e == 1 & t == 0)
  TN <- sum(e == 0 & t == 0); FN <- sum(e == 0 & t == 1)
  den <- sqrt(as.numeric(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  ifelse(den == 0, 0, (TP * TN - FP * FN) / den)
}

## ----------------------------------------------------------- 主循环 ----
# 两种扰动幅度：narrow（窄域，x* 范围小）vs wide（宽域，激发曲率）
# wide 同时加大非线性强度 b_scale，检验"宽扰动下非线性可判别"
truths <- c("glv_lin", "add_lin", "add_quad", "add_monod")
regimes <- list(
  narrow = list(lo = -0.3, hi = 0.5, b_scale = 1),
  wide   = list(lo = -1.0, hi = 2.0, b_scale = 2))
R <- 10; N <- 200; sigma <- 0.03
rows <- list()
for (regime in names(regimes)) {
  rp <- regimes[[regime]]
  for (seed in seq_len(R)) {
    sys <- make_system(p = 6, n_in = 2, seed = 400 + seed, b_scale = rp$b_scale)
    U <- matrix(runif(N * sys$p, rp$lo, rp$hi), N, sys$p); U[1, ] <- 0
    X_glv_lin <- NULL
    for (truth in truths) {
      X <- steady(sys, U, truth)
      ok <- apply(X, 1, function(r) all(is.finite(r)) && all(r > 0))
      if (truth == "glv_lin") X_glv_lin <- X
      Xo <- X + matrix(rnorm(length(X), sd = sigma), nrow(X))
      g <- gof(U[ok, , drop = FALSE], Xo[ok, , drop = FALSE])
      mcc <- mcc_recover(sys, U[ok, , drop = FALSE], Xo[ok, , drop = FALSE])
      ss_diff <- if (truth == "add_lin")
        max(abs(X[ok, ] - X_glv_lin[ok, ])) else NA
      rows[[length(rows) + 1]] <- data.frame(
        regime = regime, seed = seed, truth = truth, n_eff = sum(ok),
        relRSS = g["relRSS"], bicM2 = g["bicM2"], MCC = mcc, ss_diff = ss_diff)
    }
  }
  cat("regime", regime, "done\n")
}
df <- do.call(rbind, rows); rownames(df) <- NULL

dir.create("results/sim_results", showWarnings = FALSE, recursive = TRUE)
write.csv(df, "results/sim_results/discriminate_gof.csv", row.names = FALSE)

agg <- aggregate(cbind(relRSS, bicM2, MCC) ~ regime + truth, df, mean)
agg <- agg[order(agg$regime, match(agg$truth, truths)), ]
cat("\n===== 判别量（", R, "seeds, N=", N, ", p=6）=====\n", sep = "")
print(agg, row.names = FALSE)
for (rg in names(regimes))
  cat(sprintf("[%s] glv_lin vs add_lin 稳态最大差异: %.2e\n", rg,
              mean(df$ss_diff[df$truth == "add_lin" & df$regime == rg])))
cat("\n判读：relRSS/bicM2 越大→越偏非线性互作；glv_lin 与 add_lin 始终几乎相同（不可区分机制）。\n")
cat("Saved: results/sim_results/discriminate_gof.csv\n")

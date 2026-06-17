rm(list = ls())

################################################################################
# Fig1c_adsiht_group_lasso_scaling.R -- Fig1c node-wise benchmark
#
# Purpose: compare node-wise ADSIHT and group lasso across p = 8, 30, 100
#          under sparse additive/gLV-equivalent steady-state systems, with
#          measurement noise fixed at steady-state SNR = 30.
# Input:   none
# Output:  results/sim_results/Fig1c_adsiht_group_lasso_scaling.csv
################################################################################

suppressMessages({
  library(ADSIHT)
  library(grpreg)
})

set.seed(101)
M_ord <- 2L
p_grid <- c(8L, 30L, 100L)
k_grid <- c(4, 8, 12, 16)
R <- 30L
snr_level <- 30

make_system <- function(p, s_in, seed) {
  set.seed(seed)
  A <- matrix(0, p, p)
  for (j in seq_len(p)) {
    src <- sample(setdiff(seq_len(p), j), s_in)
    A[j, src] <- runif(s_in, 0.12, 0.32) * sample(c(-1, 1), s_in, TRUE)
  }
  gamma <- rowSums(abs(A)) + runif(p, 1.0, 1.5)
  r <- runif(p, 0.8, 1.5)
  list(p = p, s_in = s_in, A = A, gamma = gamma, r = r, adj = (A != 0) * 1L)
}

steady_linear <- function(sys, U) {
  solve_mat <- diag(sys$gamma) - sys$A
  t(apply(U, 1, function(u) as.numeric(solve(solve_mat, sys$r + u))))
}

make_basis <- function(X) {
  p <- ncol(X)
  Psi <- matrix(0, nrow(X), p * M_ord)
  for (i in seq_len(p)) {
    for (m in seq_len(M_ord)) {
      Psi[, (i - 1L) * M_ord + m] <- X[, i]^m
    }
  }
  Psi
}

standardize_design <- function(Psi) {
  Psi_bar <- colMeans(Psi)
  Psi_c <- sweep(Psi, 2, Psi_bar)
  Psi_sd <- pmax(apply(Psi_c, 2, sd), 1e-10)
  list(X = sweep(Psi_c, 2, Psi_sd, "/"), center = Psi_bar, scale = Psi_sd)
}

group_norms <- function(beta, p) {
  sapply(seq_len(p), function(i) {
    cols <- (i - 1L) * M_ord + seq_len(M_ord)
    sqrt(sum(beta[cols]^2))
  })
}

fit_adsiht_node <- function(X_cs, y, group, scale_vec, p) {
  fit <- tryCatch(ADSIHT(X_cs, matrix(y), group, ic.type = "dsic"),
                  error = function(e) NULL)
  if (is.null(fit)) return(rep(0, p * M_ord))
  beta <- fit$beta[, which.min(fit$ic)] / scale_vec
  as.numeric(beta)
}

fit_grlasso_node <- function(X_cs, y, group, scale_vec, p) {
  fit <- tryCatch(grpreg(X_cs, y, group = group, penalty = "grLasso"),
                  error = function(e) NULL)
  if (is.null(fit)) return(rep(0, p * M_ord))
  beta_path <- coef(fit)[-1, , drop = FALSE]
  rss <- colSums((matrix(y, nrow(X_cs), ncol(beta_path)) - X_cs %*% beta_path)^2)
  df <- colSums(abs(beta_path) > 1e-10)
  bic <- length(y) * log(pmax(rss, 1e-12) / length(y)) + df * log(length(y))
  beta <- beta_path[, which.min(bic)] / scale_vec
  as.numeric(beta)
}

rank_metrics <- function(score, truth) {
  ord <- order(score, decreasing = TRUE)
  y <- truth[ord]
  pos <- sum(y == 1)
  neg <- sum(y == 0)
  if (pos == 0 || neg == 0) {
    return(c(AUROC = NA_real_, AUPRC = NA_real_))
  }
  rank_pos <- which(y == 1)
  auc <- (sum(rank_pos) - pos * (pos + 1) / 2) / (pos * neg)
  auc <- 1 - auc
  tp <- cumsum(y == 1)
  fp <- cumsum(y == 0)
  recall <- tp / pos
  precision <- tp / pmax(tp + fp, 1)
  keep <- y == 1
  auprc <- sum(diff(c(0, recall[keep])) * precision[keep])
  c(AUROC = auc, AUPRC = auprc)
}

selection_metrics <- function(score, beta_mat, truth_adj, A_true, gamma_true,
                              threshold = 1e-8) {
  p <- nrow(truth_adj)
  est <- (score >= threshold) * 1L
  diag(est) <- 0L
  off <- which(row(truth_adj) != col(truth_adj))
  e <- est[off]
  t <- truth_adj[off]
  TP <- sum(e == 1 & t == 1)
  FP <- sum(e == 1 & t == 0)
  TN <- sum(e == 0 & t == 0)
  FN <- sum(e == 0 & t == 1)
  pr <- ifelse(TP + FP == 0, 0, TP / (TP + FP))
  re <- ifelse(TP + FN == 0, 0, TP / (TP + FN))
  f1 <- ifelse(pr + re == 0, 0, 2 * pr * re / (pr + re))
  den <- sqrt(as.numeric(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  mcc <- ifelse(den == 0, 0, (TP * TN - FP * FN) / den)
  rank_mets <- unname(rank_metrics(score[off], t))

  beta_true <- matrix(0, p, p * M_ord)
  for (j in seq_len(p)) {
    for (i in seq_len(p)) {
      beta_true[j, (i - 1L) * M_ord + 1L] <- if (i == j) -gamma_true[j] else A_true[j, i]
    }
  }
  coef_l2 <- mean(sqrt(rowSums((beta_mat - beta_true)^2)))
  jac_est <- matrix(0, p, p)
  for (j in seq_len(p)) {
    for (i in seq_len(p)) {
      jac_est[j, i] <- beta_mat[j, (i - 1L) * M_ord + 1L]
    }
  }
  jac_true <- A_true
  diag(jac_true) <- -gamma_true
  jac_rmse <- sqrt(mean((jac_est - jac_true)^2))
  edge_weight_rmse <- sqrt(mean((jac_est[off] - A_true[off])^2))

  sign_ok <- c()
  for (j in seq_len(p)) {
    for (i in seq_len(p)) {
      if (i != j && truth_adj[j, i] == 1 && est[j, i] == 1) {
        beta1 <- beta_mat[j, (i - 1L) * M_ord + 1L]
        sign_ok <- c(sign_ok, sign(beta1) == sign(A_true[j, i]))
      }
    }
  }
  sign_acc <- ifelse(length(sign_ok) == 0, NA_real_, mean(sign_ok))
  c(TP = TP, FP = FP, TN = TN, FN = FN, Precision = pr, Recall = re, F1 = f1,
    MCC = mcc, AUROC = rank_mets[1], AUPRC = rank_mets[2], SignAcc = sign_acc,
    CoefL2 = coef_l2, JacRMSE = jac_rmse, EdgeWeightRMSE = edge_weight_rmse)
}

infer_methods <- function(sys, U, X_obs) {
  p <- sys$p
  Psi <- make_basis(X_obs)
  std <- standardize_design(Psi)
  group <- rep(seq_len(p), each = M_ord)
  U_c <- sweep(U, 2, colMeans(U))

  beta_ads <- matrix(0, p, p * M_ord)
  beta_gl <- matrix(0, p, p * M_ord)
  score_ads <- matrix(0, p, p)
  score_gl <- matrix(0, p, p)

  time_ads <- system.time({
    for (j in seq_len(p)) {
      y <- -U_c[, j]
      beta_ads[j, ] <- fit_adsiht_node(std$X, y, group, std$scale, p)
      score_ads[j, ] <- group_norms(beta_ads[j, ], p)
    }
  })[["elapsed"]]

  time_gl <- system.time({
    for (j in seq_len(p)) {
      y <- -U_c[, j]
      beta_gl[j, ] <- fit_grlasso_node(std$X, y, group, std$scale, p)
      score_gl[j, ] <- group_norms(beta_gl[j, ], p)
    }
  })[["elapsed"]]

  diag(score_ads) <- 0
  diag(score_gl) <- 0
  list(
    ADSIHT = list(score = score_ads, beta = beta_ads, runtime = time_ads),
    GroupLasso = list(score = score_gl, beta = beta_gl, runtime = time_gl)
  )
}

rows <- list()
for (p in p_grid) {
  s_in <- if (p <= 8L) 2L else 3L
  base <- s_in * log(p)
  N_grid <- unique(ceiling(k_grid * base))
  for (seed in seq_len(R)) {
    sys <- make_system(p, s_in, seed = 1500 + 100 * p + seed)
    for (N in N_grid) {
      set.seed(9000 + 1000 * p + 10 * seed + N)
      U <- matrix(runif(N * p, -0.3, 0.5), N, p)
      U[1, ] <- 0
      X <- steady_linear(sys, U)
      ok <- apply(X, 1, function(z) all(is.finite(z)) && all(z > 0))
      if (sum(ok) < 8L) next
      X_clean <- X[ok, , drop = FALSE]
      signal_scale <- mean(apply(X_clean, 2, sd))
      sigma_x <- signal_scale / snr_level
      set.seed(9000 + 1000 * p + 100 * seed + N + round(10 * snr_level))
      X_obs <- X_clean + matrix(rnorm(length(X_clean), sd = sigma_x), nrow(X_clean), p)
      X_obs <- pmax(X_obs, 1e-6)
      fit <- infer_methods(sys, U[ok, , drop = FALSE], X_obs)
      for (method in names(fit)) {
        mets <- selection_metrics(fit[[method]]$score, fit[[method]]$beta,
                                  sys$adj, sys$A, sys$gamma)
        rows[[length(rows) + 1L]] <- data.frame(
          p = p, s_in = s_in, N = N, n_eff = sum(ok), seed = seed,
          N_over_slogp = N / base, snr = snr_level, sigma_x = sigma_x,
          signal_scale = signal_scale, method = method,
          t(mets), runtime_sec = fit[[method]]$runtime,
          stringsAsFactors = FALSE
        )
      }
      cat(sprintf("p=%3d seed=%d N=%3d done at SNR=%s\n", p, seed, N, snr_level))
    }
  }
}

df <- do.call(rbind, rows)
rownames(df) <- NULL

out <- "results/sim_results/Fig1c_adsiht_group_lasso_scaling.csv"
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
write.csv(df, out, row.names = FALSE)

cat("\nSaved:", out, "\n")
print(aggregate(cbind(MCC, AUPRC, AUROC, Precision, Recall, CoefL2, JacRMSE, EdgeWeightRMSE, runtime_sec) ~ p + method,
                df, mean, na.rm = TRUE), row.names = FALSE)

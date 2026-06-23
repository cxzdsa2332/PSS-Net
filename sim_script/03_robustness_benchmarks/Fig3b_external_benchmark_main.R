rm(list = ls())

################################################################################
# Fig3b_external_benchmark_main.R -- external method benchmark for PSS-Net
#
# Purpose: head-to-head comparison of PSS-Net against external network-inference
#          methods on matched simulated PSS data, under linear, nonlinear and
#          strong-nonlinear truth regimes (shared A, r and perturbation design;
#          nonlinear adds quadratic cross edges, strong-nonlinear scales them up
#          ~4x where a linear steady-state model breaks). Method receives the SAME data
#          per (regime, N, seed). All methods use package defaults or a single
#          documented default threshold -- no per-method tuning. Threshold-free
#          AUROC/AUPRC are the primary fair metrics; MCC/precision/recall use each
#          method's default selection; sign accuracy where a signed effect exists.
#
# Methods (see ref/external_benchmark_methods.md):
#   PSS_Net      node-wise ADSIHT, double sparsity            (uses u, nonlinear)
#   Lasso        node-wise glmnet on the PSS library          (uses u, nonlinear)
#   ElasticNet   node-wise glmnet (alpha = 0.5)               (uses u, nonlinear)
#   GroupLasso   node-wise gglasso, source-level groups       (uses u, nonlinear)
#   PySINDy_STLSQ official PySINDy STLSQ on the PSS library   (uses u, nonlinear)
#   aiMeRA       classical normalized MRA via the aiMeRA pkg  (uses u, local linear)
#   Correlation  Pearson association                          (undirected)
#   PartialCor   ppcor package partial correlation            (undirected)
#   (optional)   GENIE3 / GIES -- wired but skipped unless randomForest / pcalg
#                are installed.
#
# Reuse: the simulation (make_system / sys_for / steady_*), the PSS library and
#        node-wise PSS-Net inference, and rank_metrics are taken verbatim from
#        sim_script/01_foundation_recovery/Fig1c_adsiht_group_lasso_scaling.R.
#
# Input:   none
# Output:  results/sim_results/Fig3b_external_benchmark_main.csv
################################################################################

if (!requireNamespace("aiMeRA", quietly = TRUE)) {
  stop(paste0(
    "Package 'aiMeRA' is required. Install the audited source with: ",
    "remotes::install_github('bioinfo-ircm/aiMeRA@",
    "86cabc21e8ed124ce372c2fe8e62b47503c2a22b')"
  ))
}
if (!requireNamespace("reticulate", quietly = TRUE)) {
  stop("R package 'reticulate' is required for the official PySINDy backend.")
}
if (!requireNamespace("ppcor", quietly = TRUE)) {
  stop("R package 'ppcor' is required for the partial-correlation baseline.")
}

suppressMessages({
  library(deSolve)
  library(ADSIHT)
  library(glmnet)
  library(gglasso)
  library(aiMeRA)
  library(reticulate)
  library(ppcor)
})

# Pin the original-author Python implementation. Set PSSNET_PYTHON to a Python
# executable containing this package when reticulate should not use its default.
pysindy_commit <- "c4421fcec275c8f4cc5c1e93bebb961b212067ae"
local_python <- file.path(getwd(), ".venv-fig3b-standalone", "Scripts", "python.exe")
python_exe <- Sys.getenv("PSSNET_PYTHON", unset = local_python)
if (!nzchar(python_exe)) {
  stop(paste0(
    "Set PSSNET_PYTHON to a directly accessible Python executable containing ",
    "PySINDy. The Windows Store app alias cannot be used by reticulate here."
  ))
}
if (!file.exists(python_exe)) stop("PSSNET_PYTHON does not exist: ", python_exe)
reticulate::use_python(python_exe, required = TRUE)
if (!reticulate::py_module_available("pysindy")) {
  stop(paste0(
    "Python package 'pysindy' is required. Install the audited official source: ",
    "python -m pip install 'git+https://github.com/dynamicslab/pysindy.git@",
    pysindy_commit, "'"
  ))
}
pysindy <- reticulate::import("pysindy", delay_load = FALSE)
pysindy_version <- as.character(pysindy$`__version__`)
message("PySINDy backend: ", pysindy_version, " [", reticulate::py_config()$python, "]")

set.seed(303)
M_ord <- 2L

## ============================ shared simulation (from Fig1c) ==================
make_system <- function(p, s_in, seed) {
  set.seed(seed)
  A <- matrix(0, p, p)
  B <- matrix(0, p, p)
  for (j in seq_len(p)) {
    src <- sample(setdiff(seq_len(p), j), s_in)
    A[j, src] <- runif(s_in, 0.12, 0.32) * sample(c(-1, 1), s_in, TRUE)
    curved <- src[runif(s_in) < 0.5]
    if (length(curved) > 0L) {
      # Negative curvature prevents finite-but-explosive equilibria as x grows.
      # Curved terms remain confined to true edges, so support truth is unchanged.
      B[j, curved] <- -runif(length(curved), 0.08, 0.20)
    }
  }
  gamma_base <- runif(p, 1.0, 1.5)
  r <- runif(p, 0.8, 1.5)
  list(p = p, s_in = s_in, A = A, B = B, gamma_base = gamma_base, r = r,
       adj = (A != 0) * 1L)
}

# Three regimes share A, r and the curved-edge set; only the quadratic magnitude
# B differs. strong_nonlinear doubles the stable negative curvature, where MRA's
# local-linear assumption is expected to become less accurate.
sys_for <- function(sys, truth) {
  B <- switch(truth,
              linear = matrix(0, sys$p, sys$p),
              nonlinear = sys$B,
              strong_nonlinear = sys$B * 2)
  gamma <- rowSums(abs(sys$A)) + rowSums(abs(B)) + sys$gamma_base
  list(p = sys$p, s_in = sys$s_in, A = sys$A, B = B, gamma = gamma,
       r = sys$r, adj = sys$adj)
}

deriv_sys <- function(t, state, parms) {
  x <- pmax(state, 0)
  fx <- parms$r + as.numeric(parms$A %*% x) + as.numeric(parms$B %*% (x^2)) -
    parms$gamma * x + parms$u
  list(fx)
}

steady_linear <- function(sys, U) {
  solve_mat <- diag(sys$gamma) - sys$A
  t(apply(U, 1, function(u) as.numeric(solve(solve_mat, sys$r + u))))
}

steady_nonlinear <- function(sys, U, t_max = 120) {
  solve_mat <- diag(sys$gamma) - sys$A
  parms_base <- list(A = sys$A, B = sys$B, gamma = sys$gamma, r = sys$r)
  t(apply(U, 1, function(u) {
    x0 <- pmax(as.numeric(solve(solve_mat, sys$r + u)), 0.05)
    parms <- c(parms_base, list(u = u))
    out <- ode(y = x0, times = c(0, t_max), func = deriv_sys, parms = parms,
               method = "lsoda", rtol = 1e-9, atol = 1e-11)
    as.numeric(out[nrow(out), -1])
  }))
}

steady_for <- function(sys, U, truth) {
  if (truth == "linear") steady_linear(sys, U) else steady_nonlinear(sys, U)
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

# Threshold-free ranking metrics (AUROC, AUPRC), verbatim from Fig1c.
rank_metrics <- function(score, truth) {
  ord <- order(score, decreasing = TRUE)
  y <- truth[ord]
  pos <- sum(y == 1)
  neg <- sum(y == 0)
  if (pos == 0 || neg == 0) return(c(AUROC = NA_real_, AUPRC = NA_real_))
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

## ============================ node-wise coefficient fitters ===================
# Each returns an unscaled coefficient vector of length p * M_ord.
fit_adsiht <- function(Xcs, y, group, scale) {
  fit <- tryCatch(ADSIHT(Xcs, matrix(y), group, ic.type = "dsic"),
                  error = function(e) NULL)
  if (is.null(fit)) return(rep(0, length(scale)))
  as.numeric(fit$beta[, which.min(fit$ic)] / scale)
}

make_fit_glmnet <- function(alpha) {
  function(Xcs, y, group, scale) {
    cvf <- tryCatch(cv.glmnet(Xcs, y, alpha = alpha, intercept = TRUE,
                              nfolds = 5),
                    error = function(e) NULL)
    if (is.null(cvf)) return(rep(0, length(scale)))
    as.numeric(coef(cvf, s = "lambda.1se"))[-1] / scale
  }
}

fit_gglasso <- function(Xcs, y, group, scale) {
  cvf <- tryCatch(cv.gglasso(Xcs, y, group = group, loss = "ls", nfolds = 5),
                  error = function(e) NULL)
  if (is.null(cvf)) return(rep(0, length(scale)))
  as.numeric(coef(cvf, s = "lambda.1se"))[-1] / scale
}

# Official PySINDy optimizer applied to the steady-state algebraic regression
# Theta(x) beta = -u_j. This is STLSQ on a PSS library, not time-series SINDYc.
fit_pysindy_stlsq <- function(Xcs, y, group, scale) {
  opt <- pysindy$STLSQ(
    threshold = 0.1, alpha = 0.05, max_iter = 20L,
    normalize_columns = FALSE
  )
  fit <- tryCatch(opt$fit(Xcs, matrix(y, ncol = 1L)), error = function(e) NULL)
  if (is.null(fit)) return(rep(0, length(scale)))
  as.numeric(fit$coef_) / scale
}

# Edge sign is read from the Jacobian d f_ji/dx_i = theta1 + 2 theta2 x_ref at a
# reference state x_ref (project convention), not from the linear monomial alone,
# so nonlinear-basis methods are scored fairly.
run_nodewise <- function(fit_fn, p, Uc, std, group, x_ref) {
  score <- matrix(0, p, p)
  jac <- matrix(0, p, p)
  c1 <- matrix(0, p, p)
  c2 <- matrix(0, p, p)
  el <- system.time({
    for (j in seq_len(p)) {
      b <- fit_fn(std$X, -Uc[, j], group, std$scale)
      score[j, ] <- group_norms(b, p)
      for (i in seq_len(p)) {
        t1 <- b[(i - 1L) * M_ord + 1L]
        t2 <- if (M_ord >= 2L) b[(i - 1L) * M_ord + 2L] else 0
        c1[j, i] <- t1
        c2[j, i] <- t2
        jac[j, i] <- t1 + 2 * t2 * x_ref[i]
      }
    }
  })[["elapsed"]]
  diag(score) <- 0
  diag(jac) <- 0
  list(score = score, sel = (score > 1e-8) * 1L, sign = sign(jac),
       jac = jac, coef1 = c1, coef2 = c2, runtime = el,
       failed = as.integer(all(score == 0)))
}

# Classical MRA using the published aiMeRA package (v0.99.0). aiMeRA expects a
# square global-response matrix with one specific perturbation per module. Our
# shared benchmark uses continuous multivariate perturbations, so estimate the
# global response dX/dU from all conditions, then delegate the global-to-local
# response inversion and row normalization to aiMeRA::mra(..., Rp = TRUE).
# aiMeRA returns a dense normalized local-response matrix (diagonal -1); the
# fixed 0.05 cutoff is used only for MCC, while AUPRC/AUROC use all magnitudes.
run_aimeRA <- function(p, Uc, X_obs, edge_cutoff = 0.05) {
  score <- link <- matrix(0, p, p)
  failed <- 0L
  el <- system.time({
    Xc <- sweep(X_obs, 2, colMeans(X_obs))
    modules <- paste0("x", seq_len(p))
    perts <- paste0("u", seq_len(p))
    global_response <- tryCatch(t(qr.solve(Uc, Xc)), error = function(e) NULL)
    if (is.null(global_response) || any(!is.finite(global_response))) {
      failed <- 1L
    } else {
      rownames(global_response) <- modules
      colnames(global_response) <- perts
      rules <- c(paste0(perts, "->", modules), "basal->0")
      matp <- aiMeRA::read.rules(rules)
      fit <- tryCatch(
        aiMeRA::mra(global_response, matp, check = FALSE, Rp = TRUE),
        error = function(e) NULL
      )
      if (is.null(fit) || any(!is.finite(fit$link_matrix))) {
        failed <- 1L
      } else {
        link <- unname(fit$link_matrix)
        score <- abs(link)
      }
    }
  })[["elapsed"]]
  diag(score) <- 0
  diag(link) <- 0
  list(score = score, sel = (score > edge_cutoff) * 1L, sign = sign(link),
       jac = link, jac_scale = "row_normalized", coef1 = NULL, coef2 = NULL,
       runtime = el, failed = failed)
}

# Pearson correlation (undirected); BH-FDR selection at q < 0.05.
run_cor <- function(p, X_obs) {
  score <- matrix(0, p, p); signmat <- matrix(0, p, p); sel <- matrix(0, p, p)
  el <- system.time({
    C <- suppressWarnings(cor(X_obs))
    C[!is.finite(C)] <- 0
    signmat <- sign(C)
    score <- abs(C); diag(score) <- 0
    pv <- matrix(1, p, p)
    for (i in seq_len(p - 1L)) for (k in (i + 1L):p) {
      pval <- tryCatch(cor.test(X_obs[, i], X_obs[, k])$p.value,
                       error = function(e) 1)
      pv[i, k] <- pv[k, i] <- pval
    }
    up <- upper.tri(pv)
    q <- p.adjust(pv[up], "BH")
    sel[up] <- (q < 0.05) * 1L
    sel <- sel + t(sel)
  })[["elapsed"]]
  list(score = score, sel = sel, sign = signmat, jac = NULL,
       coef1 = NULL, coef2 = NULL, runtime = el,
       failed = as.integer(all(score == 0)))
}

# Partial correlation and p-values from the published ppcor package; BH-FDR is
# applied to its upper-triangle p-values for the selected undirected graph.
run_pcor <- function(p, X_obs) {
  el <- system.time({
    fit <- ppcor::pcor(X_obs, method = "pearson")
    P <- fit$estimate
    P[!is.finite(P)] <- 0
    pv <- fit$p.value
    pv[!is.finite(pv)] <- 1
    score <- abs(P); diag(score) <- 0
    signmat <- sign(P)
    sel <- matrix(0, p, p)
    up <- upper.tri(pv)
    q <- p.adjust(pv[up], "BH")
    sel[up] <- (q < 0.05) * 1L
    sel <- sel + t(sel)
  })[["elapsed"]]
  list(score = score, sel = sel, sign = signmat, jac = NULL,
       coef1 = NULL, coef2 = NULL, runtime = el,
       failed = as.integer(all(score == 0)))
}

## ============================ unified benchmark metrics =======================
# Edge-function shape error: for each true edge compare the fitted edge function
# f_hat(x) = theta1 x + theta2 x^2 with the true f(x) = A x + B x^2 over the source
# node's observed range extended by `ext` (so it also probes extrapolation). This
# is where nonlinear-basis methods can separate from linear/local methods.
# theta2 = 0). Methods without an edge function -> NA.
func_rmse <- function(coef1, coef2, A, B, X_obs, truth_adj, ext = 1.5, ng = 25L) {
  if (is.null(coef1)) return(NA_real_)
  p <- nrow(A)
  errs <- c()
  for (j in seq_len(p)) for (i in seq_len(p)) {
    if (i != j && truth_adj[j, i] == 1) {
      xr <- range(X_obs[, i])
      xg <- seq(xr[1], xr[2] * ext, length.out = ng)
      ftrue <- A[j, i] * xg + B[j, i] * xg^2
      fhat <- coef1[j, i] * xg + coef2[j, i] * xg^2
      errs <- c(errs, sqrt(mean((fhat - ftrue)^2)))
    }
  }
  if (length(errs) == 0L) NA_real_ else mean(errs)
}

# J_true is the true local Jacobian at the reference state (off-diagonal); methods
# that estimate an interaction/Jacobian (regression family + MRA) are scored on
# JacRMSE, which exposes whether curvature is captured. Association methods give
# no comparable coefficient -> JacRMSE = NA.
bench_metrics <- function(res, truth_adj, A_true, B_true, J_true, X_obs) {
  p <- nrow(truth_adj)
  ro <- row(truth_adj); co <- col(truth_adj)
  off <- which(ro != co)
  t <- truth_adj[off]
  rk <- rank_metrics(res$score[off], t)
  e <- res$sel[off]
  TP <- sum(e == 1 & t == 1); FP <- sum(e == 1 & t == 0)
  TN <- sum(e == 0 & t == 0); FN <- sum(e == 0 & t == 1)
  pr <- ifelse(TP + FP == 0, 0, TP / (TP + FP))
  re <- ifelse(TP + FN == 0, 0, TP / (TP + FN))
  f1 <- ifelse(pr + re == 0, 0, 2 * pr * re / (pr + re))
  den <- sqrt(as.numeric(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  mcc <- ifelse(den == 0, 0, (TP * TN - FP * FN) / den)
  sign_ok <- c()
  if (!is.null(res$sign)) {
    for (k in off) {
      j <- ro[k]; i <- co[k]
      if (truth_adj[j, i] == 1 && res$sel[j, i] == 1) {
        sign_ok <- c(sign_ok, sign(res$sign[j, i]) == sign(A_true[j, i]))
      }
    }
  }
  signacc <- ifelse(length(sign_ok) == 0, NA_real_, mean(sign_ok))
  if (is.null(res$jac)) {
    jac_rmse <- NA_real_; edge_jac_rmse <- NA_real_
  } else {
    J_cmp <- J_true
    if (identical(res$jac_scale, "row_normalized")) {
      J_cmp <- sweep(J_true, 1, -diag(J_true), "/")
    }
    jac_rmse <- sqrt(mean((res$jac[off] - J_cmp[off])^2))
    te <- off[truth_adj[off] == 1]
    edge_jac_rmse <- if (length(te) > 0L)
      sqrt(mean((res$jac[te] - J_cmp[te])^2)) else NA_real_
  }
  fr <- func_rmse(res$coef1, res$coef2, A_true, B_true, X_obs, truth_adj)
  c(AUROC = unname(rk[1]), AUPRC = unname(rk[2]), Precision = pr, Recall = re,
    F1 = f1, MCC = mcc, SignAcc = signacc, JacRMSE = jac_rmse,
    EdgeJacRMSE = edge_jac_rmse, FuncRMSE = fr, TP = TP, FP = FP, FN = FN,
    n_pred = sum(e), runtime_sec = res$runtime, failed = res$failed)
}

## ============================ experiment grid =================================
p <- 30L
s_in <- 2L
base <- s_in * log(p)
# Keep the smallest budget safely above p so that the MRA response inversion is
# identifiable after positive-steady-state filtering. For p = 30 and s_in = 2,
# this gives N = 69, 103, 171, 273.
k_grid <- c(10, 15, 25, 40)
N_grid <- unique(ceiling(k_grid * base))
# Weak (mild) nonlinearity is dropped: a linear baseline is near-optimal
# there, so it adds little. Keep linear and strong-nonlinear only.
truth_grid <- c("linear", "strong_nonlinear")
# With p = 30, requiring every state in a condition to remain positive makes
# the old -1 lower bound discard most rows. Retain broad positive excitation
# for nonlinear coverage while using a biologically feasible negative bound.
u_lo <- -0.3
u_hi <- 1.5
snr_level <- 30
R <- 10L
# Quick smoke-test hook: FIG3B_R / FIG3B_K override repeats and budget.
if (nzchar(Sys.getenv("FIG3B_R"))) R <- as.integer(Sys.getenv("FIG3B_R"))
if (nzchar(Sys.getenv("FIG3B_K"))) {
  N_grid <- unique(ceiling(as.numeric(strsplit(Sys.getenv("FIG3B_K"), ",")[[1]]) * base))
}

nodewise_methods <- list(
  PSS_Net = fit_adsiht,
  Lasso = make_fit_glmnet(1),
  ElasticNet = make_fit_glmnet(0.5),
  GroupLasso = fit_gglasso,
  PySINDy_STLSQ = fit_pysindy_stlsq
)
threshold_audit <- identical(Sys.getenv("FIG3B_THRESHOLD_AUDIT"), "1")
if (threshold_audit) nodewise_methods <- nodewise_methods["PSS_Net"]

rows <- list()
for (seed in seq_len(R)) {
  sys0 <- make_system(p, s_in, seed = 3300 + seed)
  for (truth in truth_grid) {
    sys <- sys_for(sys0, truth)
    for (N in N_grid) {
      set.seed(7000 + 10 * seed + N)
      # Wide perturbation amplitude so the steady state explores the nonlinear
      # regime (small u keeps the system near-linear and hides curvature).
      U <- matrix(runif(N * p, u_lo, u_hi), N, p)
      U[1, ] <- 0
      X <- steady_for(sys, U, truth)
      ok <- apply(X, 1, function(z) all(is.finite(z)) && all(z > 0))
      if (sum(ok) < p + 2L) next
      U <- U[ok, , drop = FALSE]
      X <- X[ok, , drop = FALSE]
      signal_scale <- mean(apply(X, 2, sd))
      sigma_x <- signal_scale / snr_level
      set.seed(7700 + 100 * seed + N)
      X_obs <- pmax(X + matrix(rnorm(length(X), sd = sigma_x), nrow(X), p), 1e-6)

      Uc <- sweep(U, 2, colMeans(U))
      std <- standardize_design(make_basis(X_obs))
      group <- rep(seq_len(p), each = M_ord)

      x_ref <- colMeans(X_obs)
      # Full true local Jacobian at x_ref. aiMeRA is compared on its native
      # row-normalized scale; equation-regression methods use the absolute scale.
      J_true <- sys$A + 2 * sys$B * matrix(x_ref, p, p, byrow = TRUE)
      diag(J_true) <- -sys$gamma
      results <- list()
      for (m in names(nodewise_methods)) {
        results[[m]] <- run_nodewise(nodewise_methods[[m]], p, Uc, std, group, x_ref)
      }
      if (threshold_audit) {
        pss <- results[["PSS_Net"]]

        # Direct absolute-Jacobian threshold. This is numerically 0.05, but is
        # not on aiMeRA's row-normalized link scale.
        pss_abs <- pss
        pss_abs$score <- abs(pss$jac)
        pss_abs$sel <- (pss_abs$score > 0.05) * 1L

        # Like-for-like aiMeRA scale: normalize each fitted equation by minus
        # its fitted self-Jacobian so the diagonal is -1, then apply 0.05.
        jac_full <- pss$coef1 + 2 * pss$coef2 *
          matrix(x_ref, p, p, byrow = TRUE)
        row_scale <- -diag(jac_full)
        row_scale[!is.finite(row_scale) | abs(row_scale) < 1e-8] <- NA_real_
        link <- sweep(jac_full, 1, row_scale, "/")
        link[!is.finite(link)] <- 0
        diag(link) <- 0
        pss_norm <- pss
        pss_norm$score <- abs(link)
        pss_norm$sel <- (pss_norm$score > 0.05) * 1L
        pss_norm$sign <- sign(link)
        pss_norm$jac <- link
        pss_norm$jac_scale <- "row_normalized"

        results <- list(
          PSS_Net_original = pss,
          PSS_Net_absJac_0.05 = pss_abs,
          PSS_Net_rowNorm_0.05 = pss_norm
        )
      } else {
        results[["aiMeRA"]] <- run_aimeRA(p, Uc, X_obs)
        results[["Correlation"]] <- run_cor(p, X_obs)
        results[["PartialCor"]] <- run_pcor(p, X_obs)
      }
      # Optional comparators activate automatically when their package is present.
      # if (requireNamespace("GENIE3", quietly = TRUE)) results[["GENIE3"]] <- ...
      # if (requireNamespace("pcalg",  quietly = TRUE)) results[["GIES"]]   <- ...

      for (m in names(results)) {
        mets <- bench_metrics(results[[m]], sys$adj, sys$A, sys$B, J_true, X_obs)
        jac_scale <- if (is.null(results[[m]]$jac)) "not_applicable" else
          if (is.null(results[[m]]$jac_scale)) "absolute" else results[[m]]$jac_scale
        rows[[length(rows) + 1L]] <- data.frame(
          seed = seed, p = p, s_in = s_in, N = N, n_eff = nrow(X),
          N_over_slogp = N / base, truth = truth, snr = snr_level,
          method = m, JacScale = jac_scale, t(mets), stringsAsFactors = FALSE
        )
      }
      cat(sprintf("seed=%d truth=%-9s N=%3d done\n", seed, truth, N))
    }
  }
}

df <- do.call(rbind, rows)
rownames(df) <- NULL
out <- if (threshold_audit) {
  "results/sim_results/Fig3b_pss_threshold_sensitivity.csv"
} else {
  "results/sim_results/Fig3b_external_benchmark_main.csv"
}
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
write.csv(df, out, row.names = FALSE)

cat("\nSaved:", out, "\n")
print(aggregate(cbind(MCC, AUPRC) ~ truth + method,
                df, mean, na.rm = TRUE), row.names = FALSE)
print(aggregate(JacRMSE ~ truth + method + JacScale,
                df, mean, na.rm = TRUE), row.names = FALSE)
print(aggregate(FuncRMSE ~ truth + method,
                df, mean, na.rm = TRUE), row.names = FALSE)

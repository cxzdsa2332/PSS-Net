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
#          method's documented selection (PSS-Net and aiMeRA both use
#          abs(row-normalized link) > 0.05 in the main comparison);
#          signed local-Jacobian and edge-function errors quantify strength
#          recovery, with linear and nonlinear Jacobian contributions separated.
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
#          (includes M_ord, EdgeJacRMSE/NRMSE, EdgeLinearJacRMSE,
#           EdgeNonlinearJacRMSE, and FuncRMSE/NRMSE)
################################################################################

# Resolve the project root from either the Rscript --file argument or the current
# working directory. This keeps relative inputs/outputs stable when the script is
# launched from RStudio, PowerShell, Terminal, or a directory outside the repo.
find_project_root <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  starts <- getwd()
  if (length(file_arg) > 0L) {
    script_file <- sub("^--file=", "", file_arg[[1L]])
    starts <- c(dirname(path.expand(script_file)), starts)
  }
  source_files <- unlist(lapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) character(0) else frame$ofile
  }), use.names = FALSE)
  if (length(source_files) > 0L) {
    starts <- c(dirname(path.expand(source_files)), starts)
  }

  for (start in unique(starts)) {
    current <- normalizePath(start, winslash = "/", mustWork = FALSE)
    repeat {
      is_root <- file.exists(file.path(
        current, "requirements", "fig3b-pysindy.txt"
      )) && file.exists(file.path(
        current, "sim_script", "03_robustness_benchmarks",
        "Fig3b_external_benchmark_main.R"
      ))
      if (is_root) return(current)
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
  }
  stop(
    "Cannot locate the PSS-Net project root. Run this script from inside the ",
    "repository or invoke it with Rscript path/to/Fig3b_external_benchmark_main.R."
  )
}

project_root <- find_project_root()
if (!identical(normalizePath(getwd(), winslash = "/", mustWork = TRUE),
               project_root)) {
  setwd(project_root)
}

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

# Pin the original-author Python implementation. Python executable discovery is
# deliberately explicit and cross-platform. Priority:
#   1. PSSNET_PYTHON (project-specific override);
#   2. RETICULATE_PYTHON (if it names an existing executable);
#   3. a project virtual environment on Windows or Unix/macOS;
#   4. python3/python on PATH.
pysindy_commit <- "c4421fcec275c8f4cc5c1e93bebb961b212067ae"
python_override <- Sys.getenv("PSSNET_PYTHON", unset = "")
reticulate_override <- Sys.getenv("RETICULATE_PYTHON", unset = "")

if (nzchar(python_override)) {
  python_override <- path.expand(python_override)
  if (!file.exists(python_override)) {
    stop("PSSNET_PYTHON does not exist: ", python_override)
  }
  python_exe <- python_override
} else {
  override_candidates <- if (nzchar(reticulate_override) &&
                             file.exists(path.expand(reticulate_override))) {
    path.expand(reticulate_override)
  } else {
    character(0)
  }
  venv_dirs <- c(".venv-fig3b-standalone", ".venv-fig3b", ".python-fig3b")
  venv_rel <- unlist(lapply(venv_dirs, function(venv) c(
    file.path(venv, "Scripts", "python.exe"), # Windows virtualenv
    file.path(venv, "bin", "python"),         # macOS/Linux virtualenv
    file.path(venv, "bin", "python3")
  )), use.names = FALSE)
  local_candidates <- file.path(project_root, venv_rel)
  path_candidates <- unname(Sys.which(c("python3", "python")))
  candidates <- unique(c(override_candidates, local_candidates, path_candidates))
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (length(candidates) == 0L) {
    stop(
      "No usable Python executable found. Create .venv-fig3b in the project ",
      "root or set PSSNET_PYTHON to its Python executable (bin/python on ",
      "macOS; Scripts/python.exe on Windows)."
    )
  }
  python_exe <- candidates[[1L]]
}
# Canonicalize only the parent directory. On macOS/Linux, venv/bin/python is a
# symlink to the base interpreter; normalizing the full path resolves that link
# and makes reticulate lose the virtual environment's site-packages.
python_exe <- file.path(
  normalizePath(dirname(python_exe), winslash = "/", mustWork = TRUE),
  basename(python_exe)
)
# RETICULATE_PYTHON has higher precedence than use_python(); align it with the
# executable selected above so an inherited setting cannot silently switch OS or
# architecture-specific environments.
Sys.setenv(RETICULATE_PYTHON = python_exe)
reticulate::use_python(python_exe, required = TRUE)
if (!reticulate::py_module_available("pysindy")) {
  requirements_file <- file.path(project_root, "requirements", "fig3b-pysindy.txt")
  stop(paste0(
    "Python package 'pysindy' (audited commit ", pysindy_commit,
    ") is not installed in ", python_exe, ". Install the audited source with:\n  ",
    shQuote(python_exe),
    " -m pip install -r ", shQuote(requirements_file)
  ))
}
pysindy <- reticulate::import("pysindy", delay_load = FALSE)
pysindy_version <- as.character(pysindy$`__version__`)
message("PySINDy backend: ", pysindy_version, " [", reticulate::py_config()$python, "]")

set.seed(303)
# Match the fitted polynomial library [x, x^2] to the maximum order used by the
# simulation truth. Higher-order library robustness is handled separately.
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
# PSS-Net additionally returns the discrete support selected by ADSIHT/DSIC.
# Other fitters return an unscaled coefficient vector of length p * M_ord and
# retain their native nonzero-coefficient selection rule in run_nodewise().
fit_adsiht <- function(Xcs, y, group, scale) {
  fit <- tryCatch(ADSIHT(Xcs, matrix(y), group, ic.type = "dsic"),
                  error = function(e) NULL)
  empty <- list(beta = rep(0, length(scale)), selected_groups = integer(0))
  if (is.null(fit) || length(fit$ic) == 0L) return(empty)
  best <- which.min(fit$ic)
  selected_variables <- fit$A_out[[best]]
  selected_groups <- if (length(selected_variables) == 0L) integer(0) else
    sort(unique(group[selected_variables]))
  list(
    beta = as.numeric(fit$beta[, best] / scale),
    selected_groups = selected_groups
  )
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

# Evaluate the fitted polynomial Jacobian at one reference state. Keeping the
# full diagonal here is necessary for the row-normalized threshold audit.
jacobian_from_coef <- function(coef, x_ref, zero_diagonal = FALSE) {
  p <- dim(coef)[1]
  jac <- matrix(0, p, p)
  powers <- seq_len(dim(coef)[3])
  for (j in seq_len(p)) for (i in seq_len(p)) {
    jac[j, i] <- sum(powers * coef[j, i, ] * x_ref[i]^(powers - 1L))
  }
  if (zero_diagonal) diag(jac) <- 0
  jac
}

# Edge sign is read from the full polynomial Jacobian at a reference state
# x_ref (project convention), not from the linear monomial alone. run_nodewise()
# retains ADSIHT's minimum-DSIC A_out; the main PSS-Net result subsequently maps
# its fitted Jacobian to the row-normalized link used for ranking and selection.
run_nodewise <- function(fit_fn, p, Uc, std, group, x_ref,
                         support_from_fit = FALSE) {
  score <- matrix(0, p, p)
  selected <- matrix(0L, p, p)
  coef <- array(0, dim = c(p, p, M_ord))
  el <- system.time({
    for (j in seq_len(p)) {
      fit_result <- fit_fn(std$X, -Uc[, j], group, std$scale)
      if (support_from_fit) {
        if (!is.list(fit_result) || is.null(fit_result$beta) ||
            is.null(fit_result$selected_groups)) {
          stop("A support-aware fitter must return beta and selected_groups.")
        }
        b <- fit_result$beta
        if (length(fit_result$selected_groups) > 0L) {
          selected[j, fit_result$selected_groups] <- 1L
        }
      } else {
        b <- fit_result
      }
      score[j, ] <- group_norms(b, p)
      for (i in seq_len(p)) {
        cols <- (i - 1L) * M_ord + seq_len(M_ord)
        coef[j, i, ] <- b[cols]
      }
    }
  })[["elapsed"]]
  if (!support_from_fit) selected <- (score > 1e-8) * 1L
  diag(score) <- 0
  diag(selected) <- 0L
  jac <- jacobian_from_coef(coef, x_ref, zero_diagonal = TRUE)
  diag(jac) <- 0
  list(score = score, sel = selected, sign = sign(jac),
       jac = jac, coef = coef, runtime = el,
       selection_rule = if (support_from_fit) "minimum_DSIC_A_out" else
         "native_nonzero_after_fit",
       score_scale = "absolute_group_norm",
       failed = as.integer(all(score == 0)))
}

# Put PSS-Net on the same local-response scale used by aiMeRA for edge ranking
# and binary selection. The fitted absolute-scale Jacobian is retained in
# res$jac for strength recovery; only score/sel/sign use the normalized link.
use_row_normalized_link <- function(res, x_ref, edge_cutoff = 0.05) {
  jac_full <- jacobian_from_coef(res$coef, x_ref)
  row_scale <- -diag(jac_full)
  invalid_scale <- !is.finite(row_scale) | abs(row_scale) < 1e-8
  row_scale[invalid_scale] <- NA_real_
  link <- sweep(jac_full, 1, row_scale, "/")
  link[!is.finite(link)] <- 0
  diag(link) <- 0

  res$score <- abs(link)
  res$sel <- (res$score > edge_cutoff) * 1L
  res$sign <- sign(link)
  res$link <- link
  res$selection_rule <- sprintf("abs(row_normalized_link)>%g", edge_cutoff)
  res$score_scale <- "row_normalized_link"
  res$invalid_row_scale <- sum(invalid_scale)
  res
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
       jac = link, jac_scale = "row_normalized", coef = NULL,
       selection_rule = sprintf("abs(row_normalized_link)>%g", edge_cutoff),
       score_scale = "row_normalized_link",
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
       coef = NULL, runtime = el,
       selection_rule = "BH_FDR<0.05", score_scale = "absolute_correlation",
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
       coef = NULL, runtime = el,
       selection_rule = "BH_FDR<0.05",
       score_scale = "absolute_partial_correlation",
       failed = as.integer(all(score == 0)))
}

## ============================ unified benchmark metrics =======================
# Edge-function shape error: for each true edge compare the full fitted
# polynomial with f_true(x) = A x + B x^2 over the observed source range,
# extended by `ext` to probe modest extrapolation. FuncNRMSE divides each edge's
# RMSE by its true RMS function magnitude before averaging across edges.
edge_function_metrics <- function(coef, A, B, X_obs, truth_adj,
                                  ext = 1.5, ng = 25L) {
  if (is.null(coef)) {
    return(c(FuncRMSE = NA_real_, FuncNRMSE = NA_real_))
  }
  p <- nrow(A)
  powers <- seq_len(dim(coef)[3])
  rmse <- nrmse <- c()
  for (j in seq_len(p)) for (i in seq_len(p)) {
    if (i != j && truth_adj[j, i] == 1) {
      xr <- range(X_obs[, i])
      xg <- seq(xr[1], xr[2] * ext, length.out = ng)
      ftrue <- A[j, i] * xg + B[j, i] * xg^2
      fhat <- vapply(xg, function(x) {
        sum(coef[j, i, ] * x^powers)
      }, numeric(1))
      edge_rmse <- sqrt(mean((fhat - ftrue)^2))
      true_rms <- sqrt(mean(ftrue^2))
      rmse <- c(rmse, edge_rmse)
      nrmse <- c(nrmse, edge_rmse / pmax(true_rms, 1e-12))
    }
  }
  if (length(rmse) == 0L) {
    c(FuncRMSE = NA_real_, FuncNRMSE = NA_real_)
  } else {
    c(FuncRMSE = mean(rmse), FuncNRMSE = mean(nrmse))
  }
}

# Decompose local edge strength at x_ref into the linear derivative theta_1 and
# the nonlinear derivative sum_{m>=2} m theta_m x_ref^(m-1). Both components are
# on the Jacobian scale, so their RMSE values can be compared and add up to the
# total fitted local effect. The simulation truth has quadratic B and zero cubic
# coefficient; a fitted cubic contribution is therefore counted as error.
edge_strength_metrics <- function(coef, A, B, x_ref, truth_adj) {
  empty <- c(EdgeLinearJacRMSE = NA_real_,
             EdgeNonlinearJacRMSE = NA_real_)
  if (is.null(coef)) return(empty)
  p <- nrow(A)
  te <- which(row(truth_adj) != col(truth_adj) & truth_adj == 1)
  if (length(te) == 0L) return(empty)

  linear_est <- coef[, , 1L]
  nonlinear_est <- matrix(0, p, p)
  if (dim(coef)[3] >= 2L) {
    for (m in 2:dim(coef)[3]) {
      nonlinear_est <- nonlinear_est +
        m * coef[, , m] * matrix(x_ref^(m - 1L), p, p, byrow = TRUE)
    }
  }
  nonlinear_true <- 2 * B * matrix(x_ref, p, p, byrow = TRUE)
  c(
    EdgeLinearJacRMSE = sqrt(mean((linear_est[te] - A[te])^2)),
    EdgeNonlinearJacRMSE = sqrt(mean(
      (nonlinear_est[te] - nonlinear_true[te])^2
    ))
  )
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
        sign_ok <- c(sign_ok, sign(res$sign[j, i]) == sign(J_true[j, i]))
      }
    }
  }
  signacc <- ifelse(length(sign_ok) == 0, NA_real_, mean(sign_ok))
  if (is.null(res$jac)) {
    jac_rmse <- NA_real_; edge_jac_rmse <- NA_real_
    edge_jac_nrmse <- NA_real_
  } else {
    J_cmp <- J_true
    if (identical(res$jac_scale, "row_normalized")) {
      J_cmp <- sweep(J_true, 1, -diag(J_true), "/")
    }
    jac_rmse <- sqrt(mean((res$jac[off] - J_cmp[off])^2))
    te <- off[truth_adj[off] == 1]
    if (length(te) > 0L) {
      edge_jac_rmse <- sqrt(mean((res$jac[te] - J_cmp[te])^2))
      edge_jac_nrmse <- edge_jac_rmse /
        pmax(sqrt(mean(J_cmp[te]^2)), 1e-12)
    } else {
      edge_jac_rmse <- edge_jac_nrmse <- NA_real_
    }
  }
  strength <- edge_strength_metrics(
    res$coef, A_true, B_true, colMeans(X_obs), truth_adj
  )
  function_error <- edge_function_metrics(
    res$coef, A_true, B_true, X_obs, truth_adj
  )
  c(AUROC = unname(rk[1]), AUPRC = unname(rk[2]), Precision = pr, Recall = re,
    F1 = f1, MCC = mcc, SignAcc = signacc, JacRMSE = jac_rmse,
    EdgeJacRMSE = edge_jac_rmse, EdgeJacNRMSE = edge_jac_nrmse,
    strength, function_error, TP = TP, FP = FP, FN = FN,
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
R <- 30L
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
cell_keys <- as.vector(outer(truth_grid, N_grid, paste, sep = "|"))
cell_counts <- setNames(integer(length(cell_keys)), cell_keys)
max_attempts <- max(20L * R, R + 100L)
attempt <- 0L
while (any(cell_counts < R)) {
  attempt <- attempt + 1L
  if (attempt > max_attempts) {
    stop("Unable to collect ", R, " valid replicates for every Fig3b cell after ",
         max_attempts, " simulation seeds.")
  }
  sys0 <- make_system(p, s_in, seed = 3300 + attempt)
  for (truth in truth_grid) {
    sys <- sys_for(sys0, truth)
    for (N in N_grid) {
      cell_key <- paste(truth, N, sep = "|")
      if (cell_counts[[cell_key]] >= R) next
      set.seed(7000 + 10 * attempt + N)
      # Wide perturbation amplitude so the steady state explores the nonlinear
      # regime (small u keeps the system near-linear and hides curvature).
      U <- matrix(runif(N * p, u_lo, u_hi), N, p)
      U[1, ] <- 0
      X <- steady_for(sys, U, truth)
      ok <- apply(X, 1, function(z) all(is.finite(z)) && all(z > 0))
      if (sum(ok) < p + 2L) {
        cat(sprintf("sim_seed=%d truth=%-16s N=%3d skipped (n_eff=%d)\n",
                    attempt, truth, N, sum(ok)))
        next
      }
      rep_id <- cell_counts[[cell_key]] + 1L
      U <- U[ok, , drop = FALSE]
      X <- X[ok, , drop = FALSE]
      signal_scale <- mean(apply(X, 2, sd))
      sigma_x <- signal_scale / snr_level
      set.seed(7700 + 100 * attempt + N)
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
        results[[m]] <- run_nodewise(
          nodewise_methods[[m]], p, Uc, std, group, x_ref,
          support_from_fit = identical(m, "PSS_Net")
        )
      }
      if (threshold_audit) {
        pss <- results[["PSS_Net"]]
        pss$selection_rule <- "minimum_DSIC_A_out"
        pss$score_scale <- "absolute_group_norm"

        # Direct absolute-Jacobian threshold. This is numerically 0.05, but is
        # not on aiMeRA's row-normalized link scale.
        pss_abs <- pss
        pss_abs$score <- abs(pss$jac)
        pss_abs$sel <- (pss_abs$score > 0.05) * 1L
        pss_abs$selection_rule <- "abs(absolute_Jacobian)>0.05"
        pss_abs$score_scale <- "absolute_Jacobian"

        # Like-for-like aiMeRA scale: normalize each fitted equation by minus
        # its fitted self-Jacobian so the diagonal is -1, then apply 0.05.
        pss_norm <- use_row_normalized_link(pss, x_ref, edge_cutoff = 0.05)
        pss_norm$jac <- pss_norm$link
        pss_norm$jac_scale <- "row_normalized"

        results <- list(
          PSS_Net_original = pss,
          PSS_Net_absJac_0.05 = pss_abs,
          PSS_Net_rowNorm_0.05 = pss_norm
        )
      } else {
        # Main PSS-Net result: same row-normalized score and 0.05 cutoff as MRA.
        # Absolute coefficients/Jacobian remain attached for recovery metrics.
        results[["PSS_Net"]] <- use_row_normalized_link(
          results[["PSS_Net"]], x_ref, edge_cutoff = 0.05
        )
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
          seed = rep_id, sim_seed = attempt, p = p, s_in = s_in, M_ord = M_ord,
          N = N, n_eff = nrow(X),
          N_over_slogp = N / base, truth = truth, snr = snr_level,
          method = m, JacScale = jac_scale,
          ScoreScale = results[[m]]$score_scale,
          SelectionRule = results[[m]]$selection_rule,
          t(mets), stringsAsFactors = FALSE
        )
      }
      cell_counts[[cell_key]] <- rep_id
      cat(sprintf("rep=%d sim_seed=%d truth=%-16s N=%3d done\n",
                  rep_id, attempt, truth, N))
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
print(aggregate(cbind(JacRMSE, EdgeJacRMSE, EdgeJacNRMSE) ~
                  truth + method + JacScale,
                df, mean, na.rm = TRUE), row.names = FALSE)
coef_df <- df[is.finite(df$EdgeLinearJacRMSE) &
                is.finite(df$EdgeNonlinearJacRMSE), ]
print(aggregate(cbind(EdgeLinearJacRMSE, EdgeNonlinearJacRMSE) ~
                  truth + method,
                coef_df, mean, na.rm = TRUE), row.names = FALSE)
func_df <- df[is.finite(df$FuncRMSE) & is.finite(df$FuncNRMSE), ]
print(aggregate(cbind(FuncRMSE, FuncNRMSE) ~ truth + method,
                func_df, mean, na.rm = TRUE), row.names = FALSE)

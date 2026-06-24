rm(list = ls())

################################################################################
# Fig1x_basis_misspecification.R -- basis/library misspecification sensitivity
#
# Foundation/identifiability companion to Figure 1: it decomposes what PSS-Net
# can recover as the FITTED basis library is varied away from the true nonlinear
# edge functions -- an identifiability statement (edge SUPPORT and SIGN recover
# across libraries; the edge FUNCTION shape needs an adequate dictionary), not
# merely a stress test. Feeds the compact Fig1g basis-robustness panel in Fig1.R.
#
# Purpose: quantify how robust PSS-Net is when the FITTED basis library does not
#          match the true nonlinear edge functions. The true data-generating
#          interaction shape is crossed with several fitted libraries, so the
#          diagonal of the grid is "matched" and the off-diagonal is misspecified.
#
#          True edge function f_ji(x) = A_ji * g(x) (+ B_ji x^2 for poly2), all
#          sharing the same support and small-signal slope A_ji:
#            poly2  : A x + B x^2        (polynomial-representable; control)
#            monod  : A * x / (1 + x)    (saturating, non-polynomial)
#            sine   : A * sin(x)         (oscillatory, non-polynomial)
#
#          Fitted no-intercept libraries (phi_m(0) = 0), node-wise ADSIHT:
#            linear : [x]
#            poly2  : [x, x^2]
#            poly3  : [x, x^2, x^3]
#            monod  : [x, x/(1+x)]       (matches monod truth)
#            fourier: [x, sin x, cos x - 1]  (matches sine truth)
#
#          Metrics, per (truth, library, seed): edge support MCC / precision /
#          recall, signed-Jacobian accuracy, and edge-function NRMSE over the
#          observed source range. The headline question is whether SUPPORT
#          recovery (MCC) is robust to library choice while FUNCTION recovery
#          (FuncNRMSE) needs an adequate library.
#
#          Reuse: ER system + deSolve steady states follow Fig3b /
#          pss_net_scalefree.R; the configurable basis follows program/
#          pss_net_demo.R. Single generous sample budget so the library -- not
#          the sample size -- is the limiting factor. Run from the repo root.
#
# Input:   none
# Output:  results/sim_results/Fig1x_basis_misspecification.csv
################################################################################

suppressMessages({
  library(deSolve)
  library(ADSIHT)
})

set.seed(303)

## ============================ basis library ==================================
# Each library is a list of no-intercept functions phi_m and their derivatives.
make_basis <- function(type) {
  switch(type,
    linear = list(
      funcs = list(function(x) x),
      dfuncs = list(function(x) rep(1, length(x))), M = 1L),
    poly2 = list(
      funcs = list(function(x) x, function(x) x^2),
      dfuncs = list(function(x) rep(1, length(x)), function(x) 2 * x), M = 2L),
    poly3 = list(
      funcs = list(function(x) x, function(x) x^2, function(x) x^3),
      dfuncs = list(function(x) rep(1, length(x)), function(x) 2 * x,
                    function(x) 3 * x^2), M = 3L),
    monod = list(
      funcs = list(function(x) x, function(x) x / (1 + x)),
      dfuncs = list(function(x) rep(1, length(x)),
                    function(x) 1 / (1 + x)^2), M = 2L),
    fourier = list(
      funcs = list(function(x) x, function(x) sin(x), function(x) cos(x) - 1),
      dfuncs = list(function(x) rep(1, length(x)), function(x) cos(x),
                    function(x) -sin(x)), M = 3L),
    stop("unknown basis: ", type))
}

build_design <- function(X, basis) {
  p <- ncol(X); M <- basis$M
  Psi <- matrix(0, nrow(X), p * M)
  for (i in seq_len(p)) for (m in seq_len(M)) {
    Psi[, (i - 1L) * M + m] <- basis$funcs[[m]](X[, i])
  }
  Psi
}

## ============================ true system + dynamics =========================
# ER topology: each target gets s_in random sources. A_ji is the shared
# small-signal slope; B_ji (poly2 only) adds stable negative curvature on half
# the incoming edges. gamma keeps the steady state contractive.
make_system <- function(p, s_in, seed) {
  set.seed(seed)
  A <- matrix(0, p, p)
  B <- matrix(0, p, p)
  for (j in seq_len(p)) {
    src <- sample(setdiff(seq_len(p), j), s_in)
    A[j, src] <- runif(s_in, 0.3, 0.6) * sample(c(-1, 1), s_in, TRUE)
    curved <- src[runif(s_in) < 0.5]
    if (length(curved) > 0L) B[j, curved] <- -runif(length(curved), 0.1, 0.25)
  }
  gamma <- rowSums(abs(A)) + rowSums(abs(B)) + runif(p, 1.0, 1.5)
  r <- runif(p, 0.8, 1.5)
  list(p = p, s_in = s_in, A = A, B = B, gamma = gamma, r = r,
       adj = (A != 0) * 1L)
}

# Per-edge true function and its derivative (used by both the ODE and metrics).
true_f <- function(truth, a, b, x) {
  xp <- pmax(x, 0)
  switch(truth,
    poly2 = a * xp + b * xp^2,
    monod = a * xp / (1 + xp),
    sine  = a * sin(xp))
}
true_df <- function(truth, a, b, x) {
  switch(truth,
    poly2 = a + 2 * b * x,
    monod = a / (1 + x)^2,
    sine  = a * cos(x))
}

deriv_sys <- function(t, state, parms) {
  x <- pmax(state, 0)
  inter <- switch(parms$truth,
    poly2 = as.numeric(parms$A %*% x) + as.numeric(parms$B %*% (x^2)),
    monod = as.numeric(parms$A %*% (x / (1 + x))),
    sine  = as.numeric(parms$A %*% sin(x)))
  list(parms$r + inter - parms$gamma * x + parms$u)
}

steady_many <- function(sys, U, truth, t_max = 120) {
  solve_mat <- diag(sys$gamma) - sys$A
  parms_base <- list(A = sys$A, B = sys$B, gamma = sys$gamma, r = sys$r,
                     truth = truth)
  t(apply(U, 1, function(u) {
    x0 <- pmax(as.numeric(solve(solve_mat, sys$r + u)), 0.05)
    out <- tryCatch(
      ode(y = x0, times = c(0, t_max), func = deriv_sys,
          parms = c(parms_base, list(u = u)), method = "lsoda",
          rtol = 1e-9, atol = 1e-11),
      error = function(e) NULL)
    if (is.null(out)) rep(NA_real_, sys$p) else as.numeric(out[nrow(out), -1])
  }))
}

## ============================ node-wise inference ============================
standardize <- function(Psi) {
  ctr <- colMeans(Psi)
  Pc <- sweep(Psi, 2, ctr)
  sdv <- pmax(apply(Pc, 2, sd), 1e-10)
  list(X = sweep(Pc, 2, sdv, "/"), scale = sdv)
}

fit_pss <- function(U, X, basis) {
  p <- ncol(X); M <- basis$M
  std <- standardize(build_design(X, basis))
  group <- rep(seq_len(p), each = M)
  rhs <- sapply(seq_len(p), function(j) -(U[, j] - mean(U[, j])))
  theta <- matrix(0, p * M, p)
  for (j in seq_len(p)) {
    fit <- tryCatch(ADSIHT(std$X, matrix(rhs[, j]), group, ic.type = "dsic"),
                    error = function(e) NULL)
    if (is.null(fit) || length(fit$ic) == 0L) next
    theta[, j] <- fit$beta[, which.min(fit$ic)] / std$scale
  }
  theta
}

## ============================ metrics ========================================
# Support MCC / precision / recall, signed-Jacobian accuracy at the reference
# state, and edge-function NRMSE over each true edge's observed source range.
eval_fit <- function(theta, basis, sys, truth, X_obs, ext = 1.2, ng = 25L) {
  p <- sys$p; M <- basis$M
  x_ref <- colMeans(X_obs)
  dref <- sapply(seq_len(M), function(m) basis$dfuncs[[m]](x_ref))  # p x M
  adj <- matrix(0L, p, p); jac <- matrix(0, p, p)
  for (j in seq_len(p)) for (i in seq_len(p)) {
    th <- theta[(i - 1L) * M + seq_len(M), j]
    if (sqrt(sum(th^2)) >= 1e-8) adj[j, i] <- 1L
    jac[j, i] <- sum(th * dref[i, ])
  }
  diag(adj) <- 0L; diag(jac) <- 0

  off <- which(row(sys$adj) != col(sys$adj))
  e <- adj[off]; t <- sys$adj[off]
  TP <- sum(e == 1 & t == 1); FP <- sum(e == 1 & t == 0)
  TN <- sum(e == 0 & t == 0); FN <- sum(e == 0 & t == 1)
  pr <- ifelse(TP + FP == 0, 0, TP / (TP + FP))
  re <- ifelse(TP + FN == 0, 0, TP / (TP + FN))
  den <- sqrt(as.numeric(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  mcc <- ifelse(den == 0, 0, (TP * TN - FP * FN) / den)

  # signed Jacobian accuracy on detected true edges
  te <- which(sys$adj == 1L & row(sys$adj) != col(sys$adj), arr.ind = TRUE)
  sgn <- c(); ferr <- c()
  for (k in seq_len(nrow(te))) {
    j <- te[k, 1]; i <- te[k, 2]
    if (adj[j, i] == 1L) {
      true_d <- true_df(truth, sys$A[j, i], sys$B[j, i], x_ref[i])
      sgn <- c(sgn, sign(jac[j, i]) == sign(true_d))
    }
    xr <- range(X_obs[, i])
    xg <- seq(xr[1], xr[2] * ext, length.out = ng)
    fhat <- rowSums(sapply(seq_len(M),
      function(m) theta[(i - 1L) * M + m, j] * basis$funcs[[m]](xg)))
    ftrue <- true_f(truth, sys$A[j, i], sys$B[j, i], xg)
    ferr <- c(ferr, sqrt(mean((fhat - ftrue)^2)) /
                pmax(sqrt(mean(ftrue^2)), 1e-12))
  }
  c(MCC = mcc, Precision = pr, Recall = re,
    SignAcc = ifelse(length(sgn) == 0, NA_real_, mean(sgn)),
    FuncNRMSE = ifelse(length(ferr) == 0, NA_real_, mean(ferr)))
}

## ============================ experiment grid ================================
p <- 20L
s_in <- 2L
base <- s_in * log(p)
N <- ceiling(20 * base)            # one generous budget: isolate library effect
snr_level <- 30
R <- 8L
u_lo <- -0.3
u_hi <- 1.5
truth_grid <- c("poly2", "monod", "sine")
lib_grid <- c("linear", "poly2", "poly3", "monod", "fourier")
# matched (truth, library) pairs, for tagging the grid diagonal
matched <- list(poly2 = "poly2", monod = "monod", sine = "fourier")
if (nzchar(Sys.getenv("FIG1X_R"))) R <- as.integer(Sys.getenv("FIG1X_R"))

rows <- list()
for (seed in seq_len(R)) {
  sys <- make_system(p, s_in, seed = 5300 + seed)
  for (truth in truth_grid) {
    set.seed(8000 + 10 * seed)
    U <- matrix(runif(N * p, u_lo, u_hi), N, p)
    U[1, ] <- 0
    X <- steady_many(sys, U, truth)
    ok <- apply(X, 1, function(z) all(is.finite(z)) && all(z > 0))
    if (sum(ok) < p + 2L) next
    U <- U[ok, , drop = FALSE]; X <- X[ok, , drop = FALSE]
    sigma <- mean(apply(X, 2, sd)) / snr_level
    set.seed(8800 + 10 * seed)
    X_obs <- pmax(X + matrix(rnorm(length(X), sd = sigma), nrow(X), p), 1e-6)
    for (lib in lib_grid) {
      basis <- make_basis(lib)
      theta <- fit_pss(U, X_obs, basis)
      m <- eval_fit(theta, basis, sys, truth, X_obs)
      rows[[length(rows) + 1L]] <- data.frame(
        seed = seed, p = p, N = nrow(X), snr = snr_level,
        truth = truth, library = lib,
        matched = identical(lib, matched[[truth]]),
        t(m), stringsAsFactors = FALSE)
    }
    cat(sprintf("seed=%d truth=%-5s done (n_eff=%d)\n", seed, truth, nrow(X)))
  }
}

df <- do.call(rbind, rows)
rownames(df) <- NULL
out <- "results/sim_results/Fig1x_basis_misspecification.csv"
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
write.csv(df, out, row.names = FALSE)
cat("\nSaved:", out, "\n\n")

# ---- summaries ----
fmt <- function(v) formatC(v, format = "f", digits = 3)
agg_mcc <- aggregate(MCC ~ truth + library, df, mean, na.rm = TRUE)
agg_fun <- aggregate(FuncNRMSE ~ truth + library, df, mean, na.rm = TRUE)
agg_sgn <- aggregate(SignAcc ~ truth + library, df, mean, na.rm = TRUE)
cat("== support MCC (truth x library; higher better) ==\n")
print(reshape(agg_mcc, idvar = "truth", timevar = "library",
              direction = "wide"), row.names = FALSE)
cat("\n== edge-function NRMSE (truth x library; lower better) ==\n")
print(reshape(agg_fun, idvar = "truth", timevar = "library",
              direction = "wide"), row.names = FALSE)
cat("\n== signed-Jacobian accuracy (truth x library) ==\n")
print(reshape(agg_sgn, idvar = "truth", timevar = "library",
              direction = "wide"), row.names = FALSE)

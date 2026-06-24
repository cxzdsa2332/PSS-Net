rm(list = ls())

################################################################################
# pss_net_demo.R -- standalone PSS-Net inference demo (first version)
#
# A self-contained reference implementation of PSS-Net network inference:
#   1. simulate a standard sparse nonlinear additive-ODE system whose edges
#      roughly follow a scale-free (hub) topology, and its perturbed
#      steady-state (PSS) data (x*, u);
#   2. run ADSIHT inference (double sparsity, no-intercept SINDy-style basis
#      library -- monomial by default, polynomial/fourier presets, or a fully
#      custom function library -- centred + scaled design, Jacobian for edge
#      signs). Default = the joint block-diagonal solve over the whole network;
#      node-wise is optional.
#   3. report recovery metrics and draw the true vs inferred directed signed
#      network on the active graphics device (no files are written).
#
# Unlike the formal sim_script/ entries, this program mixes simulation,
# inference and plotting so it can be run end-to-end as a demonstration.
#
# Run interactively (e.g. in RStudio) to see the network, or:
#   Rscript program/pss_net_demo.R
################################################################################

suppressMessages({
  library(deSolve)
  library(ADSIHT)
  library(igraph)
})

## ============================ 0. basis library (SINDy-style) =================
# A "basis" is a library of univariate functions phi_m applied to each source
# node's steady state x_i. The additive model uses
#   f_ji(x_i) = sum_m theta_jim * phi_m(x_i),
# i.e. each source enters through the SAME library (SINDy-style feature
# dictionary). Because every source shares one library of size M, the ADSIHT
# group vector is DERIVED automatically -- one group per source, each of size M
# -- so a custom library never needs a hand-built group.
#
# Constraints:
#   * every phi_m must satisfy phi_m(0) = 0 (no intercept) so the steady-state
#     identity f_ji(0) = 0 holds (centred + scaled design recovers the rest);
#   * supplying analytic derivatives dphi_m (dfuncs) gives exact Jacobian edge
#     signs; otherwise a central finite-difference derivative is used.
#
# Presets:
#   "monomial"/"polynomial": x, x^2, ..., x^order   (default order = 2)
#   "fourier"             : x, sin(k x), cos(k x)-1 for k = 1..order
#   "custom"              : user supplies funcs (+ optional dfuncs, labels)
make_basis <- function(type = c("monomial", "polynomial", "fourier", "custom"),
                       order = 2L, funcs = NULL, dfuncs = NULL, labels = NULL) {
  type <- match.arg(type)
  if (type == "monomial" || type == "polynomial") {
    pw <- seq_len(order)
    funcs  <- lapply(pw, function(m) { force(m); function(x) x^m })
    dfuncs <- lapply(pw, function(m) { force(m); function(x) m * x^(m - 1) })
    labels <- ifelse(pw == 1L, "x", paste0("x^", pw))
  } else if (type == "fourier") {
    ks <- seq_len(order)
    sin_f  <- lapply(ks, function(k) { force(k); function(x) sin(k * x) })
    sin_df <- lapply(ks, function(k) { force(k); function(x)  k * cos(k * x) })
    cos_f  <- lapply(ks, function(k) { force(k); function(x) cos(k * x) - 1 })
    cos_df <- lapply(ks, function(k) { force(k); function(x) -k * sin(k * x) })
    funcs  <- c(list(function(x) x), sin_f, cos_f)
    dfuncs <- c(list(function(x) rep(1, length(x))), sin_df, cos_df)
    labels <- c("x", paste0("sin(", ks, "x)"), paste0("cos(", ks, "x)-1"))
  } else {  # custom
    if (is.null(funcs) || length(funcs) == 0L)
      stop("custom basis needs a non-empty 'funcs' list of phi_m(x).")
  }
  M <- length(funcs)
  if (is.null(labels)) labels <- paste0("phi", seq_len(M))
  if (is.null(dfuncs)) {  # central finite-difference fallback
    dfuncs <- lapply(funcs, function(f) {
      force(f)
      function(x) { h <- 1e-6; (f(x + h) - f(x - h)) / (2 * h) }
    })
  }
  list(funcs = funcs, dfuncs = dfuncs, labels = labels, M = M, type = type)
}

# Evaluate the library on every source column: column block i (of width M) holds
# phi_1(x_i), ..., phi_M(x_i). Returns an N x (p*M) no-intercept design.
build_design <- function(X, basis) {
  p <- ncol(X); M <- basis$M
  Psi <- matrix(0, nrow(X), p * M)
  for (i in seq_len(p)) for (m in seq_len(M)) {
    Psi[, (i - 1L) * M + m] <- basis$funcs[[m]](X[, i])
  }
  Psi
}

## ============================ 1. standard nonlinear PSS simulation ===========
# Additive ODE:  dx_j/dt = r_j + sum_{i!=j} (A_ji x_i + B_ji x_i^2)
#                          - gamma_j x_j + u_j
# Topology (default "scalefree"): preferential attachment gives a power-law
# out-degree, so a few hub regulators emerge; "random" keeps a fixed in-degree.
# Half of each target's incoming edges carry a quadratic term (B != 0), so the
# system exercises within-group (which-monomial) sparsity. Diagonal dominance
# (gamma_j > sum_i |A_ji| + |B_ji|) keeps the steady state stable.
simulate_pss_nonlinear <- function(p = 10L, topology = c("scalefree", "random"),
                                   avg_in = 2, n_in = 2L, N = 200L, snr = 30,
                                   seed = 1L, u_lo = -0.6, u_hi = 0.9) {
  topology <- match.arg(topology)
  set.seed(seed)
  A <- matrix(0, p, p)
  B <- matrix(0, p, p)

  add_curvature <- function(j, src) {
    if (length(src) == 0L) return(invisible())
    curved <- src[seq_len(max(1L, floor(length(src) / 2L)))]
    B[j, curved] <<- runif(length(curved), 0.15, 0.35) * sign(A[j, curved])
  }

  if (topology == "scalefree") {
    outdeg <- rep(1, p)                                  # +1 smoothing weight
    for (j in seq_len(p)) {
      k_j <- min(1L + rpois(1, max(avg_in - 1, 0)), p - 1L)  # variable in-degree
      chosen <- integer(0)
      for (t in seq_len(k_j)) {
        w <- outdeg; w[c(j, chosen)] <- 0               # no self-loop / repeat
        i <- sample.int(p, 1, prob = w)                 # P(i) prop. to out-degree
        chosen <- c(chosen, i)
        A[j, i] <- runif(1, 0.3, 0.7) * sample(c(-1, 1), 1)
        outdeg[i] <- outdeg[i] + 1
      }
      add_curvature(j, chosen)
    }
  } else {
    for (j in seq_len(p)) {
      src <- sample(setdiff(seq_len(p), j), n_in)
      A[j, src] <- runif(n_in, 0.3, 0.7) * sample(c(-1, 1), n_in, TRUE)
      add_curvature(j, src)
    }
  }

  gamma <- rowSums(abs(A)) + rowSums(abs(B)) + runif(p, 1.0, 1.5)
  r <- runif(p, 0.8, 1.5)
  sys <- list(p = p, topology = topology, A = A, B = B, gamma = gamma, r = r,
              adj = (A != 0) * 1L, outdeg_true = colSums(A != 0))
  sys$x_wt <- steady_one(sys, rep(0, p), x0 = r / gamma)

  U <- matrix(runif(N * p, u_lo, u_hi), N, p)
  U[1, ] <- 0  # first condition is the unperturbed wild type
  X <- steady_many(sys, U)
  ok <- apply(X, 1, function(z) all(is.finite(z)) && all(z > 0))
  U <- U[ok, , drop = FALSE]
  X <- X[ok, , drop = FALSE]

  signal_scale <- mean(apply(X, 2, sd))
  sigma <- signal_scale / snr
  set.seed(seed + 99L)
  X_obs <- pmax(X + matrix(rnorm(length(X), sd = sigma), nrow(X), p), 1e-6)

  c(sys, list(U = U, X = X, X_obs = X_obs, sigma = sigma, snr = snr,
              n_eff = nrow(X)))
}

# Steady state of one perturbation condition by integrating the ODE to t_max.
steady_one <- function(sys, u, x0 = NULL, t_max = 200) {
  p <- sys$p
  if (is.null(x0)) x0 <- sys$x_wt
  deriv <- function(t, x, parms) {
    xp <- pmax(x, 0)
    inter <- as.numeric(sys$A %*% xp) + as.numeric(sys$B %*% (xp^2))
    list(sys$r + inter - sys$gamma * xp + u)
  }
  out <- tryCatch(
    ode(y = x0, times = c(0, t_max), func = deriv, parms = NULL,
        method = "lsoda", rtol = 1e-9, atol = 1e-11),
    error = function(e) NULL)
  if (is.null(out)) return(rep(NA_real_, p))
  as.numeric(out[nrow(out), -1])
}

steady_many <- function(sys, U) {
  t(apply(U, 1, function(u) steady_one(sys, u)))
}

## ============================ 2. PSS-Net inference ===========================
# Centred + scaled no-intercept design built from a SINDy-style basis library
# (see make_basis). For target j the steady-state equation gives
# -u_j = sum_i f_ji(x_i) + (self term), so we regress the centred -u_j on
# Psi(X). Group = source node (double sparsity: group selects sources,
# within-group selects which library terms phi_m are active). The group vector
# is derived from the library size M, so any custom basis works unchanged. Edge
# signs come from the Jacobian d f_ji/dx_i = sum_m theta_jim phi_m'(x_ref_i).
#
# method = "joint" (default): one ADSIHT solve over the block-diagonal design
#   X = I_p (x) Psi_cs with group = rep(1:(p*p), each = M) (CLAUDE.md rule); the
#   global DSIC criterion pools model complexity across all targets, which helps
#   on heterogeneous (hub) networks. The dense I_p (x) Psi is O(p^2 M) wide, so
#   it is slow / memory-heavy for large p (warned at p > 100).
# method = "nodewise": p independent ADSIHT solves; cheaper, the default for
#   large or homogeneous networks.
pss_net <- function(U, X, basis = make_basis(), method = c("joint", "nodewise"),
                    x_ref = NULL, edge_tol = 1e-8) {
  method <- match.arg(method)
  p <- ncol(X); M <- basis$M
  if (is.null(x_ref)) x_ref <- colMeans(X)
  if (method == "joint" && p > 100L) {
    warning(sprintf(paste0("joint block-diagonal design is O(p^2 M) wide ",
            "(here %d x %d); p > 100 may be slow or memory-heavy -- ",
            "consider method = 'nodewise'."), nrow(X) * p, p * p * M))
  }

  Psi <- build_design(X, basis)
  Psi_c <- sweep(Psi, 2, colMeans(Psi))
  sdv <- pmax(apply(Psi_c, 2, sd), 1e-10)
  Psi_cs <- sweep(Psi_c, 2, sdv, "/")
  rhs <- sapply(seq_len(p), function(j) -(U[, j] - mean(U[, j])))  # N x p

  theta <- matrix(0, p * M, p)  # column j: node j's p*M coefficients
  failed <- 0L
  if (method == "nodewise") {
    group <- rep(seq_len(p), each = M)
    for (j in seq_len(p)) {
      fit <- tryCatch(ADSIHT(Psi_cs, matrix(rhs[, j]), group, ic.type = "dsic"),
                      error = function(e) NULL)
      if (is.null(fit) || length(fit$ic) == 0L) { failed <- failed + 1L; next }
      theta[, j] <- fit$beta[, which.min(fit$ic)] / sdv
    }
  } else {
    Xbig <- kronecker(diag(p), Psi_cs)            # (N*p) x (p^2 * M)
    Ybig <- as.vector(rhs)                        # stack node responses
    group <- rep(seq_len(p * p), each = M)        # one group per (target, source)
    fit <- tryCatch(ADSIHT(Xbig, matrix(Ybig), group, ic.type = "dsic"),
                    error = function(e) NULL)
    if (is.null(fit)) {
      failed <- p
    } else {
      beta <- fit$beta[, which.min(fit$ic)]
      pM <- p * M
      for (j in seq_len(p)) theta[, j] <- beta[((j - 1L) * pM + 1):(j * pM)] / sdv
    }
  }

  adj <- matrix(0L, p, p)
  jac <- matrix(0, p, p)
  dphi <- basis$dfuncs
  for (j in seq_len(p)) for (i in seq_len(p)) {
    cols <- (i - 1L) * M + seq_len(M)
    th <- theta[cols, j]
    if (sqrt(sum(th^2)) >= edge_tol) adj[j, i] <- 1L
    jac[j, i] <- sum(vapply(seq_len(M),
                            function(m) th[m] * dphi[[m]](x_ref[i]), 0))
  }
  diag(adj) <- 0L
  diag(jac) <- 0
  list(adj = adj, jac = jac, theta = theta, x_ref = x_ref, method = method,
       failed_nodes = failed, basis = basis, M_ord = M)
}

## ============================ 3. metrics =====================================
edge_metrics <- function(est_adj, true_adj, est_jac, true_jac) {
  p <- nrow(true_adj)
  off <- which(row(true_adj) != col(true_adj))
  e <- est_adj[off]; t <- true_adj[off]
  TP <- sum(e == 1 & t == 1); FP <- sum(e == 1 & t == 0)
  TN <- sum(e == 0 & t == 0); FN <- sum(e == 0 & t == 1)
  pr <- ifelse(TP + FP == 0, 0, TP / (TP + FP))
  re <- ifelse(TP + FN == 0, 0, TP / (TP + FN))
  den <- sqrt(as.numeric(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  mcc <- ifelse(den == 0, 0, (TP * TN - FP * FN) / den)
  tp_idx <- off[e == 1 & t == 1]
  sign_acc <- if (length(tp_idx) == 0L) NA_real_ else
    mean(sign(est_jac[tp_idx]) == sign(true_jac[tp_idx]))
  list(TP = TP, FP = FP, FN = FN, Precision = pr, Recall = re, MCC = mcc,
       SignAcc = sign_acc)
}

## ============================ 4. network plotting ============================
# True vs inferred directed signed network, drawn on the ACTIVE graphics device
# (no file is written). Edge colour = activation (red) / inhibition (blue);
# inferred edges styled by status: TP solid, FP dashed, missed true edges (FN)
# dotted grey. Node size grows with true out-degree (hubs are larger).
plot_pss_networks <- function(sys, fit) {
  p <- sys$p
  vname <- paste0("x", seq_len(p))
  true_jac <- sys$A + 2 * sys$B * matrix(fit$x_ref, p, p, byrow = TRUE)

  build_graph <- function(adj, jac, ref_adj = NULL) {
    idx <- which(adj == 1L, arr.ind = TRUE)   # adj[j, i] = edge i -> j
    if (nrow(idx) == 0L) {
      df <- data.frame(from = character(0), to = character(0),
                       color = character(0), lty = numeric(0))
    } else {
      status <- if (is.null(ref_adj)) rep("TP", nrow(idx)) else
        ifelse(ref_adj[adj == 1L] == 1L, "TP", "FP")
      df <- data.frame(from = vname[idx[, "col"]], to = vname[idx[, "row"]],
                       color = ifelse(jac[adj == 1L] >= 0, "#C0392B", "#2E6F9E"),
                       lty = ifelse(status == "FP", 2, 1),
                       stringsAsFactors = FALSE)
    }
    graph_from_data_frame(df, vertices = vname, directed = TRUE)
  }

  g_true <- build_graph(sys$adj, true_jac)
  g_est <- build_graph(fit$adj, fit$jac, ref_adj = sys$adj)
  fn <- which(sys$adj == 1L & fit$adj == 0L, arr.ind = TRUE)  # missed edges
  if (nrow(fn) > 0L) {
    g_est <- add_edges(g_est, t(cbind(vname[fn[, "col"]], vname[fn[, "row"]])),
                       color = "grey70", lty = 3)
  }

  vsize <- 14 + 18 * (sys$outdeg_true / max(1, max(sys$outdeg_true)))
  lay <- layout_in_circle(g_true)
  draw <- function(g, title) {
    plot(g, layout = lay, vertex.size = vsize, vertex.color = "grey92",
         vertex.frame.color = "grey50", vertex.label.color = "grey10",
         vertex.label.cex = 0.9, edge.color = E(g)$color, edge.lty = E(g)$lty,
         edge.width = 1.8, edge.arrow.size = 0.45, edge.curved = 0.12)
    title(title, cex.main = 1.1)
  }

  op <- par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))
  on.exit(par(op))
  draw(g_true, "True network")
  draw(g_est, sprintf("PSS-Net inferred (%s)\nsolid TP, dashed FP, dotted FN",
                      fit$method))
  invisible(list(true = g_true, inferred = g_est))
}

## ============================ 5. run the demo ================================
run_pss_demo <- function(p = 10L, topology = "scalefree", N = 200L, snr = 30,
                         seed = 1L, method = "joint", basis = make_basis(),
                         plot = interactive()) {
  sys <- simulate_pss_nonlinear(p = p, topology = topology, N = N, snr = snr,
                                seed = seed)
  t0 <- system.time(fit <- pss_net(sys$U, sys$X_obs, basis = basis,
                                   method = method, x_ref = sys$x_wt))[["elapsed"]]
  true_jac <- sys$A + 2 * sys$B * matrix(sys$x_wt, sys$p, sys$p, byrow = TRUE)
  m <- edge_metrics(fit$adj, sys$adj, fit$jac, true_jac)

  cat("================ PSS-Net demo ================\n")
  cat(sprintf("system        : p = %d, %s topology, true edges = %d (density %.1f%%)\n",
              sys$p, sys$topology, sum(sys$adj),
              100 * sum(sys$adj) / (sys$p * (sys$p - 1))))
  cat(sprintf("PSS data      : N = %d conditions, SNR = %g, sigma = %.3f\n",
              sys$n_eff, sys$snr, sys$sigma))
  cat(sprintf("inference     : %s ADSIHT, %s basis [%s], %.2fs, %d node(s) empty\n",
              fit$method, fit$basis$type, paste(fit$basis$labels, collapse = ", "),
              t0, fit$failed_nodes))
  cat("---------------- edge recovery ----------------\n")
  cat(sprintf("TP=%d  FP=%d  FN=%d\n", m$TP, m$FP, m$FN))
  cat(sprintf("Precision=%.3f  Recall=%.3f  MCC=%.3f  SignAcc=%.3f\n",
              m$Precision, m$Recall, m$MCC, m$SignAcc))
  cat("---------------- hubs (out-degree) ------------\n")
  od_true <- sys$outdeg_true; od_est <- colSums(fit$adj)
  for (i in order(od_true, decreasing = TRUE)[seq_len(min(3L, p))])
    cat(sprintf("  x%-2d  true out-degree %d,  inferred %d\n",
                i, od_true[i], od_est[i]))
  cat("-----------------------------------------------\n")
  if (isTRUE(plot)) {
    plot_pss_networks(sys, fit)
  } else {
    cat("(call plot_pss_networks(sys, fit) on an active device to draw the network)\n")
  }
  invisible(list(sys = sys, fit = fit, metrics = m))
}

if (sys.nframe() == 0L) {
  invisible(run_pss_demo())
}

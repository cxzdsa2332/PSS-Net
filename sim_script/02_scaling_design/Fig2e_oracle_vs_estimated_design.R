rm(list = ls())

################################################################################
# Fig2e_oracle_vs_estimated_design.R -- package-backed D-optimal augmentation
#
# Question: does D-optimal design retain its advantage when the steady-state
# response map is estimated from a finite pilot experiment rather than supplied
# by the true system? Every strategy receives the same random pilot conditions,
# and pilot observations count toward the reported total budget.
#
# Strategies after the shared pilot:
#   random       random continuation from the common candidate pool
#   maximin      space-filling continuation conditional on the pilot inputs
#   oracle_dopt  D-optimal continuation using the true local Jacobian
#   pilot_dopt   D-optimal continuation using a ridge estimate of dx*/du from
#                the noisy pilot PSS observations
#
# Output: results/sim_results/Fig2e_oracle_vs_estimated_design.csv
#
# Terminology/literature boundary:
#   This is a two-stage, pilot-informed locally D-optimal design, not a verbatim
#   implementation of a named "pilot D-optimal" paper. Exact D-optimal
#   augmentation is delegated to AlgDesign::optFederov() with protected pilot
#   runs; only the PSS response map and candidate feature construction are
#   project-specific. Unlike a full adaptive design, H is estimated once from
#   the pilot and is not refitted after each newly observed batch. See
#   ref/pilot_doptimal_literature.md.
################################################################################

suppressMessages({
  library(ADSIHT)
  library(AlgDesign)
  library(deSolve)
})

set.seed(2206)
M_ord <- 2L

steady_one <- function(sys, u, x0 = NULL, t_max = 2000) {
  p <- sys$p
  if (is.null(x0)) x0 <- sys$x_wt
  deriv <- function(t, x, parms) {
    inter <- numeric(p)
    for (j in seq_len(p)) {
      others <- setdiff(seq_len(p), j)
      inter[j] <- sum(sys$A[j, others] * x[others] +
                        sys$Bq[j, others] * x[others]^2)
    }
    list(sys$mu - sys$gamma * x + inter + u)
  }
  out <- tryCatch(
    ode(y = x0, times = c(0, t_max), func = deriv, parms = NULL,
        method = "lsoda", rtol = 1e-9, atol = 1e-11),
    error = function(e) NULL)
  if (is.null(out)) return(rep(NA_real_, p))
  as.numeric(out[2, -1])
}

steady_many <- function(sys, U) {
  t(apply(U, 1, function(u) steady_one(sys, u)))
}

make_system <- function(p = 8L, n_in = 2L, seed = 1L) {
  set.seed(seed)
  A <- Bq <- matrix(0, p, p)
  for (j in seq_len(p)) {
    src <- sample(setdiff(seq_len(p), j), n_in)
    A[j, src] <- runif(n_in, 0.3, 0.7) * sample(c(-1, 1), n_in, TRUE)
    curved <- src[seq_len(max(1L, floor(n_in / 2L)))]
    Bq[j, curved] <- runif(length(curved), 0.10, 0.30) *
      sample(c(-1, 1), length(curved), TRUE)
  }
  sys <- list(
    p = p, A = A, Bq = Bq, gamma = runif(p, 3.0, 4.0),
    mu = runif(p, 1.0, 2.0), adj = ((A != 0) | (Bq != 0)) * 1L
  )
  sys$x_wt <- steady_one(sys, rep(0, p), x0 = sys$mu / sys$gamma)
  sys
}

jacobian_at_wt <- function(sys) {
  J <- matrix(0, sys$p, sys$p)
  for (j in seq_len(sys$p)) for (i in seq_len(sys$p)) {
    J[j, i] <- if (i == j) -sys$gamma[j] else
      sys$A[j, i] + 2 * sys$Bq[j, i] * sys$x_wt[i]
  }
  J
}

psi_row <- function(x) as.vector(sapply(x, function(z) z^seq_len(M_ord)))
# Pass only x/x^2 columns in data; the explicit AlgDesign formula `~ .` adds one
# intercept in the standard model matrix. Do not duplicate a constant column.
aug_row <- function(x) psi_row(x)

make_pool <- function(p, n) matrix(runif(n * p, -0.4, 0.8), n, p)

# Exact D-optimal augmentation via the CRAN package AlgDesign. Pilot rows are
# protected by `augment = TRUE`; the Fedorov exchange algorithm chooses the
# remaining runs from the finite candidate feature matrix. Each total budget is
# optimized separately rather than treated as a prefix of a custom Wynn order.
select_dopt_augment <- function(Phi_pool, Phi_pilot, N_total, seed) {
  n_pilot <- nrow(Phi_pilot)
  Fx <- rbind(Phi_pilot, Phi_pool)
  Fx <- sweep(Fx, 2, colMeans(Fx))
  fx_sd <- apply(Fx, 2, sd)
  fx_sd[!is.finite(fx_sd) | fx_sd < 1e-10] <- 1
  Fx <- sweep(Fx, 2, fx_sd, "/")
  # Exact equilibrium features can contain deterministic linear dependencies.
  # D-optimality is defined on the estimable column space, so pass an explicit
  # full-rank basis to AlgDesign rather than stabilizing a singular determinant
  # with a project-specific ridge term.
  qFx <- qr(Fx, tol = 1e-9)
  if (qFx$rank < ncol(Fx)) {
    keep <- sort(qFx$pivot[seq_len(qFx$rank)])
    Fx <- Fx[, keep, drop = FALSE]
  }
  set.seed(seed)
  fit <- tryCatch(
    AlgDesign::optFederov(
      frml = ~ ., data = as.data.frame(Fx), nTrials = N_total, criterion = "D",
      augment = TRUE, rows = seq_len(n_pilot),
      maxIteration = 100, nRepeats = 1
    ),
    error = function(e) stop(sprintf(
      "AlgDesign failed (pilot=%d, total=%d, rank=%d/%d, kappa=%.3g): %s",
      n_pilot, N_total, qr(Fx)$rank, ncol(Fx), kappa(Fx), conditionMessage(e)
    ))
  )
  idx <- fit$rows[fit$rows > n_pilot] - n_pilot
  if (length(idx) != N_total - n_pilot) {
    stop(sprintf("AlgDesign returned %d augmented runs; expected %d (pilot=%d, total=%d, rows=%s).",
                 length(idx), N_total - n_pilot, n_pilot, N_total,
                 paste(fit$rows, collapse = ",")))
  }
  as.integer(idx)
}

continue_maximin <- function(U_pool, U_pilot, n_add) {
  if (n_add <= 0L) return(integer(0))
  min_dist <- rep(Inf, nrow(U_pool))
  for (k in seq_len(nrow(U_pilot))) {
    d <- sqrt(rowSums((sweep(U_pool, 2, U_pilot[k, ]))^2))
    min_dist <- pmin(min_dist, d)
  }
  selected <- integer(0)
  for (k in seq_len(n_add)) {
    nxt <- which.max(min_dist)
    selected <- c(selected, nxt)
    d <- sqrt(rowSums((sweep(U_pool, 2, U_pool[nxt, ]))^2))
    min_dist <- pmin(min_dist, d)
    min_dist[selected] <- -Inf
  }
  selected
}

# Estimate the local steady-state response x*(u) = x_ref + (u-u_ref) H from
# noisy pilot observations. This is deliberately simpler than refitting the
# full nonlinear ODE and corresponds to a feasible first adaptive batch.
predict_pool_from_pilot <- function(U_pilot, X_pilot, U_pool, lambda = 0.25) {
  Uc <- sweep(U_pilot, 2, colMeans(U_pilot))
  Xc <- sweep(X_pilot, 2, colMeans(X_pilot))
  H <- solve(crossprod(Uc) + diag(lambda, ncol(Uc)), crossprod(Uc, Xc))
  sweep(U_pool, 2, colMeans(U_pilot)) %*% H +
    matrix(colMeans(X_pilot), nrow(U_pool), ncol(X_pilot), byrow = TRUE)
}

infer_network <- function(sys, U, X) {
  p <- sys$p
  Psi <- t(apply(X, 1, psi_row))
  Psi_c <- sweep(Psi, 2, colMeans(Psi))
  sdv <- apply(Psi_c, 2, sd)
  sdv[!is.finite(sdv) | sdv < 1e-10] <- 1e-10
  Psi_cs <- sweep(Psi_c, 2, sdv, "/")
  group <- rep(seq_len(p), each = M_ord)
  adj <- matrix(0L, p, p)
  failed <- 0L
  for (j in seq_len(p)) {
    rhs <- -(U[, j] - mean(U[, j]))
    fit <- tryCatch(ADSIHT(Psi_cs, matrix(rhs), group, ic.type = "dsic"),
                    error = function(e) NULL)
    if (is.null(fit) || length(fit$ic) == 0L) {
      failed <- failed + 1L
      next
    }
    th <- fit$beta[, which.min(fit$ic)] / sdv
    gn <- sapply(seq_len(p), function(i) {
      cols <- (i - 1L) * M_ord + seq_len(M_ord)
      sqrt(sum(th[cols]^2))
    })
    adj[j, gn >= 1e-8] <- 1L
  }
  diag(adj) <- 0L
  list(adj = adj, failed_nodes = failed)
}

edge_metrics <- function(est, truth) {
  off <- which(row(truth) != col(truth))
  e <- est[off]
  t <- truth[off]
  TP <- sum(e == 1 & t == 1); FP <- sum(e == 1 & t == 0)
  TN <- sum(e == 0 & t == 0); FN <- sum(e == 0 & t == 1)
  pr <- ifelse(TP + FP == 0, 0, TP / (TP + FP))
  re <- ifelse(TP + FN == 0, 0, TP / (TP + FN))
  den <- sqrt(as.numeric(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  c(Precision = pr, Recall = re,
    MCC = ifelse(den == 0, 0, (TP * TN - FP * FN) / den))
}

pilot_grid <- c(8L, 12L, 16L)
total_grid <- c(20L, 30L, 40L, 60L)
strategies <- c("random", "maximin", "oracle_dopt", "pilot_dopt")
R <- 20L
sigma <- 0.04
n_pool <- 2500L

rows <- list()
for (seed in seq_len(R)) {
  sys <- make_system(seed = 2600L + seed)
  if (any(!is.finite(sys$x_wt)) || any(sys$x_wt <= 0)) next

  set.seed(3600L + seed)
  U_pilot_max <- make_pool(sys$p, max(pilot_grid))
  U_pilot_max[1, ] <- 0
  X_pilot_true_max <- steady_many(sys, U_pilot_max)
  E_pilot <- matrix(rnorm(length(X_pilot_true_max), sd = sigma),
                    nrow(X_pilot_true_max), sys$p)
  X_pilot_obs_max <- X_pilot_true_max + E_pilot

  U_pool <- make_pool(sys$p, n_pool)
  E_pool <- matrix(rnorm(n_pool * sys$p, sd = sigma), n_pool, sys$p)
  J_true <- jacobian_at_wt(sys)
  X_pool_oracle <- t(sys$x_wt - solve(J_true, t(U_pool)))
  Phi_pool_oracle <- t(apply(X_pool_oracle, 1, aug_row))

  for (pilot_n in pilot_grid) {
    U_pilot <- U_pilot_max[seq_len(pilot_n), , drop = FALSE]
    X_pilot_true <- X_pilot_true_max[seq_len(pilot_n), , drop = FALSE]
    X_pilot_obs <- X_pilot_obs_max[seq_len(pilot_n), , drop = FALSE]
    Phi_pilot_oracle <- t(apply(X_pilot_true, 1, aug_row))
    X_pool_pilot <- predict_pool_from_pilot(U_pilot, X_pilot_obs, U_pool)
    Phi_pool_pilot <- t(apply(X_pool_pilot, 1, aug_row))
    Phi_pilot_est <- t(apply(X_pilot_obs, 1, aug_row))
    budgets <- total_grid[total_grid > pilot_n]
    n_add_max <- max(total_grid) - pilot_n

    random_order <- seq_len(n_add_max)
    maximin_order <- continue_maximin(U_pool, U_pilot, n_add_max)
    index_by_strategy <- list(
      random = setNames(lapply(budgets, function(N) {
        random_order[seq_len(N - pilot_n)]
      }), budgets),
      maximin = setNames(lapply(budgets, function(N) {
        maximin_order[seq_len(N - pilot_n)]
      }), budgets),
      oracle_dopt = setNames(lapply(budgets, function(N) {
        select_dopt_augment(Phi_pool_oracle, Phi_pilot_oracle, N,
                            seed = 460000L + 1000L * seed + 10L * pilot_n + N)
      }), budgets),
      pilot_dopt = setNames(lapply(budgets, function(N) {
        select_dopt_augment(Phi_pool_pilot, Phi_pilot_est, N,
                            seed = 560000L + 1000L * seed + 10L * pilot_n + N)
      }), budgets)
    )

    # Simulate each selected candidate at most once within this pilot setting,
    # even though Fedorov augmentation is optimized separately at each budget.
    union_idx <- unique(unlist(index_by_strategy, use.names = FALSE))
    X_union_true <- steady_many(sys, U_pool[union_idx, , drop = FALSE])
    X_union_obs <- X_union_true + E_pool[union_idx, , drop = FALSE]

    for (strategy in strategies) {
      for (N_total in budgets) {
        idx <- index_by_strategy[[strategy]][[as.character(N_total)]]
        n_add <- N_total - pilot_n
        U_all <- rbind(U_pilot, U_pool[idx, , drop = FALSE])
        X_all <- rbind(X_pilot_obs,
                       X_union_obs[match(idx, union_idx), , drop = FALSE])
        ok <- apply(X_all, 1, function(z) all(is.finite(z)) && all(z > 0))
        if (sum(ok) < 8L) {
          mets <- c(Precision = NA, Recall = NA, MCC = NA)
          failed_nodes <- sys$p
        } else {
          fit <- infer_network(sys, U_all[ok, , drop = FALSE],
                               X_all[ok, , drop = FALSE])
          mets <- edge_metrics(fit$adj, sys$adj)
          failed_nodes <- fit$failed_nodes
        }
        rows[[length(rows) + 1L]] <- data.frame(
          seed = seed, pilot_n = pilot_n, N_total = N_total,
          N_adaptive = n_add, strategy = strategy, n_eff = sum(ok),
          Precision = mets["Precision"], Recall = mets["Recall"],
          MCC = mets["MCC"], failed_nodes = failed_nodes,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  cat(sprintf("Fig2e seed %d/%d done\n", seed, R))
}

df <- do.call(rbind, rows)
rownames(df) <- NULL
out <- "results/sim_results/Fig2e_oracle_vs_estimated_design.csv"
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
write.csv(df, out, row.names = FALSE)

cat("\nSaved:", out, "\n")
print(aggregate(cbind(MCC, Precision, Recall) ~ pilot_n + N_total + strategy,
                df, mean, na.rm = TRUE), row.names = FALSE)

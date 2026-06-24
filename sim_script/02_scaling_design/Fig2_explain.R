rm(list = ls())

################################################################################
# Fig2_explain.R -- conceptual 3D perturbation-design schematic
#
# Purpose: draw a visual explanation of four perturbation-design strategies in
#          a 3D input space: random, maximin, oracle D-optimal, and D-optimal
#          based on a noisy pilot-estimated response map. All four strategies
#          share the same pilot points. This script is conceptual and does not
#          write numeric simulation results.
#
# Output:  Fig2_design_concept_random, Fig2_design_concept_maximin,
#          Fig2_design_concept_oracle, Fig2_design_concept_pilot,
#          Fig2_design_concept_dopt (oracle compatibility alias), and
#          Fig2_design_concept in the workspace.
################################################################################

suppressMessages({
  library(AlgDesign)
  library(ggplot2)
  library(grid)
})

set.seed(2026)

## ------------------------------------------------------------- Design rules ----
n_pool <- 3000L
n_design <- 28L
n_pilot <- 8L
u_lo <- -1
u_hi <- 1

pool <- matrix(runif(n_pool * 3L, u_lo, u_hi), n_pool, 3L)
colnames(pool) <- c("u1", "u2", "u3")

# A toy nonlinear steady-state response map. Oracle D-optimal knows this map;
# pilot D-optimal estimates only its local linear response from noisy pilot data.
true_response <- function(U) {
  cbind(
    0.55 + 0.70 * U[, 1] - 0.20 * U[, 3] + 0.28 * U[, 2]^2,
    0.70 - 0.35 * U[, 1] + 0.62 * U[, 2] + 0.22 * U[, 1] * U[, 3],
    0.60 + 0.18 * U[, 1] - 0.32 * U[, 2] + 0.72 * U[, 3] + 0.25 * U[, 1]^2
  )
}

state_features <- function(X) {
  # Pass x/x^2 as data columns; the explicit AlgDesign formula `~ .` adds the
  # standard intercept, so no constant column is constructed here.
  cbind(X, X^2)
}

pilot_idx <- c(which.min(rowSums(pool^2)),
               sample(setdiff(seq_len(nrow(pool)), which.min(rowSums(pool^2))),
                      n_pilot - 1L))
U_pilot <- pool[pilot_idx, , drop = FALSE]
candidate_pool <- pool[-pilot_idx, , drop = FALSE]
n_add <- n_design - n_pilot

design_random <- function(U_pilot, pool, n_add) {
  idx <- sample(seq_len(nrow(pool)), n_add)
  list(U = rbind(U_pilot, pool[idx, , drop = FALSE]), idx = idx)
}

design_maximin <- function(U_pilot, pool, n_add) {
  # Greedy maximin continuation conditional on the shared pilot inputs.
  mind <- rep(Inf, nrow(pool))
  for (i in seq_len(nrow(U_pilot))) {
    mind <- pmin(mind, sqrt(rowSums((sweep(pool, 2, U_pilot[i, ]))^2)))
  }
  idx <- integer(0)
  for (k in seq_len(n_add)) {
    nxt <- which.max(mind)
    idx <- c(idx, nxt)
    d_new <- sqrt(rowSums((sweep(pool, 2, pool[nxt, ]))^2))
    mind <- pmin(mind, d_new)
    mind[idx] <- -Inf
  }
  list(U = rbind(U_pilot, pool[idx, , drop = FALSE]), idx = idx)
}

design_dopt <- function(U_pilot, pool, Phi_pilot, Phi_pool, n_add,
                        seed = 1L) {
  # Package-backed exact D-optimal augmentation. The shared pilot rows are
  # protected from exchange; AlgDesign chooses the remaining candidate rows.
  n_pilot_local <- nrow(Phi_pilot)
  Fx <- rbind(Phi_pilot, Phi_pool)
  Fx <- sweep(Fx, 2, colMeans(Fx))
  fx_sd <- apply(Fx, 2, sd)
  fx_sd[!is.finite(fx_sd) | fx_sd < 1e-10] <- 1
  Fx <- sweep(Fx, 2, fx_sd, "/")
  qFx <- qr(Fx, tol = 1e-9)
  if (qFx$rank < ncol(Fx)) {
    keep <- sort(qFx$pivot[seq_len(qFx$rank)])
    Fx <- Fx[, keep, drop = FALSE]
  }
  set.seed(seed)
  fit <- AlgDesign::optFederov(
    frml = ~ ., data = as.data.frame(Fx),
    nTrials = n_pilot_local + n_add,
    criterion = "D", augment = TRUE,
    rows = seq_len(n_pilot_local),
    maxIteration = 100, nRepeats = 1
  )
  idx <- fit$rows[fit$rows > n_pilot_local] - n_pilot_local
  list(U = rbind(U_pilot, pool[idx, , drop = FALSE]), idx = idx)
}

X_pilot_true <- true_response(U_pilot)
X_pool_true <- true_response(candidate_pool)
Phi_pilot_oracle <- state_features(X_pilot_true)
Phi_pool_oracle <- state_features(X_pool_true)

# Pilot-estimated local response map, matching the logic of Fig2e.
X_pilot_obs <- X_pilot_true + matrix(rnorm(length(X_pilot_true), sd = 0.10),
                                     nrow(X_pilot_true), ncol(X_pilot_true))
Uc <- sweep(U_pilot, 2, colMeans(U_pilot))
Xc <- sweep(X_pilot_obs, 2, colMeans(X_pilot_obs))
H_hat <- solve(crossprod(Uc) + diag(0.25, ncol(Uc)), crossprod(Uc, Xc))
X_pool_hat <- sweep(candidate_pool, 2, colMeans(U_pilot)) %*% H_hat +
  matrix(colMeans(X_pilot_obs), nrow(candidate_pool), ncol(X_pilot_obs), byrow = TRUE)
Phi_pilot_hat <- state_features(X_pilot_obs)
Phi_pool_hat <- state_features(X_pool_hat)

designs <- list(
  "Random" = design_random(U_pilot, candidate_pool, n_add),
  "Maximin" = design_maximin(U_pilot, candidate_pool, n_add),
  "Oracle D-opt" = design_dopt(U_pilot, candidate_pool,
                                Phi_pilot_oracle, Phi_pool_oracle, n_add,
                                seed = 20261L),
  "Pilot D-opt" = design_dopt(U_pilot, candidate_pool,
                               Phi_pilot_hat, Phi_pool_hat, n_add,
                               seed = 20262L)
)

## -------------------------------------------------------------- 3D drawing ----
project_points <- function(x, y, z, theta = 38, phi = 24) {
  # Lightweight orthographic projection. This avoids base::persp(), which opens
  # a graphics device even when plot = FALSE on some systems.
  th <- theta * pi / 180
  ph <- phi * pi / 180
  xr <- x * cos(th) - y * sin(th)
  yr <- x * sin(th) * sin(ph) + y * cos(th) * sin(ph) + z * cos(ph)
  list(x = xr, y = yr)
}

cube_edges <- rbind(
  c(-1, -1, -1,  1, -1, -1), c(-1,  1, -1,  1,  1, -1),
  c(-1, -1,  1,  1, -1,  1), c(-1,  1,  1,  1,  1,  1),
  c(-1, -1, -1, -1,  1, -1), c( 1, -1, -1,  1,  1, -1),
  c(-1, -1,  1, -1,  1,  1), c( 1, -1,  1,  1,  1,  1),
  c(-1, -1, -1, -1, -1,  1), c( 1, -1, -1,  1, -1,  1),
  c(-1,  1, -1, -1,  1,  1), c( 1,  1, -1,  1,  1,  1)
)
colnames(cube_edges) <- c("x", "y", "z", "xend", "yend", "zend")

cube_df <- do.call(rbind, lapply(seq_len(nrow(cube_edges)), function(i) {
  a <- project_points(cube_edges[i, "x"], cube_edges[i, "y"], cube_edges[i, "z"])
  b <- project_points(cube_edges[i, "xend"], cube_edges[i, "yend"], cube_edges[i, "zend"])
  data.frame(x = a$x, y = a$y, xend = b$x, yend = b$y)
}))

axis_raw <- data.frame(
  label = c("u[1]", "u[2]", "u[3]"),
  x = c(-1, -1, -1), y = c(-1, -1, -1), z = c(-1, -1, -1),
  xend = c(1.18, -1, -1), yend = c(-1, 1.18, -1), zend = c(-1, -1, 1.18)
)
axis_df <- do.call(rbind, lapply(seq_len(nrow(axis_raw)), function(i) {
  a <- project_points(axis_raw$x[i], axis_raw$y[i], axis_raw$z[i])
  b <- project_points(axis_raw$xend[i], axis_raw$yend[i], axis_raw$zend[i])
  data.frame(x = a$x, y = a$y, xend = b$x, yend = b$y,
             label = axis_raw$label[i])
}))

make_panel_df <- function(design, strategy) {
  mat <- design$U
  p <- project_points(mat[, 1], mat[, 2], mat[, 3])
  data.frame(strategy = strategy, px = p$x, py = p$y,
             phase = rep(c("Shared pilot", "Selected next"),
                         c(n_pilot, nrow(mat) - n_pilot)),
             depth = mat[, 3], u1 = mat[, 1], u2 = mat[, 2], u3 = mat[, 3])
}

panel_df <- do.call(rbind, Map(make_panel_df, designs, names(designs)))
panel_df$strategy <- factor(panel_df$strategy, levels = names(designs))

draw_design_panel <- function(strategy_name) {
  d <- panel_df[panel_df$strategy == strategy_name, ]
  ggplot() +
    geom_segment(data = cube_df,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 linewidth = 0.35, color = "grey78") +
    geom_segment(data = axis_df,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 arrow = arrow(length = unit(2.1, "mm"), type = "closed"),
                 linewidth = 0.42, color = "grey35") +
    geom_text(data = axis_df,
              aes(x = xend, y = yend, label = label),
              parse = TRUE, hjust = -0.15, vjust = -0.1,
              size = 3.0, color = "grey20") +
    geom_point(data = d, aes(x = px, y = py, color = depth),
               size = 2.35, alpha = 0.95) +
    geom_point(data = d[d$phase == "Shared pilot", ],
               aes(x = px, y = py), shape = 1, size = 3.25,
               stroke = 0.8, color = "grey10") +
    scale_color_gradient2(low = "#2E6F9E", mid = "#F5F6F4", high = "#C0392B",
                          midpoint = 0, guide = "none") +
    labs(title = strategy_name) +
    coord_equal(clip = "off") +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      # Doubled left/right margins widen the gaps between the four panels.
      plot.margin = margin(7, 20, 7, 20)
    )
}

Fig2_design_concept_random <- draw_design_panel("Random")
Fig2_design_concept_maximin <- draw_design_panel("Maximin")
Fig2_design_concept_oracle <- draw_design_panel("Oracle D-opt")
Fig2_design_concept_pilot <- draw_design_panel("Pilot D-opt")
Fig2_design_concept_dopt <- Fig2_design_concept_oracle

if (requireNamespace("patchwork", quietly = TRUE)) {
  Fig2_design_concept <- Fig2_design_concept_random +
    Fig2_design_concept_maximin +
    Fig2_design_concept_oracle +
    Fig2_design_concept_pilot +
    patchwork::plot_layout(nrow = 1)
} else {
  Fig2_design_concept <- Fig2_design_concept_pilot
}

Fig2_design_concept

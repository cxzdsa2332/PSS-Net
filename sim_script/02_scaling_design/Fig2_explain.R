rm(list = ls())

################################################################################
# Fig2_explain.R -- conceptual 3D perturbation-design schematic
#
# Purpose: draw a visual explanation of three perturbation-design strategies in
#          a 3D input space: random sampling, maximin space filling, and
#          sequential D-optimal design. This script is conceptual and does not
#          write numeric simulation results.
#
# Output:  Fig2_design_concept_random, Fig2_design_concept_maximin,
#          Fig2_design_concept_dopt, and Fig2_design_concept in the workspace.
################################################################################

suppressMessages({
  library(ggplot2)
  library(grid)
})

set.seed(2026)

## ------------------------------------------------------------- Design rules ----
n_pool <- 3000L
n_design <- 28L
u_lo <- -1
u_hi <- 1

pool <- matrix(runif(n_pool * 3L, u_lo, u_hi), n_pool, 3L)
colnames(pool) <- c("u1", "u2", "u3")

feature_row <- function(u) {
  # A small nonlinear feature dictionary: intercept, main effects, pairwise
  # products and quadratic terms. D-optimal design spreads points where these
  # columns become informative, often near boundaries/corners.
  c(1, u, u[1] * u[2], u[1] * u[3], u[2] * u[3], u^2)
}

design_random <- function(pool, n) {
  pool[sample(seq_len(nrow(pool)), n), , drop = FALSE]
}

design_maximin <- function(pool, n) {
  # Greedy maximin in raw perturbation space.
  start_idx <- which.min(rowSums((sweep(pool, 2, c(u_lo, u_lo, u_lo)))^2))
  idx <- start_idx
  mind <- sqrt(rowSums((sweep(pool, 2, pool[start_idx, ]))^2))
  mind[idx] <- -Inf
  for (k in 2:n) {
    nxt <- which.max(mind)
    idx <- c(idx, nxt)
    d_new <- sqrt(rowSums((sweep(pool, 2, pool[nxt, ]))^2))
    mind <- pmin(mind, d_new)
    mind[idx] <- -Inf
  }
  pool[idx, , drop = FALSE]
}

design_dopt <- function(pool, n, lambda = 1e-2) {
  # Greedy D-optimal design in feature space, using Sherman-Morrison updates.
  Phi <- t(apply(pool, 1, feature_row))
  q <- ncol(Phi)
  idx <- which.min(rowSums(pool^2))  # include a near-baseline point first
  Minv <- diag(1 / lambda, q)
  add_point <- function(i) {
    phi <- Phi[i, ]
    Mv <- Minv %*% phi
    denom <- as.numeric(1 + phi %*% Mv)
    Minv <<- Minv - (Mv %*% t(Mv)) / denom
  }
  add_point(idx)
  avail <- rep(TRUE, nrow(pool))
  avail[idx] <- FALSE
  for (k in 2:n) {
    PM <- Phi %*% Minv
    score <- rowSums(PM * Phi)
    score[!avail] <- -Inf
    nxt <- which.max(score)
    idx <- c(idx, nxt)
    add_point(nxt)
    avail[nxt] <- FALSE
  }
  pool[idx, , drop = FALSE]
}

designs <- list(
  "Random" = design_random(pool, n_design),
  "Maximin" = design_maximin(pool, n_design),
  "D-optimal" = design_dopt(pool, n_design)
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

make_panel_df <- function(mat, strategy) {
  p <- project_points(mat[, 1], mat[, 2], mat[, 3])
  data.frame(strategy = strategy, px = p$x, py = p$y,
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
    scale_color_gradient2(low = "#2E6F9E", mid = "#F5F6F4", high = "#C0392B",
                          midpoint = 0, guide = "none") +
    labs(title = strategy_name) +
    coord_equal(clip = "off") +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      plot.margin = margin(7, 10, 7, 10)
    )
}

Fig2_design_concept_random <- draw_design_panel("Random")
Fig2_design_concept_maximin <- draw_design_panel("Maximin")
Fig2_design_concept_dopt <- draw_design_panel("D-optimal")

if (requireNamespace("patchwork", quietly = TRUE)) {
  Fig2_design_concept <- Fig2_design_concept_random +
    Fig2_design_concept_maximin +
    Fig2_design_concept_dopt +
    patchwork::plot_layout(nrow = 1) +
    patchwork::plot_annotation(
      title = "Perturbation design strategies in a 3D input space",
      subtitle = "Random samples cluster by chance; maximin fills u-space; D-optimal selects feature-informative boundary points"
    ) &
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 8.5, color = "grey30")
    )
} else {
  Fig2_design_concept <- Fig2_design_concept_dopt
}

Fig2_design_concept

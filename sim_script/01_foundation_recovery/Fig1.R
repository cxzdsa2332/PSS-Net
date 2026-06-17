rm(list = ls())

################################################################################
# Fig1.R -- Figure 1 panel objects for PSS-Net
#
# Purpose: draw Fig1a as a visual ODE-to-PSS measurement panel. The panel uses
#          the historical 8-node parameter setting from sim_script/manual/standard0.R
#          and contrasts an additive linear ODE with a multiplicative gLV ODE.
#
# Input:   none
# Output:  Fig1a, Fig1a_legend, Fig1b_function_shape, Fig1b_identifiability,
#          and Fig1b objects in the workspace for later patchwork.
################################################################################

suppressMessages({
  library(deSolve)
  library(ggplot2)
})

set.seed(42)

## ---------------------------------------------------------------- Parameters ----
n_sp <- 8
sp_names <- paste0("Sp", seq_len(n_sp))

# Historical 8-node setting from sim_script/manual/standard0.R.
r_true <- c(0.8, 1.2, 0.6, 1.0, 0.7, 1.1, 0.5, 0.9)
gam_true <- c(1.5, 1.8, 1.2, 1.6, 1.4, 1.5, 1.0, 1.6)

A_true <- matrix(0, n_sp, n_sp)
A_true[1, 3] <-  0.40; A_true[1, 5] <- -0.30
A_true[2, 1] <-  0.30; A_true[2, 4] <- -0.40; A_true[2, 7] <-  0.20
A_true[3, 2] <- -0.30; A_true[3, 6] <-  0.30
A_true[4, 1] <-  0.40; A_true[4, 3] <- -0.20; A_true[4, 8] <-  0.30
A_true[5, 2] <-  0.30; A_true[5, 6] <- -0.20
A_true[6, 4] <-  0.20; A_true[6, 5] <-  0.30
A_true[7, 3] <- -0.30; A_true[7, 8] <-  0.40
A_true[8, 1] <-  0.20; A_true[8, 6] <- -0.20

ss_linear <- function(u) {
  as.numeric(solve(diag(gam_true) - A_true, r_true + u))
}

## -------------------------------------------------------------- ODE systems ----
deriv_additive <- function(t, state, parms) {
  x <- pmax(state, 0)
  list(r_true + as.numeric(A_true %*% x) - gam_true * x + parms$u)
}

deriv_glv <- function(t, state, parms) {
  x <- pmax(state, 0)
  rate <- r_true + as.numeric(A_true %*% x) - gam_true * x + parms$u
  list(x * rate)
}

simulate_trajectory <- function(u, condition, model, deriv_fun,
                                x0, times = seq(0, 20, length.out = 220)) {
  out <- ode(y = x0, times = times, func = deriv_fun, parms = list(u = u),
             method = "lsoda", rtol = 1e-9, atol = 1e-11)
  mat <- as.data.frame(out)
  names(mat)[-1] <- sp_names
  long <- reshape(mat, varying = sp_names, v.names = "abundance",
                  timevar = "node", times = sp_names,
                  idvar = "time", direction = "long")
  rownames(long) <- NULL
  long$model <- model
  long$condition <- condition
  long$u_norm <- sqrt(sum(u^2))
  long
}

## ---------------------------------------------------------- Perturbations ----
# Three mild perturbation conditions within the historical range used by the
# manual scripts; all share the same steady-state equation.
perturbations <- list(
  "baseline u = 0" = rep(0, n_sp),
  "single-node input" = c(0.45, 0, 0, 0, 0, 0, 0, 0),
  "mixed input" = c(0.30, -0.20, 0, 0.20, 0, -0.15, 0.25, 0)
)

x0_baseline <- rep(0.35, n_sp)
rows <- list()
idx <- 1L
for (model in c("Additive ODE", "gLV ODE")) {
  deriv_fun <- if (model == "Additive ODE") deriv_additive else deriv_glv
  baseline_traj <- simulate_trajectory(perturbations[["baseline u = 0"]],
                                       "baseline u = 0", model, deriv_fun,
                                       x0_baseline)
  rows[[idx]] <- baseline_traj
  idx <- idx + 1L

  # Perturbation-response trajectories start from the model-specific baseline
  # steady state, making the vertical sampling line read as post-perturbation PSS.
  baseline_ss <- baseline_traj[baseline_traj$time == max(baseline_traj$time), ]
  baseline_ss <- baseline_ss[match(sp_names, baseline_ss$node), "abundance"]

  for (condition in setdiff(names(perturbations), "baseline u = 0")) {
    rows[[idx]] <- simulate_trajectory(perturbations[[condition]], condition,
                                       model, deriv_fun, baseline_ss)
    idx <- idx + 1L
  }
}
traj_df <- do.call(rbind, rows)
traj_df$model <- factor(traj_df$model, levels = c("Additive ODE", "gLV ODE"))
traj_df$condition <- factor(traj_df$condition, levels = names(perturbations))
traj_df$node <- factor(traj_df$node, levels = sp_names)

# Measured PSS values: the terminal observations supplied to PSS-Net.
pss_df <- traj_df[ave(traj_df$time, traj_df$model, traj_df$condition,
                      traj_df$node, FUN = max) == traj_df$time, ]
pss_df$label_x <- max(traj_df$time) + 0.45

# Plotting rule for this visual panel: model rows receive distinct pale
# backgrounds; solid curves show ODE trajectories; a black vertical dashed line
# marks the post-perturbation sampling time, and terminal dots mark x*.
bg_df <- data.frame(
  model = factor(c("Additive ODE", "gLV ODE"), levels = levels(traj_df$model)),
  xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
  bg_fill = c("#F5F8E8", "#EEF6FA")
)

## ------------------------------------------------------------------- Fig1a ----
Fig1a_legend <- ggplot() +
  annotate("text", x = 0, y = 1,
           label = "ODE trajectories start at baseline steady state for perturbations; the black dashed line marks PSS sampling time.",
           hjust = 0, vjust = 1, size = 3.3, fontface = "bold") +
  annotate("segment", x = 0.02, xend = 0.16, y = 0.42, yend = 0.42,
           linewidth = 0.55, color = "grey25") +
  annotate("text", x = 0.18, y = 0.42, label = "trajectory", hjust = 0,
           vjust = 0.5, size = 2.8, color = "grey25") +
  annotate("segment", x = 0.49, xend = 0.49, y = 0.24, yend = 0.60,
           linewidth = 0.45, linetype = "22", color = "black") +
  annotate("point", x = 0.56, y = 0.42, size = 1.8, shape = 21,
           fill = "white", color = "grey25", stroke = 0.45) +
  annotate("text", x = 0.58, y = 0.42, label = "PSS measurement", hjust = 0,
           vjust = 0.5, size = 2.8, color = "grey25") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(0, 4, 2, 4))

Fig1a <- ggplot(traj_df, aes(x = time, y = abundance, color = node, group = node)) +
  geom_rect(data = bg_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                fill = bg_fill),
            inherit.aes = FALSE, alpha = 0.9, color = NA, show.legend = FALSE) +
  geom_vline(xintercept = max(traj_df$time), linewidth = 0.35,
             linetype = "22", color = "black", alpha = 0.75) +
  geom_line(linewidth = 0.5, alpha = 0.86, show.legend = FALSE) +
  geom_point(data = pss_df, shape = 21, fill = "white", stroke = 0.5,
             size = 1.8, show.legend = FALSE) +
  geom_text(data = pss_df, aes(x = label_x, label = node), hjust = 0,
            vjust = 0.5, size = 2.35, show.legend = FALSE) +
  facet_grid(model ~ condition) +
  scale_fill_identity() +
  scale_x_continuous(limits = c(0, 22.2), breaks = c(0, 5, 10, 15, 20),
                     expand = expansion(mult = c(0.01, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.08))) +
  labs(
    x = "time after perturbation",
    y = "state abundance",
    title = "ODE dynamics generate perturbed steady-state measurements"
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey95", color = "grey82", linewidth = 0.35),
    strip.text = element_text(face = "bold", size = 8.5),
    axis.title = element_text(size = 8.5),
    axis.text = element_text(size = 7.2, color = "grey25"),
    legend.position = "none",
    panel.spacing = unit(0.9, "lines"),
    plot.margin = margin(5.5, 24, 5.5, 5.5)
  )

# Objects returned for later patchwork assembly; no files are written here.
Fig1a

## ------------------------------------------------------ Fig1b: identifiability ----
# PSS can distinguish nonlinear additive steady-state functions from standard
# linear gLV, but cannot distinguish standard gLV from an additive linear ODE
# sharing the same steady-state equation.
set.seed(43)
N_id <- 180
sigma_id <- 0.003
U_id <- matrix(runif(N_id * n_sp, -0.25, 0.55), N_id, n_sp)
U_id[1, ] <- 0

B_true <- matrix(0, n_sp, n_sp)
B_true[A_true != 0] <- 0.42 * sign(A_true[A_true != 0])
gamma_nl <- gam_true + rowSums(abs(B_true)) + 0.45

deriv_add_quad <- function(t, state, parms) {
  x <- pmax(state, 0)
  fx <- r_true + as.numeric(A_true %*% x) + as.numeric(B_true %*% (x^2)) -
    gamma_nl * x + parms$u
  list(fx)
}

integrate_terminal <- function(u, deriv_fun, gamma_ref = gam_true,
                               t_max = 80) {
  x0 <- pmax(as.numeric(solve(diag(gamma_ref) - A_true, r_true + u)), 0.05)
  out <- ode(y = x0, times = c(0, t_max), func = deriv_fun,
             parms = list(u = u), method = "lsoda", rtol = 1e-9, atol = 1e-11)
  as.numeric(out[nrow(out), -1])
}

steady_additive_linear <- function(U) {
  t(apply(U, 1, ss_linear))
}

steady_glv_linear <- function(U) {
  t(apply(U, 1, function(u) integrate_terminal(u, deriv_glv, gam_true)))
}

steady_additive_quadratic <- function(U) {
  t(apply(U, 1, function(u) integrate_terminal(u, deriv_add_quad, gamma_nl)))
}

basis_gof <- function(U, X) {
  X_obs <- X + matrix(rnorm(length(X), sd = sigma_id), nrow = nrow(X))
  X_obs <- pmax(X_obs, 1e-6)
  X1 <- sweep(X_obs, 2, colMeans(X_obs))
  X2_raw <- cbind(X_obs, X_obs^2)
  X2 <- sweep(X2_raw, 2, colMeans(X2_raw))
  out <- lapply(seq_len(n_sp), function(j) {
    y <- -(U[, j] - mean(U[, j]))
    rss1 <- sum(lm.fit(X1, y)$residuals^2)
    rss2 <- sum(lm.fit(X2, y)$residuals^2)
    bic1 <- N_id * log(rss1 / N_id) + ncol(X1) * log(N_id)
    bic2 <- N_id * log(rss2 / N_id) + ncol(X2) * log(N_id)
    data.frame(
      node = sp_names[j],
      rel_rss_drop = pmax(0, (rss1 - rss2) / rss1),
      bic_gain_per_sample = (bic1 - bic2) / N_id,
      bic_selects_quadratic = bic2 < bic1
    )
  })
  do.call(rbind, out)
}

id_rows <- list(
  "Additive linear ODE" = basis_gof(U_id, steady_additive_linear(U_id)),
  "gLV ODE" = basis_gof(U_id, steady_glv_linear(U_id)),
  "Additive nonlinear ODE" = basis_gof(U_id, steady_additive_quadratic(U_id))
)
id_df <- do.call(rbind, Map(function(d, nm) {
  d$mechanism <- nm
  d
}, id_rows, names(id_rows)))
id_df$mechanism <- factor(id_df$mechanism,
                          levels = c("Additive linear ODE", "gLV ODE",
                                     "Additive nonlinear ODE"))

# Representative steady-state function shapes for the same source effect.
edge_pick <- which(A_true != 0, arr.ind = TRUE)[1, ]
target_pick <- edge_pick[1]
source_pick <- edge_pick[2]
x_grid <- seq(0, 1.25, length.out = 160)
shape_df <- rbind(
  data.frame(
    mechanism = "Additive linear ODE",
    x = x_grid,
    effect = A_true[target_pick, source_pick] * x_grid
  ),
  data.frame(
    mechanism = "gLV ODE",
    x = x_grid,
    effect = A_true[target_pick, source_pick] * x_grid
  ),
  data.frame(
    mechanism = "Additive nonlinear ODE",
    x = x_grid,
    effect = A_true[target_pick, source_pick] * x_grid +
      B_true[target_pick, source_pick] * x_grid^2
  )
)
shape_df$mechanism <- factor(shape_df$mechanism, levels = levels(id_df$mechanism))

Fig1b_function_shape <- ggplot(shape_df,
                               aes(x = x, y = effect, color = mechanism,
                                   linetype = mechanism)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey70") +
  geom_line(linewidth = 0.8, show.legend = FALSE) +
  annotate("text", x = 0.04, y = Inf,
           label = sprintf("representative edge: %s -> %s",
                           sp_names[source_pick], sp_names[target_pick]),
           hjust = 0, vjust = 1.5, size = 2.65, color = "grey25") +
  scale_color_manual(values = c("#6D8B3D", "#3E8EAE", "#B05A8A")) +
  scale_linetype_manual(values = c("solid", "22", "solid")) +
  labs(
    x = "source abundance",
    y = "steady-state effect",
    title = "Steady-state function shape"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 10.5),
    axis.title = element_text(size = 8.4),
    axis.text = element_text(size = 7.3, color = "grey25"),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

Fig1b_identifiability <- ggplot(id_df,
                                aes(x = mechanism, y = bic_gain_per_sample,
                                    fill = mechanism)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey55") +
  geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.72,
               color = "grey30", linewidth = 0.35, show.legend = FALSE) +
  geom_point(aes(group = node), position = position_jitter(width = 0.08, height = 0),
             size = 1.45, alpha = 0.75, color = "grey20", show.legend = FALSE) +
  scale_fill_manual(values = c("#F5F8E8", "#EEF6FA", "#F6EAF2")) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.10))) +
  labs(
    x = NULL,
    y = "BIC gain after adding quadratic basis",
    title = "PSS separates nonlinear steady-state structure"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 10.5),
    axis.text.x = element_text(size = 7.6, angle = 15, hjust = 1),
    axis.text.y = element_text(size = 7.3, color = "grey25"),
    axis.title.y = element_text(size = 8.4),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

if (requireNamespace("patchwork", quietly = TRUE)) {
  Fig1b <- Fig1b_function_shape | Fig1b_identifiability
} else {
  Fig1b <- Fig1b_identifiability
}

Fig1b

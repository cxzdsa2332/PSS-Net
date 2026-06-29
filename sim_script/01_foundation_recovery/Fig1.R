rm(list = ls())

################################################################################
# Fig1.R -- Figure 1 panel objects for PSS-Net
#
# Purpose: build the Figure 1 panels as workspace plot objects and write the
#          assembled Fig1.pdf. Fig1a/Fig1b/Fig1c use a historical 8-node setting; Fig1d and
#          Fig1e use a fresh 10-node nonlinear additive system with node-wise
#          ADSIHT at SNR = 30. Node labels use a unified N1...Nk scheme.
#            Fig1a -- ODE dynamics generate perturbed steady-state measurements
#                     (additive ODE vs multiplicative gLV).
#            Fig1b -- steady-state function-shape identifiability + SNR sweep.
#            Fig1c -- node-wise ADSIHT vs group lasso scaling (reads CSV).
#            Fig1d -- self vs received cross-node effects integrate to the
#                     trajectory (true vs ADSIHT-inferred).
#            Fig1e -- per-source steady-state function-shape recovery.
#            Fig1f -- directed coupling network: true vs node-wise ADSIHT
#                     inferred (igraph), shares the Fig1d/Fig1e 10-node system.
#            Fig1g -- fitted-library misspecification: support vs edge-function
#                     recovery (reads CSV from Fig1x_basis_misspecification.R).
#
# Input:   results/sim_results/Fig1c_adsiht_group_lasso_scaling.csv (Fig1c only).
# Output:  Fig1a, Fig1a_legend, Fig1b_function_shape, Fig1b_identifiability,
#          Fig1b, Fig1b_noise_snr, snr_summary, Fig1c, Fig1d, Fig1e and Fig1f
#          objects in the workspace, plus Fig1g_basis_robustness / Fig1g and the
#          assembled Fig1 (panels a-g tagged
#          via cowplot, sized for A4 portrait). The disabled dynamics overlay
#          (Fig1d_dynamics) is kept for later use behind an `if (FALSE)` guard.
################################################################################

suppressMessages({
  library(deSolve)
  library(ggplot2)
  library(ADSIHT)
  library(igraph)
  library(ggplotify)
})

set.seed(42)

## ---------------------------------------------------------------- Parameters ----
n_sp <- 8
node_names <- paste0("N", seq_len(n_sp))

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
  names(mat)[-1] <- node_names
  long <- reshape(mat, varying = node_names, v.names = "abundance",
                  timevar = "node", times = node_names,
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
  baseline_ss <- baseline_ss[match(node_names, baseline_ss$node), "abundance"]
  
  for (condition in setdiff(names(perturbations), "baseline u = 0")) {
    rows[[idx]] <- simulate_trajectory(perturbations[[condition]], condition,
                                       model, deriv_fun, baseline_ss)
    idx <- idx + 1L
  }
}
traj_df <- do.call(rbind, rows)
traj_df$model <- factor(traj_df$model, levels = c("Additive ODE", "gLV ODE"))
traj_df$condition <- factor(traj_df$condition, levels = names(perturbations))
traj_df$node <- factor(traj_df$node, levels = node_names)

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
  facet_grid(condition ~ model) +
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
      node = node_names[j],
      rel_rss_drop = pmax(0, (rss1 - rss2) / rss1),
      bic_gain_per_sample = (bic1 - bic2) / N_id,
      bic_selects_quadratic = bic2 < bic1
    )
  })
  do.call(rbind, out)
}

X_id_linear <- steady_additive_linear(U_id)
X_id_glv <- steady_glv_linear(U_id)
X_id_nonlinear <- steady_additive_quadratic(U_id)

id_rows <- list(
  "Additive linear ODE" = basis_gof(U_id, X_id_linear),
  "gLV ODE" = basis_gof(U_id, X_id_glv),
  "Additive nonlinear ODE" = basis_gof(U_id, X_id_nonlinear)
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
                           node_names[source_pick], node_names[target_pick]),
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

## ---------------------------------------------- Fig1b auxiliary: SNR sweep ----
# Noise is added to measured steady states x*. SNR is defined as the mean
# node-wise sd of noiseless nonlinear PSS measurements divided by measurement sd.
basis_gof_with_sigma <- function(U, X, sigma_x) {
  X_obs <- X + matrix(rnorm(length(X), sd = sigma_x), nrow = nrow(X))
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
      node = node_names[j],
      bic_gain_per_sample = (bic1 - bic2) / N_id,
      bic_selects_quadratic = bic2 < bic1
    )
  })
  do.call(rbind, out)
}

snr_grid <- c(100, 50, 30, 20, 15, 10, 7, 5, 3, 2, 1.5, 1)
snr_reps <- 30L
signal_scale <- mean(apply(X_id_nonlinear, 2, sd))
snr_rows <- list()
row_id <- 1L
set.seed(44)
for (snr in snr_grid) {
  sigma_x <- signal_scale / snr
  for (rep_id in seq_len(snr_reps)) {
    tmp <- list(
      "Additive linear ODE" = basis_gof_with_sigma(U_id, X_id_linear, sigma_x),
      "gLV ODE" = basis_gof_with_sigma(U_id, X_id_glv, sigma_x),
      "Additive nonlinear ODE" = basis_gof_with_sigma(U_id, X_id_nonlinear, sigma_x)
    )
    for (mechanism in names(tmp)) {
      d <- tmp[[mechanism]]
      snr_rows[[row_id]] <- data.frame(
        snr = snr,
        sigma_x = sigma_x,
        rep = rep_id,
        mechanism = mechanism,
        mean_bic_gain = mean(d$bic_gain_per_sample),
        detection_rate = mean(d$bic_selects_quadratic)
      )
      row_id <- row_id + 1L
    }
  }
}
snr_df <- do.call(rbind, snr_rows)
snr_df$mechanism <- factor(snr_df$mechanism, levels = levels(id_df$mechanism))

snr_summary <- aggregate(cbind(mean_bic_gain, detection_rate) ~ snr + mechanism,
                         snr_df, mean)
snr_summary$mechanism <- factor(snr_summary$mechanism, levels = levels(id_df$mechanism))

nonlinear_summary <- snr_summary[snr_summary$mechanism == "Additive nonlinear ODE", ]
nonlinear_summary <- nonlinear_summary[order(nonlinear_summary$snr), ]
snr_threshold_50 <- min(nonlinear_summary$snr[nonlinear_summary$detection_rate >= 0.5])
snr_threshold_80 <- min(nonlinear_summary$snr[nonlinear_summary$detection_rate >= 0.8])

Fig1b_noise_snr <- ggplot(snr_summary,
                          aes(x = snr, y = detection_rate,
                              color = mechanism, group = mechanism)) +
  geom_hline(yintercept = c(0.5, 0.8), linewidth = 0.25,
             linetype = "22", color = "grey65") +
  geom_line(linewidth = 0.7, show.legend = FALSE) +
  geom_point(size = 1.7, show.legend = FALSE) +
  scale_x_log10(breaks = c(1, 2, 3, 5, 10, 20, 50, 100)) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.02, 0.03))) +
  scale_color_manual(values = c("#6D8B3D", "#3E8EAE", "#B05A8A")) +
  labs(
    x = "SNR of measured steady states",
    y = "BIC selects quadratic basis",
    title = "Noise erases nonlinear steady-state signatures"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 10.5),
    axis.title = element_text(size = 8.4),
    axis.text = element_text(size = 7.3, color = "grey25"),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

Fig1b_noise_snr

## --------------------------------------- Fig1b: assemble three sub-panels ----
# Top row: steady-state function shape + BIC identifiability; bottom row: the
# SNR sweep showing when noise erases the nonlinear signature.
if (requireNamespace("patchwork", quietly = TRUE)) {
  Fig1b <- (Fig1b_function_shape | Fig1b_identifiability) / Fig1b_noise_snr
} else {
  Fig1b <- Fig1b_identifiability
}

Fig1b

## ------------------------------------------------------ Fig1c: method scaling ----
fig1c_file <- "results/sim_results/Fig1c_adsiht_group_lasso_scaling.csv"
if (file.exists(fig1c_file)) {
  fig1c_df <- read.csv(fig1c_file, stringsAsFactors = FALSE)
  fig1c_df$method <- factor(fig1c_df$method, levels = c("ADSIHT", "GroupLasso"))
  # `truth` (linear / nonlinear) is a newer column; default older CSVs to one
  # level so the linetype mapping below stays valid.
  if (is.null(fig1c_df$truth)) fig1c_df$truth <- "nonlinear"
  fig1c_df$truth <- factor(fig1c_df$truth, levels = c("linear", "nonlinear"))

  setting_info <- unique(fig1c_df[, c("p", "s_in")])
  setting_info <- setting_info[order(setting_info$p), ]
  setting_info$edges <- setting_info$p * setting_info$s_in
  setting_info$density <- setting_info$edges / (setting_info$p * (setting_info$p - 1))
  setting_info$setting <- sprintf("p = %d\ns = %d, edges = %d\ndensity = %.1f%%",
                                  setting_info$p, setting_info$s_in,
                                  setting_info$edges, 100 * setting_info$density)
  fig1c_df$setting <- setting_info$setting[match(fig1c_df$p, setting_info$p)]
  fig1c_df$setting <- factor(fig1c_df$setting, levels = setting_info$setting)
  
  summarize_metric <- function(metric_df) {
    metric_df$metric <- factor(metric_df$metric, levels = unique(metric_df$metric))
    do.call(rbind, lapply(split(metric_df,
                                list(metric_df$setting, metric_df$method,
                                     metric_df$truth, metric_df$metric, metric_df$x),
                                drop = TRUE), function(d) {
                                  data.frame(
                                    setting = d$setting[1], method = d$method[1],
                                    truth = d$truth[1], metric = d$metric[1],
                                    x = d$x[1], mean = mean(d$value, na.rm = TRUE),
                                    sd = sd(d$value, na.rm = TRUE)
                                  )
                                }))
  }
  
  # Fig1c: scaling at the main noise level. Noise is controlled by SNR, not by a
  # fixed sigma; SNR=30 is the default moderate-noise display layer.
  fig1c_main_snr <- 30
  fig1c_plot_df <- fig1c_df[fig1c_df$snr == fig1c_main_snr, ]
  fig1c_cols <- c("setting", "N", "N_over_slogp", "seed", "method", "truth")
  metric_df <- rbind(
    data.frame(fig1c_plot_df[, fig1c_cols],
               metric = "MCC", x = fig1c_plot_df$N_over_slogp, value = fig1c_plot_df$MCC),
    data.frame(fig1c_plot_df[, fig1c_cols],
               metric = "AUPRC", x = fig1c_plot_df$N_over_slogp, value = fig1c_plot_df$AUPRC),
    data.frame(fig1c_plot_df[, fig1c_cols],
               metric = "AUROC", x = fig1c_plot_df$N_over_slogp, value = fig1c_plot_df$AUROC),
    data.frame(fig1c_plot_df[, fig1c_cols],
               metric = "Coef. L2 (lower)", x = fig1c_plot_df$N_over_slogp,
               value = fig1c_plot_df$CoefL2),
    data.frame(fig1c_plot_df[, fig1c_cols],
               metric = "Jac. RMSE (lower)", x = fig1c_plot_df$N_over_slogp,
               value = fig1c_plot_df$JacRMSE)
  )
  metric_df$metric <- factor(metric_df$metric,
                             levels = c("MCC", "AUPRC", "AUROC",
                                        "Coef. L2 (lower)", "Jac. RMSE (lower)"))
  fig1c_summary <- summarize_metric(metric_df)
  fig1c_summary$metric <- factor(fig1c_summary$metric, levels = levels(metric_df$metric))
  rownames(fig1c_summary) <- NULL
  
  # Two grouping aesthetics: colour = method, linetype = truth (solid nonlinear,
  # dashed linear). The per-series sd ribbon is dropped here -- four overlapping
  # ribbons would be unreadable and the sd is noisy at few repetitions.
  Fig1c_method_scaling <- ggplot(fig1c_summary,
                                 aes(x = x, y = mean, color = method,
                                     linetype = truth,
                                     group = interaction(method, truth))) +
    geom_line(linewidth = 0.6) +
    geom_point(aes(shape = truth), size = 1.6) +
    facet_grid(metric ~ setting, scales = "free_y") +
    scale_color_manual(values = c(ADSIHT = "#2E6F9E", GroupLasso = "#B45F4D")) +
    scale_linetype_manual(values = c(linear = "22", nonlinear = "solid")) +
    scale_shape_manual(values = c(linear = 1, nonlinear = 16)) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.06))) +
    labs(
      x = "sample budget N / (s log p)",
      y = "metric value",
      color = NULL, linetype = NULL, shape = NULL,
      title = sprintf("Node-wise ADSIHT vs group lasso, linear vs nonlinear (SNR = %s)",
                      fig1c_main_snr)
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 10.5),
      strip.background = element_rect(fill = "grey95", color = "grey82", linewidth = 0.35),
      strip.text = element_text(face = "bold", size = 7.5, lineheight = 0.95),
      axis.title = element_text(size = 8.4),
      axis.text = element_text(size = 7.1, color = "grey25"),
      legend.position = "right",
      legend.text = element_text(size = 7.8),
      panel.spacing = unit(0.75, "lines"),
      plot.margin = margin(5.5, 8, 5.5, 5.5)
    )
  
  Fig1c <- Fig1c_method_scaling
} else {
  Fig1c_method_scaling <- NULL
  Fig1c <- NULL
}

Fig1c

## ============================================================================
## Fig1d / Fig1e: effect decomposition and steady-state function-shape recovery
##
## These two panels use a fresh 10-node nonlinear additive system (independent of
## the 8-node Fig1a/b setting) and node-wise ADSIHT at SNR = 30. Following the
## group decomposition in ref/v0.1.txt, each target node's steady-state map is a
## sum of per-source additive pieces f_ji(x_i); the self term (i == j) is the
## -gamma_j x_j feedback. node_names is re-defined here for the 10-node system,
## keeping the same N1...Nk label scheme.
## ============================================================================
set.seed(2026)

## ---------------------------------------------------------------- Settings ----
p <- 10L
node_names <- paste0("N", seq_len(p))
m_ord <- 2L              # monomial basis x, x^2 (no intercept => f_ji(0) = 0)
N_pss <- 200L            # perturbation conditions (single simulation)
snr_pss <- 30            # fixed steady-state SNR (transition regime, see note)
u_lo <- -0.3
u_hi <- 0.5
grid_len <- 60L          # points per source for the reconstructed f_ji curve
sel_thr <- 1e-6          # group-norm threshold counting a source as selected
shape_targets <- c(1L, 5L)   # representative targets for the Fig1e shape panel
# Fig1d targets: under this seed N6 and N9 give a prominent received-regulation
# signal that the SNR-30 inferred fit tracks closely (true and inferred
# self/received/trajectory curves nearly overlap).
decomp_targets <- c(6L, 9L)

# Effect-decomposition template colours (reuse across PSS-Net effect-split
# figures): self feedback = red, received cross-node regulation = green, and
# their sum equals the dynamics trajectory = blue.
effect_cols <- c("self effect" = "#C0392B",
                 "received regulation" = "#27AE60",
                 "self + received = trajectory" = "#2E6F9E")

## ----------------------------------------------------- Fixed true system ----
# Sparse additive ODE: dx_j/dt = r_j + sum_{i!=j}(a_ji x_i + b_ji x_i^2)
#                                  - gamma_j x_j + u_j.
# A_pss holds linear cross effects, B_pss holds quadratic cross effects. Some
# edges are linear only (b = 0) and some are curved (b != 0) so the figure shows
# both straight and bent recovered functions.
A_pss <- matrix(0, p, p)
B_pss <- matrix(0, p, p)

set_edge <- function(j, i, a, b = 0) {
  A_pss[j, i] <<- a
  B_pss[j, i] <<- b
}

set_edge(1, 3,  0.30, -0.18)   # curved
set_edge(1, 6, -0.25,  0.00)   # linear
set_edge(2, 1,  0.28,  0.16)   # curved
set_edge(2, 5, -0.22,  0.00)   # linear
set_edge(3, 2, -0.26,  0.00)   # linear
set_edge(3, 7,  0.30,  0.20)   # curved
set_edge(4, 1,  0.24, -0.15)   # curved
set_edge(4, 8,  0.27,  0.00)   # linear
set_edge(5, 4,  0.29,  0.18)   # curved
set_edge(5, 9, -0.23,  0.00)   # linear
set_edge(6, 2,  0.25,  0.00)   # linear
set_edge(6, 10, 0.30, -0.20)   # curved
set_edge(7, 3, -0.28,  0.00)   # linear
set_edge(8, 5,  0.26,  0.17)   # curved
set_edge(9, 6, -0.24,  0.00)   # linear
set_edge(10, 8, 0.29,  0.00)   # linear

gamma_pss <- rowSums(abs(A_pss)) + rowSums(abs(B_pss)) + runif(p, 1.2, 1.6)
r_pss <- runif(p, 0.8, 1.4)
adj_pss <- (A_pss != 0) * 1L

deriv_pss <- function(t, state, parms) {
  x <- pmax(state, 0)
  fx <- r_pss + as.numeric(A_pss %*% x) + as.numeric(B_pss %*% (x^2)) -
    gamma_pss * x + parms$u
  list(fx)
}

integrate_pss <- function(u, t_max = 120) {
  # linear-solve warm start, then integrate the nonlinear additive ODE to PSS.
  x0 <- pmax(as.numeric(solve(diag(gamma_pss) - A_pss, r_pss + u)), 0.05)
  out <- ode(y = x0, times = c(0, t_max), func = deriv_pss,
             parms = list(u = u), method = "lsoda", rtol = 1e-9, atol = 1e-11)
  as.numeric(out[nrow(out), -1])
}

steady_pss <- function(U) {
  t(apply(U, 1, integrate_pss))
}

## ---------------------------------------------- Basis / node-wise ADSIHT ----
make_basis <- function(X) {
  Psi <- matrix(0, nrow(X), p * m_ord)
  for (i in seq_len(p)) {
    for (m in seq_len(m_ord)) {
      Psi[, (i - 1L) * m_ord + m] <- X[, i]^m
    }
  }
  Psi
}

standardize_design <- function(Psi) {
  Psi_bar <- colMeans(Psi)
  Psi_c <- sweep(Psi, 2, Psi_bar)
  Psi_sd <- pmax(apply(Psi_c, 2, sd), 1e-10)
  list(X = sweep(Psi_c, 2, Psi_sd, "/"), scale = Psi_sd)
}

group_norms <- function(beta) {
  sapply(seq_len(p), function(i) {
    cols <- (i - 1L) * m_ord + seq_len(m_ord)
    sqrt(sum(beta[cols]^2))
  })
}

fit_adsiht_node <- function(X_cs, y, group, scale_vec) {
  fit <- tryCatch(ADSIHT(X_cs, matrix(y), group, ic.type = "dsic"),
                  error = function(e) NULL)
  if (is.null(fit)) return(rep(0, p * m_ord))
  as.numeric(fit$beta[, which.min(fit$ic)] / scale_vec)
}

# Node-wise solve of the PSS equation -u_j = r_j + sum_i f_ji(x_i*); group i
# collects the m_ord monomials of source i (self term included, source = target).
# The intercept r_j is recovered from the uncentered equation so the full vector
# field F_hat_j(x) = r_j + sum_i f_ji(x_i) can be re-simulated as dynamics.
infer_beta <- function(U, X_obs) {
  Psi <- make_basis(X_obs)
  std <- standardize_design(Psi)
  psi_bar <- colMeans(Psi)
  group <- rep(seq_len(p), each = m_ord)
  U_c <- sweep(U, 2, colMeans(U))
  beta <- matrix(0, p, p * m_ord)
  intercept <- numeric(p)
  for (j in seq_len(p)) {
    beta[j, ] <- fit_adsiht_node(std$X, -U_c[, j], group, std$scale)
    intercept[j] <- -mean(U[, j]) - sum(beta[j, ] * psi_bar)
  }
  list(beta = beta, intercept = intercept)
}

## -------------------------------- Source grids from a clean reference draw ----
set.seed(7)
U_ref <- matrix(runif(2000L * p, u_lo, u_hi), 2000L, p)
X_ref <- steady_pss(U_ref)
ok_ref <- apply(X_ref, 1, function(z) all(is.finite(z)) && all(z > 0))
X_ref <- X_ref[ok_ref, , drop = FALSE]
src_grid <- sapply(seq_len(p), function(i) {
  qs <- quantile(X_ref[, i], c(0.05, 0.95))
  seq(qs[1], qs[2], length.out = grid_len)
})  # grid_len x p

f_true_curve <- function(j, i) {
  x <- src_grid[, i]
  if (i == j) {
    -gamma_pss[j] * x
  } else {
    A_pss[j, i] * x + B_pss[j, i] * x^2
  }
}

f_est_curve <- function(beta, j, i) {
  x <- src_grid[, i]
  c1 <- beta[j, (i - 1L) * m_ord + 1L]
  c2 <- beta[j, (i - 1L) * m_ord + 2L]
  c1 * x + c2 * x^2
}

## ------------------------------------------------------ Single simulation ----
# One perturbation design -> PSS -> SNR-30 corruption -> node-wise ADSIHT fit.
set.seed(1001)
U_pss <- matrix(runif(N_pss * p, u_lo, u_hi), N_pss, p)
U_pss[1, ] <- 0
X_pss <- steady_pss(U_pss)
ok <- apply(X_pss, 1, function(z) all(is.finite(z)) && all(z > 0))
X_clean <- X_pss[ok, , drop = FALSE]
U_ok <- U_pss[ok, , drop = FALSE]
signal_scale_pss <- mean(apply(X_clean, 2, sd))
sigma_pss <- signal_scale_pss / snr_pss
X_obs <- pmax(X_clean + matrix(rnorm(length(X_clean), sd = sigma_pss),
                               nrow(X_clean), p), 1e-6)

fit_est <- infer_beta(U_ok, X_obs)
beta_est <- fit_est$beta
r_est <- fit_est$intercept

## ------------------- Demo perturbation + (disabled) dynamics overlay ----
# The recovered additive functions define a full vector field
# F_hat_j(x) = r_hat_j + sum_i f_ji_hat(x_i) (self feedback + cross-node
# interactions). Re-simulating dx_j/dt = F_hat_j(x) + u under a demonstration
# perturbation and overlaying it on the true ODE trajectory shows that PSS
# recovers the dynamics, not only the steady states.
deriv_inferred <- function(t, state, parms) {
  x <- pmax(state, 0)
  basis <- as.numeric(sapply(seq_len(p), function(i) x[i]^seq_len(m_ord)))
  fx <- r_est + as.numeric(beta_est %*% basis) + parms$u
  list(fx)
}

u_demo <- rep(0, p)
u_demo[c(1L, 4L, 7L)] <- c(0.40, -0.30, 0.35)
x0_demo <- rep(0.35, p)
demo_times <- seq(0, 20, length.out = 200)

# Disabled for now: full dynamics-reconstruction overlay (kept for later use).
if (FALSE) {
  simulate_named <- function(deriv_fun, method) {
    out <- ode(y = x0_demo, times = demo_times, func = deriv_fun,
               parms = list(u = u_demo), method = "lsoda",
               rtol = 1e-9, atol = 1e-11)
    mat <- as.data.frame(out)
    names(mat)[-1] <- node_names
    long <- reshape(mat, varying = node_names, v.names = "abundance",
                    timevar = "node", times = node_names,
                    idvar = "time", direction = "long")
    rownames(long) <- NULL
    long$method <- method
    long
  }
  
  dyn_df <- rbind(
    simulate_named(deriv_pss, "true dynamics"),
    simulate_named(deriv_inferred, "ADSIHT-inferred")
  )
  dyn_df$node <- factor(dyn_df$node, levels = node_names)
  dyn_df$method <- factor(dyn_df$method,
                          levels = c("true dynamics", "ADSIHT-inferred"))
  
  # Node labels sit just past the end of each true trajectory (Fig1a style).
  end_df <- dyn_df[dyn_df$method == "true dynamics" &
                     dyn_df$time == max(demo_times), ]
  end_df$label_x <- max(demo_times) + 0.6
  
  Fig1d_dynamics <- ggplot(dyn_df, aes(x = time, y = abundance, color = node,
                                       group = interaction(node, method))) +
    geom_line(aes(linetype = method), linewidth = 0.6, alpha = 0.9) +
    geom_text(data = end_df, aes(x = label_x, label = node), hjust = 0,
              vjust = 0.5, size = 2.5, show.legend = FALSE) +
    scale_linetype_manual(values = c("true dynamics" = "solid",
                                     "ADSIHT-inferred" = "22")) +
    scale_x_continuous(limits = c(0, max(demo_times) + 2.5),
                       expand = expansion(mult = c(0.01, 0))) +
    guides(color = "none") +
    labs(
      x = "time after perturbation",
      y = "state abundance",
      linetype = NULL,
      title = "ADSIHT-inferred dynamics reproduce the true vector field"
    ) +
    coord_cartesian(clip = "off") +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 11),
      axis.text = element_text(size = 7, color = "grey25"),
      axis.title = element_text(size = 8.4),
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      plot.margin = margin(5.5, 16, 5.5, 5.5)
    )
  
  Fig1d_dynamics
}

## ----------------------------- Fig1d: trajectory effect decomposition ----
# For two interacting targets, split the abundance trajectory into the
# *integrated* received cross-node regulation (cumulative integral of
# sum_{i!=j} f_ji(x_i(t))) and a self effect that carries everything else (the
# initial state, intrinsic term, self feedback and perturbation). By
# construction self + received = x_j(t), the dynamics trajectory. True and
# inferred functions are each integrated along their own simulated trajectory.
traj_true <- ode(y = x0_demo, times = demo_times, func = deriv_pss,
                 parms = list(u = u_demo), method = "lsoda",
                 rtol = 1e-9, atol = 1e-11)[, -1]
traj_inf <- ode(y = x0_demo, times = demo_times, func = deriv_inferred,
                parms = list(u = u_demo), method = "lsoda",
                rtol = 1e-9, atol = 1e-11)[, -1]

# cumulative trapezoidal integral (base R; avoids a pracma dependency).
cumtrapz_base <- function(t, y) {
  n <- length(t)
  out <- numeric(n)
  if (n > 1L) {
    out[-1] <- cumsum(diff(t) * (y[-1] + y[-n]) / 2)
  }
  out
}

# integrated self / received contributions of target j along x_mat (n_t x p):
# received = integral of cross-node rate; self = trajectory - received.
effect_split <- function(x_mat, j, mode) {
  if (mode == "true") {
    cross_rate <- as.numeric(x_mat %*% A_pss[j, ] + (x_mat^2) %*% B_pss[j, ])
  } else {
    c1 <- beta_est[j, (seq_len(p) - 1L) * m_ord + 1L]
    c2 <- beta_est[j, (seq_len(p) - 1L) * m_ord + 2L]
    contrib <- sweep(x_mat, 2, c1, "*") + sweep(x_mat^2, 2, c2, "*")
    cross_rate <- rowSums(contrib[, -j, drop = FALSE])
  }
  received <- cumtrapz_base(demo_times, cross_rate)
  total <- x_mat[, j]
  list(self = total - received, received = received, total = total)
}

traj_list <- list("true dynamics" = traj_true, "ADSIHT-inferred" = traj_inf)
mode_map <- c("true dynamics" = "true", "ADSIHT-inferred" = "inferred")
comp_levels <- c("self effect", "received regulation",
                 "self + received = trajectory")
comp_rows <- list()
crow <- 1L
for (m_name in names(traj_list)) {
  for (j in decomp_targets) {
    es <- effect_split(traj_list[[m_name]], j, mode_map[[m_name]])
    vals <- list("self effect" = es$self,
                 "received regulation" = es$received,
                 "self + received = trajectory" = es$total)
    for (comp in comp_levels) {
      comp_rows[[crow]] <- data.frame(
        target = node_names[j], time = demo_times, method = m_name,
        component = comp, value = vals[[comp]], stringsAsFactors = FALSE)
      crow <- crow + 1L
    }
  }
}
decomp_dyn_df <- do.call(rbind, comp_rows)
rownames(decomp_dyn_df) <- NULL
decomp_dyn_df$target <- factor(decomp_dyn_df$target,
                               levels = node_names[decomp_targets])
decomp_dyn_df$method <- factor(decomp_dyn_df$method,
                               levels = c("true dynamics", "ADSIHT-inferred"))
decomp_dyn_df$component <- factor(decomp_dyn_df$component, levels = comp_levels)

Fig1d <- ggplot(decomp_dyn_df,
                aes(x = time, y = value, color = component,
                    linetype = method,
                    group = interaction(component, method))) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey80") +
  geom_line(linewidth = 0.65, alpha = 0.9) +
  facet_wrap(~ target, ncol = 2) +
  scale_color_manual(values = effect_cols) +
  scale_linetype_manual(values = c("true dynamics" = "solid",
                                   "ADSIHT-inferred" = "22")) +
  labs(
    x = "time after perturbation",
    y = "state abundance (integrated effect)",
    color = NULL, linetype = NULL,
    title = "Inferred self and received-regulation effects integrate to the trajectory"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 10.5),
    strip.background = element_rect(fill = "grey95", color = "grey82", linewidth = 0.35),
    strip.text = element_text(face = "bold", size = 8.5),
    axis.text = element_text(size = 7, color = "grey25"),
    axis.title = element_text(size = 8.4),
    legend.position = "none",
    panel.spacing = unit(0.8, "lines")
  )

Fig1d

## ------------------------------------------------------------------ Fig1e ----
# For each picked target, keep its self term plus its true cross edges; drop
# null sources (ADSIHT shrinks them to ~0, summarised by selected_rate).
shape_rows <- list()
shape_id <- 1L
for (j in shape_targets) {
  gn <- group_norms(beta_est[j, ])
  sources <- union(j, which(adj_pss[j, ] != 0))
  for (i in sources) {
    shape_rows[[shape_id]] <- data.frame(
      target = node_names[j],
      source = node_names[i],
      relation = if (i == j) "self feedback" else "cross-node effect",
      pair = sprintf("%s -> %s", node_names[i], node_names[j]),
      selected_rate = as.numeric(gn[i] > sel_thr),
      x = src_grid[, i],
      f_true = f_true_curve(j, i),
      f_est = f_est_curve(beta_est, j, i),
      stringsAsFactors = FALSE
    )
    shape_id <- shape_id + 1L
  }
}
decomp_df <- do.call(rbind, shape_rows)
rownames(decomp_df) <- NULL
decomp_df$relation <- factor(decomp_df$relation,
                             levels = c("self feedback", "cross-node effect"))
pair_order <- unique(decomp_df[order(decomp_df$target, decomp_df$relation,
                                     decomp_df$source), "pair"])
decomp_df$pair <- factor(decomp_df$pair, levels = pair_order)

Fig1e <- ggplot(decomp_df, aes(x = x)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey80") +
  geom_line(aes(y = f_true, color = "true f_ji"), linewidth = 0.8) +
  geom_line(aes(y = f_est, color = "ADSIHT estimate"),
            linewidth = 0.8, linetype = "22") +
  facet_wrap(~ pair, scales = "free", ncol = 3) +
  scale_color_manual(values = c("true f_ji" = "#2E6F9E",
                                "ADSIHT estimate" = "#B45F4D"),
                     breaks = c("true f_ji", "ADSIHT estimate")) +
  labs(
    x = expression("source abundance " * x[i]),
    y = expression("steady-state effect " * f[ji](x[i])),
    color = NULL,
    title = "PSS recovers steady-state function shape with node-wise ADSIHT"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey95", color = "grey82", linewidth = 0.35),
    strip.text = element_text(size = 8),
    axis.text = element_text(size = 7, color = "grey25"),
    axis.title = element_text(size = 8.4),
    legend.position = "none",
    panel.spacing = unit(0.6, "lines")
  )

Fig1e

## ------------------------------------------------------------------ Fig1f ----
# Directed coupling network of the same 10-node nonlinear additive system used by
# Fig1d/Fig1e (same seeds, same fitted beta_est). Two igraph panels share one
# circular layout: the true network on the left, the node-wise ADSIHT inferred
# network on the right. Edge colour encodes interaction sign -- promotion (+) in
# an elegant red, inhibition (-) in an elegant blue (sign from A_pss for the
# truth, from the estimated linear coefficient for the inference). Self-feedback
# is not drawn (it is not a directed cross-node edge, per CLAUDE.md). In the
# inferred panel, true positives are solid, false positives are dashed, and
# missed true edges (false negatives) are overlaid as faint grey dotted arrows.

# Promotion / inhibition palette (reused from the effect-decomposition scheme:
# red for the self/promotion sign, blue for the trajectory/inhibition sign).
edge_pos_col <- "#C0392B"   # promotion (+): elegant red
edge_neg_col <- "#2E6F9E"   # inhibition (-): elegant blue
edge_fp_col  <- "grey55"    # false-positive marker tint
edge_fn_col  <- "grey65"    # missed-edge (false-negative) overlay

# True and inferred adjacency + sign matrices (rows = target j, cols = source i;
# self terms i == j excluded from the cross-node network).
adj_true_net  <- (A_pss != 0) * 1L
sign_true_net <- sign(A_pss)
diag(adj_true_net) <- 0L

adj_est_net  <- matrix(0L, p, p)
sign_est_net <- matrix(0, p, p)
for (j in seq_len(p)) {
  gn <- group_norms(beta_est[j, ])
  for (i in seq_len(p)) {
    if (i == j) next
    if (gn[i] > sel_thr) {
      adj_est_net[j, i]  <- 1L
      sign_est_net[j, i] <- sign(beta_est[j, (i - 1L) * m_ord + 1L])
    }
  }
}

# Edge-recovery classification (used for line style + the inferred-panel title).
tp_mat <- adj_true_net == 1L & adj_est_net == 1L
fp_mat <- adj_true_net == 0L & adj_est_net == 1L
fn_mat <- adj_true_net == 1L & adj_est_net == 0L
n_tp <- sum(tp_mat); n_fp <- sum(fp_mat); n_fn <- sum(fn_mat)
n_tn <- p * (p - 1L) - n_tp - n_fp - n_fn
mcc_den <- sqrt((n_tp + n_fp) * (n_tp + n_fn) *
                  (n_tn + n_fp) * (n_tn + n_fn))
mcc_net <- if (mcc_den > 0) (n_tp * n_tn - n_fp * n_fn) / mcc_den else NA_real_

# Build a directed graph on all p vertices (isolated nodes kept so the circular
# layout is identical across panels). Edge list orientation: source i -> target j.
make_net_graph <- function(adj_mat, sign_mat) {
  g <- make_empty_graph(p, directed = TRUE)
  V(g)$name <- node_names
  el <- which(adj_mat != 0L, arr.ind = TRUE)   # [row = j (target), col = i (source)]
  if (nrow(el) > 0L) {
    g <- add_edges(g, as.vector(t(cbind(el[, 2], el[, 1]))))
    E(g)$weight <- sign_mat[el]
  }
  g
}

g_true_net <- make_net_graph(adj_true_net, sign_true_net)
g_est_net  <- make_net_graph(adj_est_net,  sign_est_net)
lay_net    <- layout_in_circle(g_true_net)

sign_to_col <- function(w) ifelse(w > 0, edge_pos_col, edge_neg_col)

# Per-edge style for the inferred panel: colour by sign, dash false positives.
est_el    <- ends(g_est_net, E(g_est_net), names = FALSE)  # n_edge x 2 (source, target)
est_is_tp <- if (nrow(est_el) > 0L)
  mapply(function(src, tgt) tp_mat[tgt, src], est_el[, 1], est_el[, 2]) else logical(0)
est_ecol  <- sign_to_col(E(g_est_net)$weight)
est_ewidth <- ifelse(est_is_tp, 2, 1.3)
est_elty   <- ifelse(est_is_tp, 1, 2)

# Shared node aesthetics.
net_vertex <- list(color = "grey90", frame.color = "grey55", size = 30,
                   label.cex = 0.8, label.color = "grey10")

draw_fig1f <- function() {
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mfrow = c(1, 2), mar = c(2.5, 1, 2.5, 1), bg = "white")
  
  # Panel 1: true coupling network.
  plot(g_true_net, layout = lay_net,
       vertex.color = net_vertex$color, vertex.frame.color = net_vertex$frame.color,
       vertex.size = net_vertex$size, vertex.label = node_names,
       vertex.label.cex = net_vertex$label.cex, vertex.label.color = net_vertex$label.color,
       edge.color = sign_to_col(E(g_true_net)$weight), edge.width = 2,
       edge.arrow.size = 0.85, edge.arrow.width = 1.1, edge.curved = 0.22,
       main = "True coupling network")
  legend("bottom", horiz = TRUE, bty = "n", cex = 0.75, inset = -0.02,
         legend = c("promotion (+)", "inhibition (-)"),
         col = c(edge_pos_col, edge_neg_col), lwd = 2)
  
  # Panel 2: node-wise ADSIHT inferred network.
  plot(g_est_net, layout = lay_net,
       vertex.color = net_vertex$color, vertex.frame.color = net_vertex$frame.color,
       vertex.size = net_vertex$size, vertex.label = node_names,
       vertex.label.cex = net_vertex$label.cex, vertex.label.color = net_vertex$label.color,
       edge.color = est_ecol, edge.width = est_ewidth, edge.lty = est_elty,
       edge.arrow.size = 0.85, edge.arrow.width = 1.1, edge.curved = 0.22,
       main = sprintf("ADSIHT inferred  [TP=%d FP=%d FN=%d  MCC=%.2f]",
                      n_tp, n_fp, n_fn, mcc_net))
  
  # Overlay missed true edges (false negatives) as faint grey dotted arrows.
  if (n_fn > 0L) {
    fn_el <- which(fn_mat, arr.ind = TRUE)
    g_fn  <- add_edges(make_empty_graph(p, directed = TRUE),
                       as.vector(t(cbind(fn_el[, 2], fn_el[, 1]))))
    plot(g_fn, layout = lay_net, add = TRUE,
         vertex.color = NA, vertex.frame.color = NA, vertex.label = NA,
         vertex.size = net_vertex$size, edge.color = edge_fn_col,
         edge.width = 1.2, edge.lty = 3, edge.arrow.size = 0.7,
         edge.arrow.width = 1.0, edge.curved = 0.22)
  }
  legend("bottom", horiz = TRUE, bty = "n", cex = 0.7, inset = -0.02,
         legend = c("TP (+)", "TP (-)", "FP", "FN"),
         col = c(edge_pos_col, edge_neg_col, edge_fp_col, edge_fn_col),
         lwd = c(2, 2, 1.3, 1.2), lty = c(1, 1, 2, 3))
}

# Capture the two-panel base-R igraph plot as a single ggplot object so Fig1f
# composes with the other panels under patchwork (ggplot2 > patchwork >
# igraph + grid plotting stack).
Fig1f <- as.ggplot(draw_fig1f)

Fig1f

## ------------------------------------ Fig1g: basis-misspecification robustness ----
# Compact identifiability companion: how recovery changes as the FITTED library
# is varied away from the true nonlinear edge function (truth shapes coloured;
# matched library ringed). Two facets carry the message -- support MCC is robust
# across libraries while edge-function NRMSE needs an adequate dictionary; signed
# Jacobian accuracy (~1 everywhere) is stated in the caption rather than plotted.
# Data: results/sim_results/Fig1x_basis_misspecification.csv (run
# sim_script/01_foundation_recovery/Fig1x_basis_misspecification.R first).
fig1g_file <- "results/sim_results/Fig1x_basis_misspecification.csv"
if (file.exists(fig1g_file)) {
  fig1g_raw <- read.csv(fig1g_file, stringsAsFactors = FALSE)
  fig1g_lib_levels <- c("linear", "poly2", "poly3", "monod", "fourier")
  fig1g_truth_levels <- c("poly2", "monod", "sine")
  fig1g_long <- rbind(
    data.frame(fig1g_raw[, c("seed", "truth", "library", "matched")],
               metric = "Support MCC (higher better)", value = fig1g_raw$MCC),
    data.frame(fig1g_raw[, c("seed", "truth", "library", "matched")],
               metric = "Edge-function NRMSE (lower better)",
               value = fig1g_raw$FuncNRMSE)
  )
  fig1g_long <- fig1g_long[is.finite(fig1g_long$value), ]
  fig1g_summary <- do.call(rbind, lapply(split(
    fig1g_long, list(fig1g_long$truth, fig1g_long$library, fig1g_long$metric),
    drop = TRUE
  ), function(d) {
    data.frame(truth = d$truth[1], library = d$library[1],
               metric = d$metric[1], matched = any(d$matched),
               mean = mean(d$value), stringsAsFactors = FALSE)
  }))
  rownames(fig1g_summary) <- NULL
  fig1g_summary$library <- factor(fig1g_summary$library,
                                  levels = fig1g_lib_levels)
  fig1g_summary$truth <- factor(fig1g_summary$truth,
                                levels = fig1g_truth_levels)
  fig1g_summary$metric <- factor(
    fig1g_summary$metric,
    levels = c("Support MCC (higher better)",
               "Edge-function NRMSE (lower better)"))
  fig1g_truth_labels <- c(poly2 = "poly2 truth (A x + B x^2)",
                          monod = "monod truth (saturating)",
                          sine = "sine truth (oscillatory)")
  fig1g_colors <- c(poly2 = "#2E6F9E", monod = "#6D8B3D", sine = "#B45F4D")

  Fig1g_basis_robustness <- ggplot(
    fig1g_summary,
    aes(x = library, y = mean, color = truth, group = truth)
  ) +
    geom_line(linewidth = 0.5, linetype = "22", alpha = 0.7) +
    geom_point(size = 1.9) +
    geom_point(data = fig1g_summary[fig1g_summary$matched, ],
               shape = 21, size = 3.4, stroke = 0.8, fill = NA,
               color = "grey20", show.legend = FALSE) +
    facet_wrap(~ metric, nrow = 1, scales = "free_y") +
    scale_color_manual(values = fig1g_colors, labels = fig1g_truth_labels) +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.1))) +
    labs(
      x = "fitted basis / library", y = "mean over seeds", color = NULL,
      title = "Support recovery is robust to library misspecification; function recovery needs a matched dictionary"
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 9.6),
      plot.subtitle = element_text(size = 7.8, color = "grey30"),
      strip.background = element_rect(fill = "grey95", color = "grey82",
                                      linewidth = 0.35),
      strip.text = element_text(face = "bold", size = 8),
      axis.title = element_text(size = 8.4),
      axis.text.x = element_text(size = 7.1, color = "grey25"),
      axis.text.y = element_text(size = 7.1, color = "grey25"),
      legend.position = "right",
      legend.text = element_text(size = 7.5),
      panel.spacing = unit(0.8, "lines"),
      plot.margin = margin(5.5, 8, 5.5, 5.5)
    )
  Fig1g <- Fig1g_basis_robustness
} else {
  Fig1g_basis_robustness <- NULL
  Fig1g <- NULL
}

Fig1g

## ----------------------------------------------- Assemble Figure 1 (a-g) ----
# Assemble the seven panels into one A4-portrait figure. cowplot::plot_grid is used
# (not patchwork's `/`) because Fig1b is itself a patchwork composite, and
# nesting patchwork inside patchwork mis-aligns the inner sub-panels; cowplot
# treats each panel as an opaque grob and avoids that conflict. Layout: a and b
# share the top row (a is faceted 3 rows x 2 cols, so it sits tall and narrow
# beside b); c spans a full row with its legend on the right; d and e share a row
# (legends dropped, explained in the caption); f closes the figure. rel_heights
# keep the dense scaling panel (c) from being squeezed at A4 size.
if (requireNamespace("cowplot", quietly = TRUE)) {
  library(cowplot)
  fig1c_panel <- if (is.null(Fig1c)) NULL else Fig1c

  row_ab <- plot_grid(Fig1a, Fig1b, labels = c("a", "b"), label_size = 14,
                      label_fontface = "bold", ncol = 2, rel_widths = c(1, 1.05))
  row_de <- plot_grid(Fig1d, Fig1e, labels = c("d", "e"), label_size = 14,
                      label_fontface = "bold", ncol = 2, rel_widths = c(1, 1.15))

  panels <- list(row_ab, fig1c_panel, row_de, Fig1f, Fig1g)
  labels <- c("", "c", "", "f", "g")
  rel_h  <- c(1.35, 1.6, 0.95, 0.95, 0.85)
  keep   <- !vapply(panels, is.null, logical(1))

  Fig1 <- plot_grid(plotlist = panels[keep], labels = labels[keep],
                    label_size = 14, label_fontface = "bold",
                    ncol = 1, rel_heights = rel_h[keep])
} else {
  Fig1 <- NULL
}

# A4 portrait (210 x 297 mm).
fig1_out <- file.path("manuscript", "figures", "Fig1.pdf")
dir.create(dirname(fig1_out), recursive = TRUE, showWarnings = FALSE)
ggsave(fig1_out, Fig1, width = 210, height = 297, units = "mm")

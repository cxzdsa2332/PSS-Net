rm(list = ls())

################################################################################
# Fig2.R -- Figure 2 panel objects for PSS-Net
#
# Purpose: build Figure 2 scaling/design panel objects from existing simulation
#          result files.
#            Fig2a -- high-dimensional sample complexity of node-wise ADSIHT
#                     under sparse PSS systems.
#            Fig2b -- optimal perturbation design in a linear steady-state
#                     system: random vs maximin vs sequential D-optimal.
#
# Input:   results/sim_results/Fig2a_highdim_sample_complexity.csv
#          results/sim_results/Fig2b_design_linear.csv
# Output:  Fig2a_highdim_sample_complexity, Fig2a, Fig2b_design_linear, Fig2b
#          and the Fig2b_gain companion (delta-MCC vs random) in the workspace.
################################################################################

suppressMessages({
  library(ggplot2)
})

## ------------------------------------------------ Fig2a: sample complexity ----
fig2a_file <- "results/sim_results/Fig2a_highdim_sample_complexity.csv"
if (!file.exists(fig2a_file)) {
  stop("Missing ", fig2a_file,
       ". Run sim_script/02_scaling_design/Fig2a_highdim_sample_complexity.R first.")
}

fig2a_raw <- read.csv(fig2a_file, stringsAsFactors = FALSE)
required_cols <- c("p", "s", "N", "N_over_slogp", "seed", "Pr", "Re", "MCC")
missing_cols <- setdiff(required_cols, names(fig2a_raw))
if (length(missing_cols) > 0L) {
  stop("Missing required columns in Fig2a_highdim_sample_complexity.csv: ",
       paste(missing_cols, collapse = ", "))
}

fig2a_info <- unique(fig2a_raw[, c("p", "s")])
fig2a_info <- fig2a_info[order(fig2a_info$p), ]
fig2a_info$edges <- fig2a_info$p * fig2a_info$s
fig2a_info$density <- fig2a_info$edges / (fig2a_info$p * (fig2a_info$p - 1))
fig2a_info$setting <- sprintf("p = %d, s = %d, edges = %d, density = %.1f%%",
                              fig2a_info$p, fig2a_info$s, fig2a_info$edges,
                              100 * fig2a_info$density)
fig2a_raw$setting <- fig2a_info$setting[match(fig2a_raw$p, fig2a_info$p)]
fig2a_raw$setting <- factor(fig2a_raw$setting, levels = fig2a_info$setting)
fig2a_raw$p_label <- factor(paste0("p = ", fig2a_raw$p),
                            levels = paste0("p = ", fig2a_info$p))

fig2a_metric_df <- data.frame(
  fig2a_raw[, c("p", "p_label", "setting", "s", "N", "N_over_slogp", "seed")],
  value = fig2a_raw$MCC
)

fig2a_summary <- do.call(rbind, lapply(split(
  fig2a_metric_df,
  list(fig2a_metric_df$p_label, fig2a_metric_df$setting,
       fig2a_metric_df$N_over_slogp),
  drop = TRUE
), function(d) {
  data.frame(
    p = d$p[1],
    p_label = d$p_label[1],
    setting = d$setting[1],
    N_over_slogp = d$N_over_slogp[1],
    N = d$N[1],
    mean = mean(d$value, na.rm = TRUE),
    sd = sd(d$value, na.rm = TRUE),
    n_seed = length(unique(d$seed)),
    stringsAsFactors = FALSE
  )
}))
rownames(fig2a_summary) <- NULL
fig2a_summary$p_label <- factor(fig2a_summary$p_label,
                                levels = paste0("p = ", fig2a_info$p))
fig2a_summary$setting <- factor(fig2a_summary$setting,
                                levels = fig2a_info$setting)

threshold_mcc <- 0.8
threshold_rows <- lapply(split(fig2a_summary, fig2a_summary$p_label), function(d) {
  d <- d[order(d$N_over_slogp), ]
  hit <- which(d$mean >= threshold_mcc)
  if (length(hit) == 0L) {
    return(data.frame(
      p = d$p[1],
      p_label = d$p_label[1],
      threshold_x = NA_real_,
      label_x = max(d$N_over_slogp),
      label_y = min(0.96, max(d$mean) + 0.08),
      label = sprintf("%s: MCC < %.1f at max budget", d$p_label[1], threshold_mcc),
      stringsAsFactors = FALSE
    ))
  }
  k <- hit[1]
  if (k == 1L) {
    x_thr <- d$N_over_slogp[k]
  } else {
    x0 <- d$N_over_slogp[k - 1L]; y0 <- d$mean[k - 1L]
    x1 <- d$N_over_slogp[k]; y1 <- d$mean[k]
    x_thr <- x0 + (threshold_mcc - y0) * (x1 - x0) / (y1 - y0)
  }
  data.frame(
    p = d$p[1],
    p_label = d$p_label[1],
    threshold_x = x_thr,
    label_x = x_thr,
    label_y = threshold_mcc + 0.08,
    label = sprintf("%s: %.1f", d$p_label[1], x_thr),
    stringsAsFactors = FALSE
  )
})
fig2a_threshold_df <- do.call(rbind, threshold_rows)
fig2a_threshold_df$p_label <- factor(fig2a_threshold_df$p_label,
                                     levels = paste0("p = ", fig2a_info$p))

fig2a_label_df <- fig2a_info
fig2a_label_df$x <- 1.05
fig2a_label_df$y <- c(0.95, 0.85, 0.75)[seq_len(nrow(fig2a_label_df))]
fig2a_label_df$p_label <- factor(paste0("p = ", fig2a_label_df$p),
                                 levels = paste0("p = ", fig2a_info$p))

Fig2a_highdim_sample_complexity <- ggplot(
  fig2a_summary,
  aes(x = N_over_slogp, y = mean, color = p_label, fill = p_label,
      group = p_label)
) +
  geom_vline(xintercept = c(1, 2, 5), linewidth = 0.28,
             linetype = "22", color = "grey72") +
  geom_hline(yintercept = threshold_mcc, linewidth = 0.35,
             linetype = "33", color = "grey45") +
  geom_ribbon(aes(ymin = pmax(0, mean - sd), ymax = pmin(1, mean + sd)),
              alpha = 0.13, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 1.9) +
  geom_segment(data = fig2a_threshold_df[!is.na(fig2a_threshold_df$threshold_x), ],
               aes(x = threshold_x, xend = threshold_x, y = 0,
                   yend = threshold_mcc, color = p_label),
               inherit.aes = FALSE, linewidth = 0.45, linetype = "42",
               show.legend = FALSE) +
  geom_label(data = fig2a_threshold_df,
             aes(x = label_x, y = label_y, label = label, color = p_label),
             inherit.aes = FALSE, fill = "white", linewidth = 0.18,
             size = 2.35, label.padding = unit(0.12, "lines"),
             show.legend = FALSE) +
  geom_text(data = fig2a_label_df,
            aes(x = x, y = y, label = setting, color = p_label),
            inherit.aes = FALSE, hjust = 0, vjust = 0.5,
            size = 2.45, lineheight = 0.95, show.legend = FALSE) +
  scale_color_manual(values = c("p = 8" = "#2E6F9E",
                                "p = 50" = "#6D8B3D",
                                "p = 100" = "#B45F4D")) +
  scale_fill_manual(values = c("p = 8" = "#2E6F9E",
                               "p = 50" = "#6D8B3D",
                               "p = 100" = "#B45F4D")) +
  scale_x_continuous(breaks = c(1, 2, 3, 5, 8, 10, 12, 16, 20, 24, 30),
                     expand = expansion(mult = c(0.03, 0.08))) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.03, 0.04))) +
  labs(
    x = expression("sample budget " * N / (s * log(p))),
    y = "MCC",
    color = NULL,
    title = "Node-wise ADSIHT sample complexity follows the sparse budget scale",
    subtitle = "Vertical guides mark N = 1, 2, and 5 x s log(p); labels mark the interpolated budget for MCC > 0.8"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 10.8),
    plot.subtitle = element_text(size = 8, color = "grey30"),
    axis.title = element_text(size = 8.5),
    axis.text = element_text(size = 7.2, color = "grey25"),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    panel.spacing = unit(0.9, "lines"),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

Fig2a <- Fig2a_highdim_sample_complexity

Fig2a

## ----------------------------------------------- Fig2b: perturbation design ----
# Linear 8-node additive-ODE benchmark, node-wise ADSIHT, monomial M = 2, noise
# at steady-state SNR = 30 (aligned with Fig1c/Fig2a). Three design strategies
# are compared on the same rescaled budget axis N / (s log p) as Fig2a: random
# (baseline), maximin space-filling, and sequential D-optimal (the method). The
# design gain is read as a leftward shift -- the budget needed to reach a target
# MCC -- and is quantified as the sample-savings of D-optimal over random.
fig2b_file <- "results/sim_results/Fig2b_design_linear.csv"
if (!file.exists(fig2b_file)) {
  stop("Missing ", fig2b_file,
       ". Run sim_script/02_scaling_design/Fig2b_design_linear.R first.")
}

fig2b_raw <- read.csv(fig2b_file, stringsAsFactors = FALSE)
fig2b_required <- c("strategy", "N", "MCC", "seed")
fig2b_missing <- setdiff(fig2b_required, names(fig2b_raw))
if (length(fig2b_missing) > 0L) {
  stop("Missing required columns in Fig2b_design_linear.csv: ",
       paste(fig2b_missing, collapse = ", "))
}
# Rescaled budget axis, matching Fig2a. Fall back to N / (s log p) if absent.
if (is.null(fig2b_raw$N_over_slogp)) {
  fig2b_raw$N_over_slogp <- fig2b_raw$N / (fig2b_raw$s * log(fig2b_raw$p))
}

fig2b_labels <- c(random = "Random",
                  maximin = "Maximin (space-filling)",
                  dopt = "Sequential D-optimal")
fig2b_cols_raw <- c(random = "#9AA0A6", maximin = "#6D8B3D", dopt = "#2E6F9E")
fig2b_cols <- setNames(fig2b_cols_raw, fig2b_labels[names(fig2b_cols_raw)])

fig2b_summary <- do.call(rbind, lapply(split(
  fig2b_raw, list(fig2b_raw$strategy, fig2b_raw$N_over_slogp), drop = TRUE
), function(d) {
  data.frame(
    strat = d$strategy[1],
    x = d$N_over_slogp[1],
    N = d$N[1],
    mean = mean(d$MCC, na.rm = TRUE),
    se = sd(d$MCC, na.rm = TRUE) / sqrt(sum(!is.na(d$MCC))),
    n_seed = length(unique(d$seed)),
    stringsAsFactors = FALSE
  )
}))
rownames(fig2b_summary) <- NULL
fig2b_summary$strategy <- factor(fig2b_labels[fig2b_summary$strat],
                                 levels = fig2b_labels)

# Interpolated budget at which each strategy first reaches the target MCC; the
# horizontal gap between strategies is the sample-savings (design gain).
fig2b_target <- 0.8
cross_budget <- function(d) {
  d <- d[order(d$x), ]
  hit <- which(d$mean >= fig2b_target)
  if (length(hit) == 0L) return(NA_real_)
  k <- hit[1]
  if (k == 1L) return(d$x[1])
  x0 <- d$x[k - 1L]; y0 <- d$mean[k - 1L]
  x1 <- d$x[k]; y1 <- d$mean[k]
  x0 + (fig2b_target - y0) * (x1 - x0) / (y1 - y0)
}
fig2b_cross <- do.call(rbind, lapply(split(fig2b_summary, fig2b_summary$strat),
                                     function(d) {
  data.frame(strat = d$strat[1], strategy = d$strategy[1],
             cross_x = cross_budget(d), stringsAsFactors = FALSE)
}))
rownames(fig2b_cross) <- NULL
x_rand <- fig2b_cross$cross_x[fig2b_cross$strat == "random"]
x_dopt <- fig2b_cross$cross_x[fig2b_cross$strat == "dopt"]
x_mmin <- fig2b_cross$cross_x[fig2b_cross$strat == "maximin"]
gain_dopt <- if (length(x_rand) && length(x_dopt) && !is.na(x_rand) && !is.na(x_dopt))
  100 * (x_rand - x_dopt) / x_rand else NA_real_
fig2b_gain_text <- if (!is.na(gain_dopt)) sprintf(
  "D-optimal reaches MCC %.1f at ~%.0f%% fewer conditions than random (%.1f vs %.1f x s log p)",
  fig2b_target, gain_dopt, x_dopt, x_rand) else
  sprintf("no strategy reaches MCC %.1f within the budget range", fig2b_target)

Fig2b_design_linear <- ggplot(
  fig2b_summary,
  aes(x = x, y = mean, color = strategy, fill = strategy, group = strategy)
) +
  geom_hline(yintercept = fig2b_target, linewidth = 0.35,
             linetype = "33", color = "grey45") +
  geom_vline(xintercept = c(1, 2, 5), linewidth = 0.28,
             linetype = "22", color = "grey78") +
  geom_segment(data = fig2b_cross[!is.na(fig2b_cross$cross_x), ],
               aes(x = cross_x, xend = cross_x, y = 0, yend = fig2b_target,
                   color = strategy),
               inherit.aes = FALSE, linewidth = 0.45, linetype = "42",
               show.legend = FALSE) +
  geom_ribbon(aes(ymin = pmax(0, mean - se), ymax = pmin(1, mean + se)),
              alpha = 0.14, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 1.9) +
  scale_color_manual(values = fig2b_cols) +
  scale_fill_manual(values = fig2b_cols, guide = "none") +
  scale_x_continuous(breaks = c(1, 2, 3, 5, 8),
                     expand = expansion(mult = c(0.03, 0.05))) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.03, 0.04))) +
  labs(
    x = expression("sample budget " * N / (s * log(p))),
    y = "MCC",
    color = NULL,
    title = "Active perturbation design lowers the budget for network recovery",
    subtitle = fig2b_gain_text
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 10.8),
    plot.subtitle = element_text(size = 8, color = "grey30"),
    axis.title = element_text(size = 8.5),
    axis.text = element_text(size = 7.2, color = "grey25"),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

Fig2b <- Fig2b_design_linear

Fig2b

## ------------------------------- Fig2b companion: design gain over random ----
# Non-curve view of the same result: MCC improvement of each structured design
# relative to random at matched budget. This makes the (modest) gain legible and
# shows it concentrates at small budgets and decays once data is plentiful.
fig2b_gain_rows <- do.call(rbind, lapply(
  split(fig2b_summary, fig2b_summary$x), function(d) {
    base <- d$mean[d$strat == "random"]
    if (length(base) == 0L) return(NULL)
    out <- d[d$strat != "random", c("strat", "strategy", "x")]
    out$delta <- d$mean[match(out$strat, d$strat)] - base
    out
  }))
rownames(fig2b_gain_rows) <- NULL

fig2b_gain_rows$x_label <- factor(sprintf("%.1f", fig2b_gain_rows$x),
                                  levels = sprintf("%.1f", sort(unique(fig2b_gain_rows$x))))

Fig2b_gain <- ggplot(fig2b_gain_rows,
                     aes(x = x_label, y = delta, fill = strategy)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey55") +
  geom_col(position = position_dodge(width = 0.7), width = 0.62) +
  scale_fill_manual(values = fig2b_cols) +
  labs(
    x = expression("sample budget " * N / (s * log(p))),
    y = expression(Delta * "MCC vs random"),
    fill = NULL,
    title = "Design gain is largest when data are scarce"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 10.8),
    axis.title = element_text(size = 8.5),
    axis.text = element_text(size = 7.2, color = "grey25"),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

Fig2b_gain

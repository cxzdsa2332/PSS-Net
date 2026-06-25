rm(list = ls())

################################################################################
# Fig2.R -- Figure 2 panel objects for PSS-Net
#
# Purpose: build Figure 2 scaling/design panel objects from existing simulation
#          result files.
#            Fig2a -- high-dimensional sample complexity of node-wise PSS-Net
#                     under sparse PSS systems.
#            Fig2b -- conceptual 3D input-space design mechanism.
#            Fig2c -- regime-by-budget map of MCC gain over random (D-optimal
#                     shown as the oracle upper bound; Fig2e gives the discount).
#            Fig2d -- budget required to reach target recovery levels (oracle
#                     D-optimal upper bound).
#            Fig2e -- pilot-informed exact D-optimal augmentation and oracle
#                     regret under a fair total experimental budget.
#
# Note: panels b--e deliberately separate mechanism, performance gain, sample
#       savings and feasibility instead of repeating four learning curves.
#
# Input:   results/sim_results/Fig2a_highdim_sample_complexity.csv
#          results/sim_results/Fig2b_design_linear.csv
#          results/sim_results/design_nl_comparison.csv             (Fig2c)
#          results/sim_results/design_nl_seq_comparison.csv         (Fig2d)
#          results/sim_results/Fig2e_oracle_vs_estimated_design.csv (Fig2e)
# Output:  Fig2a..Fig2e panel objects, descriptive aliases, and assembled Fig2.
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
    title = "Node-wise PSS-Net sample complexity follows the sparse budget scale"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 10.8),
    plot.subtitle = element_text(size = 8, color = "grey30"),
    axis.title = element_text(size = 8.5),
    axis.text = element_text(size = 7.2, color = "grey25"),
    legend.position = "none",
    panel.spacing = unit(0.9, "lines"),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

Fig2a <- Fig2a_highdim_sample_complexity

Fig2a

## ------------------------------- Fig2b: perturbation-design mechanism ----
# Load the conceptual u-space schematic in an isolated environment because the
# standalone script intentionally clears its own workspace.
concept_env <- new.env(parent = globalenv())
sys.source("sim_script/02_scaling_design/Fig2_explain.R", envir = concept_env)
Fig2b_design_mechanism <- concept_env$Fig2_design_concept
Fig2b <- Fig2b_design_mechanism

Fig2b

## ------------------------------------ shared strategy labels / colours (c-d) ----
# One palette so the quantitative design panels read as one comparison.
fig2_labels <- c(random = "Random",
                 maximin = "Maximin (space-filling)",
                 dopt = "D-optimal (oracle)")
fig2_cols_raw <- c(random = "#9AA0A6", maximin = "#6D8B3D", dopt = "#2E6F9E")
fig2_cols <- setNames(fig2_cols_raw, fig2_labels[names(fig2_cols_raw)])

# p = 8, s = 2 design benchmark; same N / (s log p) budget axis as Fig2a.
fig2_design_p <- 8L
fig2_design_s <- 2L
fig2_slogp <- fig2_design_s * log(fig2_design_p)

# Load one design-simulation CSV and retain seed-level results for paired
# comparisons against random at the same system seed and budget.
load_design <- function(file) {
  if (!file.exists(file)) {
    stop("Missing ", file,
         ". Run the matching sim_script/02_scaling_design design sim first.")
  }
  raw <- read.csv(file, stringsAsFactors = FALSE)
  required <- c("strategy", "N", "MCC", "seed")
  missing <- setdiff(required, names(raw))
  if (length(missing) > 0L) {
    stop("Missing required columns in ", file, ": ",
         paste(missing, collapse = ", "))
  }
  raw$x <- raw$N / fig2_slogp
  raw
}

# Per (strategy, N) mean and SE of MCC for budget-to-target calculations.
summarize_design <- function(raw) {
  summ <- do.call(rbind, lapply(split(
    raw, list(raw$strategy, raw$N), drop = TRUE
  ), function(d) {
    data.frame(
      strat = d$strategy[1],
      x = d$x[1],
      N = d$N[1],
      mean = mean(d$MCC, na.rm = TRUE),
      se = sd(d$MCC, na.rm = TRUE) / sqrt(sum(!is.na(d$MCC))),
      n_seed = length(unique(d$seed)),
      stringsAsFactors = FALSE
    )
  }))
  rownames(summ) <- NULL
  summ$strategy <- factor(fig2_labels[summ$strat], levels = fig2_labels)
  summ[order(summ$strat, summ$x), ]
}

## ---------------- Fig2c: where structured perturbation design adds value ----
# Compress the three former budget curves into one regime-by-budget map. Each
# tile is the mean MCC gain over random at the same budget; a dot marks the
# better structured strategy within a regime and budget.
regime_summaries <- list(
  "Linear" = load_design("results/sim_results/Fig2b_design_linear.csv"),
  "Nonlinear" = load_design("results/sim_results/design_nl_comparison.csv"),
  "Strong nonlinear" = load_design("results/sim_results/design_nl_seq_comparison.csv")
)
regime_raw <- regime_summaries
regime_summaries <- lapply(regime_raw, summarize_design)

fig2c_gain <- do.call(rbind, lapply(names(regime_raw), function(regime) {
  d <- regime_raw[[regime]]
  base <- d[d$strategy == "random", c("seed", "N", "MCC")]
  names(base)[3] <- "MCC_random"
  paired <- merge(d[d$strategy != "random", ], base,
                  by = c("seed", "N"), all = FALSE)
  paired$delta_i <- paired$MCC - paired$MCC_random
  out <- do.call(rbind, lapply(split(
    paired, list(paired$strategy, paired$N), drop = TRUE
  ), function(z) {
    data.frame(
      strat = z$strategy[1], N = z$N[1], x = z$x[1],
      delta = mean(z$delta_i, na.rm = TRUE),
      delta_se = sd(z$delta_i, na.rm = TRUE) / sqrt(sum(!is.na(z$delta_i))),
      win_rate = mean(z$delta_i > 0, na.rm = TRUE),
      n_pair = sum(!is.na(z$delta_i)), stringsAsFactors = FALSE
    )
  }))
  out$regime <- regime
  out
}))
rownames(fig2c_gain) <- NULL
# The single dopt strategy here uses the true response map -> label it as the
# oracle upper bound, consistent with Fig2e (oracle_dopt vs pilot_dopt).
fig2c_gain$method <- ifelse(fig2c_gain$strat == "dopt",
                            "D-optimal (oracle)", "Maximin")
fig2c_gain$method <- factor(fig2c_gain$method,
                            levels = c("Maximin", "D-optimal (oracle)"))
fig2c_gain$regime <- factor(fig2c_gain$regime,
                            levels = c("Linear", "Nonlinear", "Strong nonlinear"))
fig2c_gain$best <- ave(fig2c_gain$delta,
                       interaction(fig2c_gain$regime, fig2c_gain$x),
                       FUN = function(z) z == max(z))
fig2c_gain$cell_label <- sprintf("%+.2f\n%.0f%%", fig2c_gain$delta,
                                 100 * fig2c_gain$win_rate)

Fig2c_design_gain_map <- ggplot(
  fig2c_gain, aes(x = factor(N), y = method, fill = delta)
) +
  geom_tile(width = 0.72, height = 0.82, color = "white", linewidth = 0.45) +
  geom_text(aes(label = cell_label), size = 2.15, lineheight = 0.88,
            color = "grey15") +
  geom_point(data = fig2c_gain[fig2c_gain$best, ],
             shape = 21, size = 1.5, stroke = 0.35, fill = "white",
             color = "grey15", position = position_nudge(y = 0.27)) +
  scale_fill_gradient2(low = "#B45F4D", mid = "#F7F7F5", high = "#2E6F9E",
                       midpoint = 0, name = expression(Delta * "MCC\nvs random")) +
  scale_x_discrete(expand = expansion(mult = c(0.03, 0.04))) +
  facet_grid(. ~ regime, scales = "free_x", space = "free_x") +
  labs(
    x = "total perturbation conditions N", y = NULL,
    title = "Designed perturbations are compared directly with random inputs"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 10.8),
    plot.subtitle = element_text(size = 7.8, color = "grey35"),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95", color = "grey82"),
    strip.text = element_text(face = "bold", size = 8),
    axis.text.y = element_text(size = 7.4, color = "grey20"),
    axis.text.x = element_text(size = 6.8),
    legend.position = "right",
    legend.title = element_text(size = 7.5),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

Fig2c <- Fig2c_design_gain_map

Fig2c

## ------------------------- Fig2d: experimental budget saved at target MCC ----
threshold_budget <- function(s, targets = c(0.5, 0.6)) {
  do.call(rbind, lapply(split(s, s$strat), function(d) {
    d <- d[order(d$x), ]
    do.call(rbind, lapply(targets, function(target) {
      hit <- which(d$mean >= target)
      x_req <- NA_real_
      censored <- FALSE
      if (length(hit) > 0L) {
        k <- hit[1]
        if (k == 1L) x_req <- d$x[k] else {
          x0 <- d$x[k - 1L]; x1 <- d$x[k]
          y0 <- d$mean[k - 1L]; y1 <- d$mean[k]
          x_req <- x0 + (target - y0) * (x1 - x0) / (y1 - y0)
        }
      } else {
        x_req <- max(d$x) + 0.8
        censored <- TRUE
      }
      data.frame(strat = d$strat[1], strategy = d$strategy[1],
                 target = target, x_required = x_req, censored = censored)
    }))
  }))
}

fig2d_budget <- do.call(rbind, lapply(names(regime_summaries), function(regime) {
  out <- threshold_budget(regime_summaries[[regime]])
  out$regime <- regime
  out
}))
rownames(fig2d_budget) <- NULL
fig2d_budget$regime <- factor(fig2d_budget$regime,
                              levels = c("Linear", "Nonlinear", "Strong nonlinear"))
fig2d_budget$target_label <- factor(sprintf("Target MCC = %.1f", fig2d_budget$target),
                                    levels = c("Target MCC = 0.5", "Target MCC = 0.6"))
fig2d_random <- fig2d_budget[fig2d_budget$strat == "random",
                             c("regime", "target_label", "x_required")]
names(fig2d_random)[3] <- "x_random"
fig2d_budget <- merge(fig2d_budget, fig2d_random,
                      by = c("regime", "target_label"), all.x = TRUE)

Fig2d_budget_to_target <- ggplot(
  fig2d_budget, aes(y = regime, x = x_required, color = strategy)
) +
  geom_segment(data = fig2d_budget[fig2d_budget$strat != "random", ],
               aes(x = x_random, xend = x_required, yend = regime),
               color = "grey75", linewidth = 0.65,
               arrow = arrow(length = unit(1.5, "mm"), type = "closed")) +
  geom_point(aes(shape = censored), size = 2.4,
             position = position_dodge(width = 0.42)) +
  geom_text(data = fig2d_budget[fig2d_budget$censored, ],
            aes(label = "> max"), nudge_y = -0.18, size = 2.3,
            show.legend = FALSE) +
  facet_wrap(~ target_label, nrow = 1, scales = "free_x") +
  scale_color_manual(values = fig2_cols) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 24), guide = "none") +
  labs(
    x = expression("required budget " * N / (s * log(p))), y = NULL,
    color = NULL,
    title = "Designed inputs change the budget required for network recovery"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 10.8),
    plot.subtitle = element_text(size = 7.8, color = "grey35"),
    strip.background = element_rect(fill = "grey95", color = "grey82"),
    strip.text = element_text(face = "bold", size = 8),
    axis.text = element_text(size = 7.2, color = "grey25"),
    legend.position = "bottom",
    legend.text = element_text(size = 7.5),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

Fig2d <- Fig2d_budget_to_target

Fig2d

## ---------------- Fig2e: package-backed pilot-informed D-optimal ----------
fig2e_file <- "results/sim_results/Fig2e_oracle_vs_estimated_design.csv"
if (!file.exists(fig2e_file)) {
  stop("Missing ", fig2e_file,
       ". Run sim_script/02_scaling_design/Fig2e_oracle_vs_estimated_design.R first.")
}
fig2e_raw <- read.csv(fig2e_file, stringsAsFactors = FALSE)
fig2e_required <- c("seed", "pilot_n", "N_total", "strategy", "MCC")
fig2e_missing <- setdiff(fig2e_required, names(fig2e_raw))
if (length(fig2e_missing) > 0L) {
  stop("Missing required Fig2e columns: ", paste(fig2e_missing, collapse = ", "))
}

fig2e_summary <- aggregate(MCC ~ pilot_n + N_total + strategy,
                           fig2e_raw, mean, na.rm = TRUE)
fig2e_wide <- reshape(fig2e_summary, idvar = c("pilot_n", "N_total"),
                      timevar = "strategy", direction = "wide")
fig2e_wide$oracle_regret <- fig2e_wide$MCC.oracle_dopt - fig2e_wide$MCC.pilot_dopt
fig2e_wide$gain_vs_random <- fig2e_wide$MCC.pilot_dopt - fig2e_wide$MCC.random
fig2e_wide$gain_vs_maximin <- fig2e_wide$MCC.pilot_dopt - fig2e_wide$MCC.maximin
fig2e_wide$label <- sprintf("%.2f\n(%+.2f)", fig2e_wide$oracle_regret,
                            fig2e_wide$gain_vs_random)
fig2e_wide$beats_random <- fig2e_wide$gain_vs_random > 0

Fig2e_regret <- ggplot(
  fig2e_wide, aes(x = factor(N_total), y = factor(pilot_n), fill = oracle_regret)
) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = label), size = 2.55, lineheight = 0.9, color = "grey15") +
  geom_point(aes(shape = beats_random), size = 1.8,
             position = position_nudge(x = 0.34, y = 0.31), color = "grey15") +
  scale_fill_gradient2(low = "#2E6F9E", mid = "#F7F7F5", high = "#B45F4D",
                       midpoint = 0, name = "Oracle regret") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 4),
                     labels = c(`TRUE` = "Pilot D-opt > random",
                                `FALSE` = "Pilot D-opt <= random"),
                     name = NULL) +
  labs(
    x = "Total experimental budget (pilot included)",
    y = "Pilot conditions",
    title = "Pilot data determine how much oracle design value is recoverable"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 10.2),
    plot.subtitle = element_text(size = 7.6, color = "grey35"),
    panel.grid = element_blank(),
    axis.text = element_text(size = 7.2, color = "grey25"),
    legend.position = "bottom",
    legend.title = element_text(size = 7.5),
    legend.text = element_text(size = 7.2),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

design_boxes <- data.frame(
  x = c(1, 3, 5, 7), y = 1,
  label = c("Shared random\npilot PSS", "Estimate local\nsteady-state map",
            "Build candidate\nPSS features", "AlgDesign exact\nD-opt augmentation")
)
design_arrows <- data.frame(x = c(1.75, 3.75, 5.75), xend = c(2.25, 4.25, 6.25),
                          y = 1, yend = 1)
Fig2e_design_flow <- ggplot() +
  geom_label(data = design_boxes, aes(x = x, y = y, label = label),
             fill = c("#EEF2F5", "#F5EFE5", "#EAF0E2", "#E8EFF5"),
             color = "grey20", linewidth = 0.25, size = 3.0,
             label.padding = unit(0.45, "lines")) +
  geom_segment(data = design_arrows,
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(2.2, "mm"), type = "closed"),
               linewidth = 0.55, color = "grey40") +
  coord_cartesian(xlim = c(0.2, 7.8), ylim = c(0.65, 1.45), clip = "off") +
  labs(title = "Two-stage design: PSS mapping + existing exact D-optimal backend") +
  theme_void(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 10.2),
        plot.margin = margin(5.5, 8, 0, 5.5))

if (requireNamespace("patchwork", quietly = TRUE)) {
  Fig2e_oracle_vs_estimated_design <-
    Fig2e_design_flow / Fig2e_regret + patchwork::plot_layout(heights = c(0.36, 1))
} else {
  Fig2e_oracle_vs_estimated_design <- Fig2e_regret
}

Fig2e <- Fig2e_oracle_vs_estimated_design

Fig2e

## ------------------------------------------------------- assembled Figure 2 ----
# Panel descriptions live in the separate caption file
# sim_script/02_scaling_design/Fig2_Caption.md; the caption is intentionally kept
# out of the rendered figure (subtitles were also removed from the panels).
if (requireNamespace("patchwork", quietly = TRUE)) {
  Fig2 <- patchwork::wrap_elements(Fig2a) /
    patchwork::wrap_elements(Fig2b) /
    (patchwork::wrap_elements(Fig2c) | patchwork::wrap_elements(Fig2d)) /
    patchwork::wrap_elements(Fig2e) +
    patchwork::plot_layout(heights = c(1.0, 0.72, 0.9, 1.15)) +
    patchwork::plot_annotation(tag_levels = "a")
} else {
  Fig2 <- NULL
}

Fig2

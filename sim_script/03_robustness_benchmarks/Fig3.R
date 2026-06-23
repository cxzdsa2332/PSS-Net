rm(list = ls())

################################################################################
# Fig3.R -- Figure 3 panel objects for PSS-Net
#
# Purpose: build Figure 3 robustness/benchmark panel objects.
#            Fig3a -- method capability matrix (positioning panel; no data file,
#                     hand-curated from ref/external_benchmark_methods.md).
#            Fig3e -- compositional-data limitations for PSS network recovery
#                     (historically named Fig3c; object name kept for now).
#          Benchmark panels Fig3b/Fig3c read the external-benchmark CSV and are
#          built in analysis_script (see simulation_plan.md).
#
# Input:   results/sim_results/Fig3c_compositional_data_limitation.csv
# Output:  Fig3a_method_capability_matrix / Fig3a, and
#          Fig3c_compositional_data_limitation / Fig3c objects in workspace.
################################################################################

suppressMessages({
  library(ggplot2)
})

## ------------------------------------- Fig3a: method capability matrix ----
# Qualitative positioning panel: which capabilities each method class has, so the
# unique PSS-Net cell (directed + uses u + nonlinear edge functions on steady-state
# data) is visible. Hand-curated from ref/external_benchmark_methods.md.
# Encoding: 1 = yes, 0.5 = partial / conditional, 0 = no.
fig3a_methods <- c(
  "PSS-Net (ADSIHT)", "MRA", "Lasso / Elastic net", "Group lasso",
  "PySINDy STLSQ (PSS)", "GENIE3 / GRNBoost2", "GIES (interventional)",
  "Correlation / partial cor.", "SPIEC-EASI / SparCC"
)
fig3a_caps <- c("Uses\nperturbation u", "Directed\nedges", "Signed\nedges",
                "Nonlinear\nedge functions", "Built-in sparse\nselection",
                "Steady-state\n(no time series)", "Compositional\naware")
# Rows follow fig3a_methods; columns follow fig3a_caps. "Built-in sparse
# selection" = produces a sparse edge set directly, without a post-hoc threshold
# (so the dense MRA / correlation baselines are 0; they need thresholding).
fig3a_values <- rbind(
  c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0),  # PSS-Net (group + within-group sparsity)
  c(1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0),  # MRA (dense OLS + t-test threshold)
  c(1.0, 1.0, 1.0, 0.5, 1.0, 1.0, 0.0),  # Lasso / EN (nonlinear via library)
  c(1.0, 1.0, 1.0, 0.5, 1.0, 1.0, 0.0),  # Group lasso
  c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0),  # PySINDy STLSQ on steady-state PSS library
  c(0.0, 1.0, 0.0, 1.0, 0.5, 1.0, 0.5),  # GENIE3 (state-only, ranked importances)
  c(1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0),  # GIES (interventional, DAG, linear-Gaussian)
  c(0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0),  # Correlation / partial cor. (dense)
  c(0.0, 0.0, 0.5, 0.0, 0.5, 1.0, 1.0)   # SPIEC-EASI (sparse) / SparCC (dense)
)

fig3a_df <- data.frame(
  method = factor(rep(fig3a_methods, times = length(fig3a_caps)),
                  levels = rev(fig3a_methods)),
  capability = factor(rep(fig3a_caps, each = length(fig3a_methods)),
                      levels = fig3a_caps),
  value = as.vector(fig3a_values)
)
fig3a_df$glyph <- c("·", "~", "✓")[match(fig3a_df$value, c(0, 0.5, 1))]
fig3a_df$val_f <- factor(fig3a_df$value, levels = c(0, 0.5, 1))

# Highlight the PSS-Net row (top of the reversed factor).
fig3a_hi <- length(fig3a_methods)

Fig3a_method_capability_matrix <- ggplot(
  fig3a_df, aes(x = capability, y = method, fill = val_f)
) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = glyph), size = 3.1, color = "grey15") +
  annotate("rect", xmin = 0.5, xmax = length(fig3a_caps) + 0.5,
           ymin = fig3a_hi - 0.5, ymax = fig3a_hi + 0.5,
           fill = NA, color = "#2E6F9E", linewidth = 0.9) +
  scale_fill_manual(values = c("0" = "#EEEEEE", "0.5" = "#F3D9B0",
                               "1" = "#BCD3B6"), guide = "none") +
  scale_x_discrete(position = "top") +
  labs(
    x = NULL, y = NULL,
    title = "PSS-Net occupies a unique cell among network-inference methods"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 10.8),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 6.8, color = "grey20", lineheight = 0.9),
    axis.text.y = element_text(size = 7.6, color = "grey20"),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

Fig3a <- Fig3a_method_capability_matrix

Fig3a

## ----------------------------------------- Fig3c: compositional limitation ----
fig3c_file <- "results/sim_results/Fig3c_compositional_data_limitation.csv"
if (!file.exists(fig3c_file)) {
  stop("Missing ", fig3c_file,
       ". Run sim_script/03_robustness_benchmarks/Fig3c_compositional_data_limitation.R first.")
}

fig3c_raw <- read.csv(fig3c_file, stringsAsFactors = FALSE)
required_cols <- c("seed", "input", "n_eff", "Pr", "Re", "MCC")
missing_cols <- setdiff(required_cols, names(fig3c_raw))
if (length(missing_cols) > 0L) {
  stop("Missing required columns in Fig3c_compositional_data_limitation.csv: ",
       paste(missing_cols, collapse = ", "))
}

input_levels <- c("abs", "rel", "clr", "rel_x_T")
input_labels <- c(
  abs = "Absolute\nabundance",
  rel = "Relative\nabundance",
  clr = "CLR\ntransform",
  rel_x_T = "Relative x\nnoisy total"
)
fig3c_raw$input <- factor(fig3c_raw$input, levels = input_levels)

fig3c_metric_df <- rbind(
  data.frame(fig3c_raw[, c("seed", "input", "n_eff")],
             metric = "MCC", value = fig3c_raw$MCC),
  data.frame(fig3c_raw[, c("seed", "input", "n_eff")],
             metric = "Precision", value = fig3c_raw$Pr),
  data.frame(fig3c_raw[, c("seed", "input", "n_eff")],
             metric = "Recall", value = fig3c_raw$Re)
)
fig3c_metric_df$metric <- factor(fig3c_metric_df$metric,
                                 levels = c("MCC", "Precision", "Recall"))

fig3c_summary <- do.call(rbind, lapply(split(
  fig3c_metric_df,
  list(fig3c_metric_df$input, fig3c_metric_df$metric),
  drop = TRUE
), function(d) {
  data.frame(
    input = d$input[1],
    metric = d$metric[1],
    mean = mean(d$value, na.rm = TRUE),
    sd = sd(d$value, na.rm = TRUE),
    n_seed = length(unique(d$seed)),
    stringsAsFactors = FALSE
  )
}))
rownames(fig3c_summary) <- NULL
fig3c_summary$input <- factor(fig3c_summary$input, levels = input_levels)
fig3c_summary$metric <- factor(fig3c_summary$metric,
                               levels = c("MCC", "Precision", "Recall"))

fig3c_mcc <- fig3c_summary[fig3c_summary$metric == "MCC", ]
oracle_mcc <- fig3c_mcc$mean[fig3c_mcc$input == "abs"]
fig3c_mcc$loss_vs_abs <- oracle_mcc - fig3c_mcc$mean

Fig3c_compositional_data_limitation <- ggplot(
  fig3c_summary,
  aes(x = input, y = mean, fill = input)
) +
  geom_hline(data = fig3c_summary[fig3c_summary$metric == "MCC", ],
             aes(yintercept = oracle_mcc), inherit.aes = FALSE,
             linewidth = 0.3, linetype = "22", color = "grey45") +
  geom_col(width = 0.62, color = "grey30", linewidth = 0.25, show.legend = FALSE) +
  geom_errorbar(aes(ymin = pmax(0, mean - sd), ymax = pmin(1, mean + sd)),
                width = 0.16, linewidth = 0.35, color = "grey20") +
  geom_point(data = fig3c_metric_df,
             aes(x = input, y = value),
             inherit.aes = FALSE, position = position_jitter(width = 0.08, height = 0),
             size = 1.25, alpha = 0.48, color = "grey15") +
  facet_wrap(~ metric, nrow = 1) +
  scale_x_discrete(labels = input_labels) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.02, 0.05))) +
  scale_fill_manual(values = c(
    abs = "#2E6F9E",
    rel = "#B45F4D",
    clr = "#8B6BAE",
    rel_x_T = "#6D8B3D"
  )) +
  labs(
    x = NULL,
    y = "network recovery",
    title = "Compositional measurements degrade PSS network recovery",
    subtitle = "Dashed line marks oracle absolute-abundance MCC; points are simulation seeds"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 10.8),
    plot.subtitle = element_text(size = 8, color = "grey30"),
    strip.background = element_rect(fill = "grey95", color = "grey82", linewidth = 0.35),
    strip.text = element_text(face = "bold", size = 8.5),
    axis.text.x = element_text(size = 7.2, color = "grey25"),
    axis.text.y = element_text(size = 7.2, color = "grey25"),
    axis.title = element_text(size = 8.5),
    panel.spacing = unit(0.9, "lines"),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

Fig3c <- Fig3c_compositional_data_limitation

Fig3c

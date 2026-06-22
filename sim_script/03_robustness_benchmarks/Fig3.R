rm(list = ls())

################################################################################
# Fig3.R -- Figure 3 panel objects for PSS-Net
#
# Purpose: build Figure 3 robustness/benchmark panel objects from existing
#          simulation result files. Current implementation focuses on Fig3c:
#          compositional-data limitations for PSS network recovery.
#
# Input:   results/sim_results/Fig3c_compositional_data_limitation.csv
# Output:  Fig3c_compositional_data_limitation and Fig3c objects in workspace.
################################################################################

suppressMessages({
  library(ggplot2)
})

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

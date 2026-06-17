################################################################################
# plot_design_curves.R
#
# Purpose: Read perturbation-design simulation CSV files and generate MCC curves.
# Input:   results/sim_results/design_comparison.csv
#          results/sim_results/design_nl_comparison.csv
#          results/sim_results/design_nl_seq_comparison.csv
# Output:  results/figure/design_mcc_curve.pdf
#          results/figure/design_nl_mcc_curve.pdf
#          results/figure/design_nl_seq_mcc_curve.pdf
################################################################################

rm(list = ls())

suppressMessages({
  library(ggplot2)
})

dir.create("results/figure", showWarnings = FALSE, recursive = TRUE)

plot_design_curve <- function(input_csv, output_pdf, title, subtitle,
                              dopt_label = "Sequential D-optimal") {
  if (!file.exists(input_csv)) {
    stop("Missing input file: ", input_csv,
         "\nRun the corresponding sim_script/02_scaling_design/pss_net_design*.R first.")
  }

  df <- read.csv(input_csv, stringsAsFactors = FALSE)
  required <- c("strategy", "N", "MCC")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("Input file ", input_csv, " is missing columns: ",
         paste(missing, collapse = ", "))
  }

  sm <- aggregate(MCC ~ strategy + N, df, function(z) {
    c(m = mean(z), se = sd(z) / sqrt(length(z)))
  })
  plot_df <- data.frame(
    strategy = sm$strategy,
    N = sm$N,
    MCC = sm$MCC[, "m"],
    se = sm$MCC[, "se"]
  )

  lab <- c(
    random = "Random",
    maximin = "Maximin (space-filling)",
    dopt = dopt_label
  )
  plot_df$strategy <- factor(lab[plot_df$strategy], levels = lab)

  p_curve <- ggplot(plot_df, aes(N, MCC, color = strategy, fill = strategy)) +
    geom_ribbon(aes(ymin = MCC - se, ymax = MCC + se),
                alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    scale_color_manual(values = c("grey50", "steelblue3", "tomato3")) +
    scale_fill_manual(values = c("grey50", "steelblue3", "tomato3")) +
    labs(title = title, subtitle = subtitle,
         x = "Number of perturbation conditions N", y = "MCC",
         color = NULL, fill = NULL) +
    theme_bw() +
    theme(legend.position = "bottom")

  ggsave(output_pdf, p_curve, width = 7, height = 5)
  cat("Saved:", output_pdf, "\n")
}

plot_design_curve(
  input_csv = "results/sim_results/design_comparison.csv",
  output_pdf = "results/figure/design_mcc_curve.pdf",
  title = "Optimal perturbation design improves network recovery",
  subtitle = "8-node linear gLV benchmark, monomial M=2, ADSIHT"
)

plot_design_curve(
  input_csv = "results/sim_results/design_nl_comparison.csv",
  output_pdf = "results/figure/design_nl_mcc_curve.pdf",
  title = "Optimal perturbation design under nonlinear interactions",
  subtitle = "8-node nonlinear benchmark with quadratic terms, ADSIHT",
  dopt_label = "Sequential D-optimal (active)"
)

plot_design_curve(
  input_csv = "results/sim_results/design_nl_seq_comparison.csv",
  output_pdf = "results/figure/design_nl_seq_mcc_curve.pdf",
  title = "Optimal perturbation design under strong nonlinear interactions",
  subtitle = "8-node nonlinear benchmark with strong quadratic terms, ADSIHT"
)

################################################################################
# summarize_mcc_comparison.R
#
# Purpose: Read per-run method-comparison metrics and write formatted tables.
# Input:   results/sim_results/mcc_comparison.csv
# Output:  results/table/mcc_comparison.txt
################################################################################

rm(list = ls())

input_csv <- "results/sim_results/mcc_comparison.csv"
output_txt <- "results/table/mcc_comparison.txt"

if (!file.exists(input_csv)) {
  stop("Missing input file: ", input_csv,
       "\nRun sim_script/pss_net_compare.R first.")
}

df_all <- read.csv(input_csv, stringsAsFactors = FALSE)
required <- c("version", "method", "seed", "MCC", "F1", "Pr", "Re",
              "TP", "FP", "FN", "CoefL2", "JacRMSE")
missing <- setdiff(required, names(df_all))
if (length(missing) > 0) {
  stop("Input file ", input_csv, " is missing columns: ",
       paste(missing, collapse = ", "))
}

summary_tab <- function(df, ver, meth) {
  sub <- df[df$version == ver & df$method == meth, ]
  c(version = ver, method = meth,
    MCC_mean = round(mean(sub$MCC), 3),
    MCC_sd = round(sd(sub$MCC), 3),
    MCC_min = round(min(sub$MCC), 3),
    MCC_max = round(max(sub$MCC), 3),
    F1_mean = round(mean(sub$F1), 3),
    Pr_mean = round(mean(sub$Pr), 3),
    Re_mean = round(mean(sub$Re), 3),
    TP_mean = round(mean(sub$TP), 2),
    FP_mean = round(mean(sub$FP), 2),
    FN_mean = round(mean(sub$FN), 2),
    CoefL2_mean = round(mean(sub$CoefL2), 4),
    JacRMSE_mean = round(mean(sub$JacRMSE), 4))
}

combos <- list(
  c("v3_no_smooth", "ADSIHT"),
  c("v3_no_smooth", "grLasso"),
  c("v1_pre_smooth", "ADSIHT"),
  c("v1_pre_smooth", "grLasso")
)

sum_df <- do.call(rbind, lapply(combos, function(x) {
  as.data.frame(t(summary_tab(df_all, x[1], x[2])),
                stringsAsFactors = FALSE)
}))

seeds <- sort(unique(df_all$seed))
n_runs <- length(seeds)
n_edges <- 18
M <- 2
n_nodes <- 8

write_summary <- function(con = stdout()) {
  cat("================================================================\n", file = con)
  cat(sprintf("  PSS-Net: v3 (no smooth) vs v1 (pre-smooth)  |  %d runs\n",
              n_runs), file = con)
  cat(sprintf("  %d nodes | %d true edges | M=%d | N_cond=300 | sigma=0.03\n",
              n_nodes, n_edges, M), file = con)
  cat("  Pre-smoothing: lm(x_j ~ bs(u_j, df=6)) per node\n", file = con)
  cat("================================================================\n", file = con)
  cat(sprintf("  %-18s %-8s  %6s+/%-4s  [%5s,%5s]  F1=%5s  Pr=%5s  Re=%5s  TP=%4s FP=%4s FN=%4s\n",
              "Version", "Method", "MCC", "sd", "min", "max",
              "mean", "mean", "mean", "avg", "avg", "avg"), file = con)
  cat("  ", strrep("-", 100), "\n", sep = "", file = con)

  for (i in seq_len(nrow(sum_df))) {
    r <- sum_df[i, ]
    if (i == 3) cat("  ", strrep(".", 100), "\n", sep = "", file = con)
    cat(sprintf("  %-18s %-8s  %6s+/%-4s  [%5s,%5s]  F1=%5s  Pr=%5s  Re=%5s  TP=%4s FP=%4s FN=%4s\n",
                r$version, r$method,
                r$MCC_mean, r$MCC_sd,
                r$MCC_min, r$MCC_max,
                r$F1_mean, r$Pr_mean, r$Re_mean,
                r$TP_mean, r$FP_mean, r$FN_mean), file = con)
  }
  cat("================================================================\n\n", file = con)

  cat("-- Per-run MCC -------------------------------------------------\n", file = con)
  cat(sprintf("  %5s  %10s  %10s  %10s  %10s  | d_ADS  d_GL\n",
              "Seed", "v3-ADSIHT", "v1-ADSIHT",
              "v3-grLasso", "v1-grLasso"), file = con)
  cat("  ", strrep("-", 72), "\n", sep = "", file = con)

  get_mcc <- function(seed, version, method) {
    df_all$MCC[df_all$seed == seed &
                 df_all$version == version &
                 df_all$method == method][1]
  }

  for (seed in seeds) {
    mcc_v3_a <- get_mcc(seed, "v3_no_smooth", "ADSIHT")
    mcc_v1_a <- get_mcc(seed, "v1_pre_smooth", "ADSIHT")
    mcc_v3_g <- get_mcc(seed, "v3_no_smooth", "grLasso")
    mcc_v1_g <- get_mcc(seed, "v1_pre_smooth", "grLasso")
    cat(sprintf("  %5d  %10.3f  %10.3f  %10.3f  %10.3f  | %+6.3f %+6.3f\n",
                seed, mcc_v3_a, mcc_v1_a, mcc_v3_g, mcc_v1_g,
                mcc_v1_a - mcc_v3_a, mcc_v1_g - mcc_v3_g), file = con)
  }
  cat("  ", strrep("-", 72), "\n", sep = "", file = con)
}

dir.create(dirname(output_txt), showWarnings = FALSE, recursive = TRUE)
write_summary(stdout())
out_con <- file(output_txt, open = "wt")
on.exit(close(out_con), add = TRUE)
write_summary(out_con)
cat("Saved:", output_txt, "\n")

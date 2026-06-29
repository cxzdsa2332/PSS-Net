rm(list = ls())

################################################################################
# Fig3.R -- Figure 3 panel objects for PSS-Net
#
# Purpose: build Figure 3 robustness/benchmark panel objects.
#            Fig3a -- method capability matrix (positioning panel; no data file,
#                     hand-curated from ref/external_benchmark_methods.md).
#            Fig3b -- external-method edge-recovery benchmark at a representative
#                     sample budget (default-selection MCC and ranked-edge AUPRC).
#                     This is the single benchmark panel; the former integrated
#                     MCC/AUPRC/FuncNRMSE "Fig3c" panel was a near-duplicate of
#                     Fig3b and has been removed.
#            Fig3c -- topology schematic: scale-free vs ER example networks
#                     drawn as a top/bottom comparison (node size = out-degree).
#            Fig3d -- node-wise vs joint PSS-Net structure dependence on those
#                     topologies (paired scatter vs the identity line).
#            Fig3e -- compositional-data limitations for PSS network recovery
#                     (simulation file retains its historical Fig3c name).
#            (diagnostic, not a main panel) Fig3c_benchmark_recovery_curves --
#                     edge-strength (Jacobian NRMSE) and nonlinear edge-function
#                     (function NRMSE) recovery vs sample budget, regression
#                     family only.
#
# Input:   results/sim_results/Fig3b_external_benchmark_main.csv
#          results/sim_results/Fig3c_structure_dependence.csv
#          results/sim_results/Fig3c_compositional_data_limitation.csv
# Output:  Fig3a_method_capability_matrix / Fig3a,
#          Fig3b_benchmark_edge_recovery / Fig3b,
#          Fig3c_topology_schematic / Fig3c,
#          Fig3d_structure_results / Fig3d,
#          Fig3e_compositional_data_limitation / Fig3e objects in workspace.
################################################################################

suppressMessages({
  library(ggplot2)
  library(patchwork)
})

## ------------------------------------- Fig3a: method capability matrix ----
# Qualitative positioning panel: which capabilities each method class has.
# "Uses perturbation u" and "steady-state (no time series)" were merged because
# every method in the latter column was marked yes, so that column added no
# discrimination. The combined column specifically means that both X* and the
# matched perturbation input u are used in a steady-state estimating equation.
# Hand-curated from ref/external_benchmark_methods.md.
# Encoding: 1 = yes, 0.5 = partial / conditional, 0 = no.
fig3a_methods <- c(
  "PSS-Net (ADSIHT)", "MRA", "Lasso / Elastic net", "Group lasso",
  "PySINDy STLSQ (PSS)", "GENIE3 / GRNBoost2", "GIES (interventional)",
  "Correlation / partial cor.", "SPIEC-EASI / SparCC"
)
fig3a_caps <- c("Perturbed steady-state\ninput (X*, u)", "Directed\nedges", "Signed\nedges",
                "Nonlinear\nedge functions", "Built-in sparse\nselection",
                "Supports scale-free /\nhub topology", "Compositional\naware")
# Rows follow fig3a_methods; columns follow fig3a_caps. "Built-in sparse
# selection" = produces a sparse edge set directly, without a post-hoc threshold
# (so the dense MRA / correlation baselines are 0; they need thresholding).
# "Supports scale-free / hub topology" means that the method does not impose a
# homogeneous-degree topology. A partial mark indicates an important restriction:
# GIES only learns DAGs, while association methods do not identify directed hubs.
fig3a_values <- rbind(
  c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0),  # PSS-Net (group + within-group sparsity)
  c(1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0),  # MRA: signed normalized local responses
  c(1.0, 1.0, 1.0, 0.5, 1.0, 1.0, 0.0),  # Lasso / EN (nonlinear via library)
  c(1.0, 1.0, 1.0, 0.5, 1.0, 1.0, 0.0),  # Group lasso
  c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0),  # PySINDy STLSQ on steady-state PSS library
  c(0.0, 1.0, 0.0, 1.0, 0.5, 1.0, 0.5),  # GENIE3 (state-only, ranked importances)
  c(1.0, 1.0, 0.0, 0.0, 1.0, 0.5, 0.0),  # GIES: scale-free DAGs only; no feedback cycles
  c(0.0, 0.0, 1.0, 0.0, 0.0, 0.5, 0.0),  # Correlation / partial cor.: undirected hubs only
  c(0.0, 0.0, 0.5, 0.0, 0.5, 0.5, 1.0)   # SPIEC-EASI / SparCC: undirected hubs only
)

fig3a_df <- data.frame(
  method = factor(rep(fig3a_methods, times = length(fig3a_caps)),
                  levels = rev(fig3a_methods)),
  capability = factor(rep(fig3a_caps, each = length(fig3a_methods)),
                      levels = fig3a_caps),
  value = as.vector(fig3a_values)
)
fig3a_df$glyph <- c("No", "Partial", "Yes")[match(fig3a_df$value, c(0, 0.5, 1))]
fig3a_df$val_f <- factor(fig3a_df$value, levels = c(0, 0.5, 1))

# Highlight the PSS-Net row (top of the reversed factor).
fig3a_hi <- length(fig3a_methods)

Fig3a_method_capability_matrix <- ggplot(
  fig3a_df, aes(x = capability, y = method, fill = val_f)
) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = glyph), size = 2.35, color = "grey15") +
  annotate("rect", xmin = 0.5, xmax = length(fig3a_caps) + 0.5,
           ymin = fig3a_hi - 0.5, ymax = fig3a_hi + 0.5,
           fill = NA, color = "#2E6F9E", linewidth = 0.9) +
  scale_fill_manual(values = c("0" = "#EEEEEE", "0.5" = "#F3D9B0",
                               "1" = "#BCD3B6"), guide = "none") +
  scale_x_discrete(position = "top") +
  labs(
    x = NULL, y = NULL,
    title = "Capability profile of network-inference methods"
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

## -------------------------------- Fig3b: external edge-recovery benchmark ----
fig3b_file <- "results/sim_results/Fig3b_external_benchmark_main.csv"
if (!file.exists(fig3b_file)) {
  stop("Missing ", fig3b_file,
       ". Run sim_script/03_robustness_benchmarks/",
       "Fig3b_external_benchmark_main.R first.")
}

fig3b_raw <- read.csv(fig3b_file, stringsAsFactors = FALSE)
if (!"M_ord" %in% names(fig3b_raw) || anyNA(fig3b_raw$M_ord) ||
    any(fig3b_raw$M_ord != 2L)) {
  warning(
    "Fig3b CSV predates the current M_ord = 2 benchmark schema; rerun ",
    "Fig3b_external_benchmark_main.R before using this panel as a current result."
  )
}
required_cols <- c(
  "seed", "N", "N_over_slogp", "truth", "method", "MCC", "AUPRC"
)
missing_cols <- setdiff(required_cols, names(fig3b_raw))
if (length(missing_cols) > 0L) {
  stop("Missing required columns in Fig3b_external_benchmark_main.csv: ",
       paste(missing_cols, collapse = ", "))
}

# Use the planned representative budget N / (s log(p)) = 25. Selecting the
# nearest available budget keeps the plotting script compatible with a rerun
# whose integer N differs slightly because p or s changed.
fig3b_target_budget <- 25
fig3b_budgets <- aggregate(N_over_slogp ~ N, fig3b_raw, median, na.rm = TRUE)
fig3b_representative_N <- fig3b_budgets$N[
  which.min(abs(fig3b_budgets$N_over_slogp - fig3b_target_budget))
]
fig3b_representative_ratio <- fig3b_budgets$N_over_slogp[
  fig3b_budgets$N == fig3b_representative_N
][1]
fig3b_plot_raw <- fig3b_raw[fig3b_raw$N == fig3b_representative_N, ]

if (nrow(fig3b_plot_raw) == 0L) {
  stop("No rows available at the representative Fig3b sample budget.")
}

fig3b_method_order <- c(
  "PSS_Net", "aiMeRA", "Lasso", "ElasticNet", "GroupLasso",
  "PySINDy_STLSQ", "Correlation", "PartialCor"
)
fig3b_present_methods <- unique(fig3b_plot_raw$method)
fig3b_method_order <- c(
  intersect(fig3b_method_order, fig3b_present_methods),
  setdiff(fig3b_present_methods, fig3b_method_order)
)
# Reverse factor levels so the first (PSS-Net) method is displayed at the top.
fig3b_plot_raw$method <- factor(
  fig3b_plot_raw$method, levels = rev(fig3b_method_order)
)

fig3b_truth_order <- c("linear", "strong_nonlinear")
fig3b_truth_order <- c(
  intersect(fig3b_truth_order, unique(fig3b_plot_raw$truth)),
  setdiff(unique(fig3b_plot_raw$truth), fig3b_truth_order)
)
fig3b_plot_raw$truth <- factor(fig3b_plot_raw$truth,
                               levels = fig3b_truth_order)

fig3b_key_cols <- c("seed", "N", "N_over_slogp", "truth", "method")
fig3b_metric_df <- rbind(
  data.frame(fig3b_plot_raw[, fig3b_key_cols],
             metric = "MCC", value = fig3b_plot_raw$MCC),
  data.frame(fig3b_plot_raw[, fig3b_key_cols],
             metric = "AUPRC", value = fig3b_plot_raw$AUPRC)
)
fig3b_metric_df <- fig3b_metric_df[is.finite(fig3b_metric_df$value), ]
fig3b_metric_df$metric <- factor(fig3b_metric_df$metric,
                                 levels = c("MCC", "AUPRC"))

fig3b_summary <- do.call(rbind, lapply(split(
  fig3b_metric_df,
  list(fig3b_metric_df$truth, fig3b_metric_df$method,
       fig3b_metric_df$metric),
  drop = TRUE
), function(d) {
  data.frame(
    truth = d$truth[1],
    method = d$method[1],
    metric = d$metric[1],
    mean = mean(d$value),
    sd = sd(d$value),
    n_seed = length(unique(d$seed)),
    stringsAsFactors = FALSE
  )
}))
rownames(fig3b_summary) <- NULL
fig3b_summary$truth <- factor(fig3b_summary$truth,
                              levels = fig3b_truth_order)
fig3b_summary$method <- factor(fig3b_summary$method,
                               levels = rev(fig3b_method_order))
fig3b_summary$metric <- factor(fig3b_summary$metric,
                               levels = c("MCC", "AUPRC"))
fig3b_summary$lower <- pmax(0, fig3b_summary$mean - fig3b_summary$sd)
fig3b_summary$upper <- pmin(1, fig3b_summary$mean + fig3b_summary$sd)

fig3b_method_labels <- setNames(fig3b_method_order, fig3b_method_order)
fig3b_label_updates <- c(
  PSS_Net = "PSS-Net", aiMeRA = "aiMeRA (MRA)",
  ElasticNet = "Elastic net", GroupLasso = "Group lasso",
  PySINDy_STLSQ = "PySINDy STLSQ", PartialCor = "Partial correlation"
)
label_keys <- intersect(names(fig3b_label_updates),
                        names(fig3b_method_labels))
fig3b_method_labels[label_keys] <- fig3b_label_updates[label_keys]

fig3b_method_colors <- c(
  PSS_Net = "#2E6F9E", aiMeRA = "#B45F4D", Lasso = "#D18B47",
  ElasticNet = "#C4A33B", GroupLasso = "#6D8B3D",
  PySINDy_STLSQ = "#8B6BAE", Correlation = "#737373",
  PartialCor = "#A0A0A0"
)
extra_methods <- setdiff(fig3b_method_order, names(fig3b_method_colors))
if (length(extra_methods) > 0L) {
  fig3b_method_colors[extra_methods] <- "#555555"
}

fig3b_subtitle <- sprintf(
  "M = 2; N = %d (N / [s log(p)] = %.1f). PSS-Net/MRA: |row-normalized link| > 0.05; points: seeds; bars: mean ± SD",
  fig3b_representative_N, fig3b_representative_ratio
)
Fig3b_benchmark_edge_recovery <- ggplot() +
  geom_point(
    data = fig3b_metric_df,
    aes(x = value, y = method),
    position = position_jitter(width = 0, height = 0.08),
    size = 1.2, alpha = 0.32, color = "grey25"
  ) +
  geom_segment(
    data = fig3b_summary,
    aes(x = lower, xend = upper, y = method, yend = method,
        color = method),
    linewidth = 0.65, show.legend = FALSE
  ) +
  geom_point(
    data = fig3b_summary,
    aes(x = mean, y = method, color = method),
    size = 2.35, show.legend = FALSE
  ) +
  facet_grid(
    metric ~ truth,
    labeller = labeller(
      truth = c(linear = "Linear truth",
                strong_nonlinear = "Nonlinear truth"),
      metric = c(MCC = "MCC\n(default selection)",
                 AUPRC = "AUPRC\n(edge ranking)")
    )
  ) +
  scale_x_continuous(
    limits = c(-0.02, 1.02), breaks = seq(0, 1, by = 0.25),
    expand = expansion(mult = 0)
  ) +
  scale_y_discrete(labels = fig3b_method_labels) +
  scale_color_manual(values = fig3b_method_colors) +
  labs(
    x = "edge-recovery score", y = NULL,
    title = "External benchmark separates edge selection from edge ranking",
    subtitle = fig3b_subtitle
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 10.8),
    plot.subtitle = element_text(size = 7.8, color = "grey30"),
    strip.background = element_rect(fill = "grey95", color = "grey82",
                                    linewidth = 0.35),
    strip.text = element_text(face = "bold", size = 8, lineheight = 0.95),
    axis.title = element_text(size = 8.5),
    axis.text.x = element_text(size = 7.2, color = "grey25"),
    axis.text.y = element_text(size = 7.2, color = "grey20"),
    panel.spacing = unit(0.8, "lines"),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

Fig3b <- Fig3b_benchmark_edge_recovery

Fig3b

## ---- Diagnostic (supplement): edge-strength / function recovery curves ----
# Not a main Fig3 panel (object kept as Fig3c_benchmark_recovery_curves for the
# supplement / simulation_plan reference); the Fig3c panel is now the topology
# schematic below.
fig3c_required_cols <- c(
  "seed", "N", "N_over_slogp", "truth", "method", "JacScale",
  "EdgeJacNRMSE", "FuncNRMSE"
)
fig3c_missing_cols <- setdiff(fig3c_required_cols, names(fig3b_raw))
if (length(fig3c_missing_cols) > 0L) {
  stop(
    "Fig3b CSV lacks the edge-strength columns required by Fig3c: ",
    paste(fig3c_missing_cols, collapse = ", "),
    ". Rerun Fig3b_external_benchmark_main.R."
  )
}

# Only equation-regression methods estimate absolute-scale edge functions.
# aiMeRA is intentionally excluded: its Jacobian is row-normalized and it does
# not return a comparable f_ji(x). Association methods return neither quantity.
fig3c_method_order <- c(
  "PSS_Net", "Lasso", "ElasticNet", "GroupLasso", "PySINDy_STLSQ"
)
fig3c_method_order <- intersect(fig3c_method_order, unique(fig3b_raw$method))
fig3c_recovery_raw <- fig3b_raw[
  fig3b_raw$method %in% fig3c_method_order &
    fig3b_raw$JacScale == "absolute",
]

fig3c_key_cols <- c("seed", "N", "N_over_slogp", "truth", "method")
fig3c_recovery_df <- rbind(
  data.frame(
    fig3c_recovery_raw[, fig3c_key_cols],
    metric = "EdgeJacNRMSE", value = fig3c_recovery_raw$EdgeJacNRMSE
  ),
  data.frame(
    fig3c_recovery_raw[, fig3c_key_cols],
    metric = "FuncNRMSE", value = fig3c_recovery_raw$FuncNRMSE
  )
)
fig3c_recovery_df <- fig3c_recovery_df[
  is.finite(fig3c_recovery_df$value),
]
fig3c_recovery_df$truth <- factor(
  fig3c_recovery_df$truth, levels = fig3b_truth_order
)
fig3c_recovery_df$method <- factor(
  fig3c_recovery_df$method, levels = fig3c_method_order
)
fig3c_recovery_df$metric <- factor(
  fig3c_recovery_df$metric,
  levels = c("EdgeJacNRMSE", "FuncNRMSE")
)

fig3c_recovery_summary <- do.call(rbind, lapply(split(
  fig3c_recovery_df,
  list(
    fig3c_recovery_df$truth, fig3c_recovery_df$method,
    fig3c_recovery_df$N, fig3c_recovery_df$metric
  ),
  drop = TRUE
), function(d) {
  data.frame(
    truth = d$truth[1],
    method = d$method[1],
    N = d$N[1],
    N_over_slogp = median(d$N_over_slogp),
    metric = d$metric[1],
    mean = mean(d$value),
    sd = sd(d$value),
    n_seed = length(unique(d$seed)),
    stringsAsFactors = FALSE
  )
}))
rownames(fig3c_recovery_summary) <- NULL
fig3c_recovery_summary$truth <- factor(
  fig3c_recovery_summary$truth, levels = fig3b_truth_order
)
fig3c_recovery_summary$method <- factor(
  fig3c_recovery_summary$method, levels = fig3c_method_order
)
fig3c_recovery_summary$metric <- factor(
  fig3c_recovery_summary$metric,
  levels = c("EdgeJacNRMSE", "FuncNRMSE")
)
fig3c_recovery_summary$lower <- pmax(
  0, fig3c_recovery_summary$mean - fig3c_recovery_summary$sd
)
fig3c_recovery_summary$upper <-
  fig3c_recovery_summary$mean + fig3c_recovery_summary$sd

fig3c_method_labels <- fig3b_method_labels[fig3c_method_order]
fig3c_method_colors <- fig3b_method_colors[fig3c_method_order]
fig3c_x_breaks <- sort(unique(fig3c_recovery_summary$N_over_slogp))

Fig3c_benchmark_function_recovery <- ggplot(
  fig3c_recovery_summary,
  aes(
    x = N_over_slogp, y = mean, color = method, fill = method,
    group = method
  )
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.07, color = NA, show.legend = FALSE
  ) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  facet_grid(
    metric ~ truth, scales = "free_y",
    labeller = labeller(
      truth = c(linear = "Linear truth",
                strong_nonlinear = "Nonlinear truth"),
      metric = c(
        EdgeJacNRMSE = "Local edge strength\nJacobian NRMSE",
        FuncNRMSE = "Full edge function\nfunction NRMSE"
      )
    )
  ) +
  scale_x_continuous(
    breaks = fig3c_x_breaks,
    labels = function(x) format(round(x, 1), nsmall = 1),
    expand = expansion(mult = c(0.03, 0.05))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.08))) +
  scale_color_manual(
    values = fig3c_method_colors, labels = fig3c_method_labels
  ) +
  scale_fill_manual(
    values = fig3c_method_colors, labels = fig3c_method_labels
  ) +
  labs(
    x = expression("sample budget " * N / (s * log(p))),
    y = "relative recovery error (lower is better)",
    color = NULL, fill = NULL,
    title = "PSS-based methods recover local edge strength and nonlinear edge functions",
    subtitle = paste0(
      "Polynomial library order M = 2; curves are means and ribbons are ± SD. ",
      "MRA is excluded because its response matrix is row-normalized and has no edge function."
    )
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 10.8),
    plot.subtitle = element_text(size = 7.8, color = "grey30"),
    strip.background = element_rect(
      fill = "grey95", color = "grey82", linewidth = 0.35
    ),
    strip.text = element_text(face = "bold", size = 8, lineheight = 0.95),
    axis.title = element_text(size = 8.5),
    axis.text = element_text(size = 7.2, color = "grey25"),
    legend.position = "bottom",
    legend.text = element_text(size = 7.5),
    legend.key.width = unit(1.1, "lines"),
    panel.spacing = unit(0.8, "lines"),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

Fig3c_benchmark_recovery_curves <- Fig3c_benchmark_function_recovery

## ------------------------- Fig3c: topology schematic (scale-free vs ER) ----
# Fig3c is a standalone schematic of the two contrasting topologies, drawn as a
# TOP/BOTTOM comparison (scale-free above, ER below): small example graphs on a
# circle with node size / colour growing with out-degree, so the scale-free hubs
# are visibly large and the ER nodes look uniform. The matched node-wise vs joint
# recovery on these topologies is the separate Fig3d panel that follows.

fig3d_topo_labels <- c(scalefree = "Scale-free (hub)",
                       er = "Erdos-Renyi (homogeneous)")

# -- small illustrative networks (same generative rule as the simulation) --
# pa_power sharpens the preferential attachment for the SCHEMATIC only: the
# source is chosen with probability proportional to out-degree^pa_power, so a
# couple of dominant hubs emerge and the contrast with the uniform ER graph is
# unmistakable. The actual simulation (Fig3c_structure_dependence.R) uses the
# milder pa_power = 2.
fig3d_demo_net <- function(topology, p = 18L, avg_in = 2L, seed = 4L,
                           pa_power = 3) {
  set.seed(seed)
  A <- matrix(0L, p, p)
  if (topology == "scalefree") {
    outdeg <- rep(1, p)                       # +1 smoothing = preferential weight
    for (j in seq_len(p)) {
      k_j <- min(1L + rpois(1, avg_in - 1), p - 1L)
      chosen <- integer(0)
      for (t in seq_len(k_j)) {
        w <- outdeg^pa_power; w[c(j, chosen)] <- 0
        i <- sample.int(p, 1, prob = w)       # P(source) prop. to out-degree^power
        chosen <- c(chosen, i)
        A[j, i] <- 1L
        outdeg[i] <- outdeg[i] + 1
      }
    }
  } else {
    for (j in seq_len(p)) {
      src <- sample(setdiff(seq_len(p), j), avg_in)
      A[j, src] <- 1L
    }
  }
  A
}

fig3d_circle <- function(p) {
  ang <- seq(0, 2 * pi, length.out = p + 1L)[seq_len(p)]
  data.frame(node = seq_len(p), x = cos(ang), y = sin(ang))
}

fig3d_nodes <- list()
fig3d_edges <- list()
for (tp in c("scalefree", "er")) {
  A <- fig3d_demo_net(tp)
  p_demo <- nrow(A)
  lay <- fig3d_circle(p_demo)
  nd <- lay
  nd$outdeg <- colSums(A)                     # source out-degree (hub size)
  nd$topology <- tp
  fig3d_nodes[[tp]] <- nd
  idx <- which(A == 1L, arr.ind = TRUE)       # A[j, i] = 1 means edge i -> j
  if (nrow(idx) > 0L) {
    fig3d_edges[[tp]] <- data.frame(
      x = lay$x[idx[, "col"]], y = lay$y[idx[, "col"]],     # source
      xend = lay$x[idx[, "row"]], yend = lay$y[idx[, "row"]], # target
      topology = tp
    )
  }
}
fig3d_nodes_df <- do.call(rbind, fig3d_nodes)
fig3d_edges_df <- do.call(rbind, fig3d_edges)
fig3d_nodes_df$topology <- factor(fig3d_nodes_df$topology,
                                  levels = c("scalefree", "er"))
fig3d_edges_df$topology <- factor(fig3d_edges_df$topology,
                                  levels = c("scalefree", "er"))

Fig3c_topology_schematic <- ggplot() +
  geom_curve(
    data = fig3d_edges_df,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "grey70", linewidth = 0.3, alpha = 0.7, curvature = 0.2,
    arrow = arrow(length = unit(0.05, "inches"), type = "closed")
  ) +
  geom_point(
    data = fig3d_nodes_df,
    aes(x = x, y = y, size = outdeg, fill = outdeg),
    shape = 21, color = "grey35", stroke = 0.3
  ) +
  facet_wrap(~ topology, ncol = 1,
             labeller = labeller(topology = fig3d_topo_labels)) +
  scale_size_continuous(range = c(1.5, 11), guide = "none") +
  scale_fill_gradient(low = "#DCE6EF", high = "#1F4E79", guide = "none") +
  coord_equal(clip = "off") +
  labs(title = "Two contrasting network topologies") +
  theme_void(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 10.8, hjust = 0,
                              margin = margin(b = 8)),
    strip.text = element_text(face = "bold", size = 8.5, color = "grey20",
                              margin = margin(b = 2, t = 4)),
    panel.spacing = unit(0.6, "lines"),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

Fig3c <- Fig3c_topology_schematic

Fig3c

## ------------------------- Fig3d: node-wise vs joint structure dependence ----
fig3d_file <- "results/sim_results/Fig3c_structure_dependence.csv"
if (!file.exists(fig3d_file)) {
  stop("Missing ", fig3d_file,
       ". Run sim_script/03_robustness_benchmarks/",
       "Fig3c_structure_dependence.R first.")
}
fig3d_raw <- read.csv(fig3d_file, stringsAsFactors = FALSE)
# Three metric rows: overall edge recovery (MCC) plus TWO hub-mining metrics
# (out-degree Spearman rho and top-k hub hit-rate), so the panel foregrounds how
# well each scheme recovers the hubs, not just the edge set.
fig3d_metric_levels <- c("Edge MCC", "Hub out-degree rho", "Top-k hub hit-rate")
fig3d_keys <- c("topology", "N", "seed")
fig3d_long <- rbind(
  data.frame(fig3d_raw[, fig3d_keys], scheme = "node-wise",
             metric = "Edge MCC", value = fig3d_raw$MCC_nodewise),
  data.frame(fig3d_raw[, fig3d_keys], scheme = "joint",
             metric = "Edge MCC", value = fig3d_raw$MCC_joint),
  data.frame(fig3d_raw[, fig3d_keys], scheme = "node-wise",
             metric = "Hub out-degree rho", value = fig3d_raw$rho_nodewise),
  data.frame(fig3d_raw[, fig3d_keys], scheme = "joint",
             metric = "Hub out-degree rho", value = fig3d_raw$rho_joint),
  data.frame(fig3d_raw[, fig3d_keys], scheme = "node-wise",
             metric = "Top-k hub hit-rate", value = fig3d_raw$topk_nodewise),
  data.frame(fig3d_raw[, fig3d_keys], scheme = "joint",
             metric = "Top-k hub hit-rate", value = fig3d_raw$topk_joint)
)
fig3d_long <- fig3d_long[is.finite(fig3d_long$value), ]

fig3d_summary <- do.call(rbind, lapply(split(
  fig3d_long,
  list(fig3d_long$topology, fig3d_long$N, fig3d_long$scheme,
       fig3d_long$metric),
  drop = TRUE
), function(d) {
  data.frame(
    topology = d$topology[1], N = d$N[1], scheme = d$scheme[1],
    metric = d$metric[1], mean = mean(d$value), sd = sd(d$value),
    stringsAsFactors = FALSE
  )
}))
rownames(fig3d_summary) <- NULL
fig3d_summary$topology <- factor(fig3d_summary$topology,
                                 levels = c("scalefree", "er"))
fig3d_summary$scheme <- factor(fig3d_summary$scheme,
                               levels = c("node-wise", "joint"))
fig3d_summary$metric <- factor(fig3d_summary$metric,
                               levels = fig3d_metric_levels)
fig3d_summary$lower <- pmax(-1, fig3d_summary$mean - fig3d_summary$sd)
fig3d_summary$upper <- pmin(1, fig3d_summary$mean + fig3d_summary$sd)

fig3d_scheme_colors <- c("node-wise" = "#B45F4D", "joint" = "#2E6F9E")

# Paired scatter against the identity line (not a line / dumbbell / bar chart):
# each point is one (topology, N) cell, x = node-wise recovery, y = joint
# recovery, so a point ABOVE the dashed y = x line means the joint solve wins.
# Colour = topology, point size = sample budget N, facet = metric. Scale-free
# points (especially on the hub metrics) sit clearly above the diagonal while
# the ER points hug it, making the structure-dependent advantage obvious in a
# chart type used nowhere else in Fig3.
fig3d_nw <- fig3d_summary[fig3d_summary$scheme == "node-wise",
                          c("topology", "N", "metric", "mean")]
fig3d_jt <- fig3d_summary[fig3d_summary$scheme == "joint",
                          c("topology", "N", "metric", "mean")]
names(fig3d_nw)[4] <- "nodewise"
names(fig3d_jt)[4] <- "joint"
fig3d_db <- merge(fig3d_nw, fig3d_jt, by = c("topology", "N", "metric"))
fig3d_db$topology <- factor(fig3d_db$topology, levels = c("scalefree", "er"))
fig3d_db$metric <- factor(fig3d_db$metric, levels = fig3d_metric_levels)

fig3d_topo_colors <- c(scalefree = "#1F4E79", er = "#C0883B")
fig3d_lim <- c(0, max(fig3d_db$nodewise, fig3d_db$joint) * 1.06)

Fig3d_structure_results <- ggplot(
  fig3d_db, aes(x = nodewise, y = joint, color = topology)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "22",
              color = "grey55", linewidth = 0.4) +
  geom_point(aes(size = N), alpha = 0.85) +
  facet_wrap(~ metric, nrow = 1) +
  scale_color_manual(values = fig3d_topo_colors, labels = fig3d_topo_labels) +
  scale_size_continuous(range = c(1.6, 4.6),
                        breaks = sort(unique(fig3d_db$N))) +
  coord_equal(xlim = fig3d_lim, ylim = fig3d_lim, clip = "off") +
  labs(
    x = "node-wise recovery", y = "joint recovery",
    color = NULL, size = "N",
    title = "Joint vs node-wise PSS-Net (above dashed line = joint better)"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 10.4),
    strip.background = element_rect(fill = "grey95", color = "grey82",
                                    linewidth = 0.35),
    strip.text = element_text(face = "bold", size = 8, lineheight = 0.95),
    axis.title = element_text(size = 8.5),
    axis.text = element_text(size = 7.2, color = "grey25"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 7.5),
    panel.spacing = unit(0.8, "lines"),
    plot.margin = margin(2, 8, 5.5, 5.5)
  )

Fig3d <- Fig3d_structure_results

Fig3d

## ----------------------------------------- Fig3e: compositional limitation ----
fig3c_file <- "results/sim_results/Fig3c_compositional_data_limitation.csv"
if (!file.exists(fig3c_file)) {
  stop("Missing ", fig3c_file,
       ". Run sim_script/03_robustness_benchmarks/Fig3e_compositional_data_limitation.R first.")
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
  geom_errorbar(aes(ymin = pmax(-0.4, mean - sd), ymax = pmin(1, mean + sd)),
                width = 0.16, linewidth = 0.35, color = "grey20") +
  geom_point(data = fig3c_metric_df,
             aes(x = input, y = value),
             inherit.aes = FALSE, position = position_jitter(width = 0.08, height = 0),
             size = 1.25, alpha = 0.48, color = "grey15") +
  facet_wrap(~ metric, nrow = 1) +
  scale_x_discrete(labels = input_labels) +
  scale_y_continuous(limits = c(-0.4, 1), expand = expansion(mult = c(0.02, 0.05))) +
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

Fig3e_compositional_data_limitation <- Fig3c_compositional_data_limitation
Fig3e <- Fig3e_compositional_data_limitation

Fig3e

## ------------------------------------------------------- assembled Figure 3 ----
if (requireNamespace("patchwork", quietly = TRUE)) {
  Fig3 <- (patchwork::wrap_elements(Fig3a) | patchwork::wrap_elements(Fig3b)) /
    (patchwork::wrap_elements(Fig3c) | patchwork::wrap_elements(Fig3d)) /
    patchwork::wrap_elements(Fig3e) +
    patchwork::plot_layout(heights = c(0.8, 1.0, 0.75)) +
    patchwork::plot_annotation(tag_levels = "a")
} else {
  Fig3 <- NULL
}

Fig3

fig3_out <- file.path("manuscript", "figures", "Fig3.pdf")
dir.create(dirname(fig3_out), recursive = TRUE, showWarnings = FALSE)
if (!is.null(Fig3)) {
  ggsave(fig3_out, Fig3, width = 210, height = 297, units = "mm")
}

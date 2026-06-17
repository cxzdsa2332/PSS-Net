#!/usr/bin/env Rscript

# Purpose: draw a conceptual figure for PSS-Net.
# Inputs: none.
# Outputs:
#   results/figure/pss_concept_diagram.pdf
#   results/figure/pss_concept_diagram.png

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

out_dir <- file.path("results", "figure")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Palette: restrained, manuscript-friendly, and not tied to a single hue.
pal <- list(
  ink = "#24343D",
  muted = "#5B6770",
  light = "#F6F7F4",
  panel = "#FFFFFF",
  border = "#C9D0D3",
  green = "#3B8F6A",
  teal = "#2F7F8F",
  blue = "#4267AC",
  red = "#B84A4A",
  gold = "#D89A2B",
  grey = "#8D969B",
  pale_grey = "#E8ECEE"
)

rect_df <- data.frame(
  xmin = c(0.05, 0.28, 0.53, 0.77, 0.32),
  xmax = c(0.23, 0.48, 0.72, 0.96, 0.68),
  ymin = c(0.18, 0.18, 0.18, 0.18, 0.805),
  ymax = c(0.78, 0.78, 0.78, 0.78, 0.91),
  fill = c(pal$panel, pal$panel, pal$panel, pal$panel, "#F2F3F3"),
  color = c(pal$border, pal$border, pal$border, pal$border, "#BFC6CA")
)

panel_labels <- data.frame(
  x = c(0.14, 0.38, 0.625, 0.865, 0.50),
  y = c(0.73, 0.73, 0.73, 0.73, 0.875),
  label = c("Perturbed\nsteady states", "Interpretable\ndynamic model",
            "Sparse system\ninference", "Core nodes &\nregulatory effects",
            "Black-box alternative"),
  color = c(pal$teal, pal$blue, pal$green, pal$red, pal$grey)
)

# Mini data matrix: rows are perturbation conditions, columns are variables/OTUs.
grid_pts <- expand.grid(row = seq_len(7), col = seq_len(5))
grid_pts$x <- 0.085 + (grid_pts$col - 1) * 0.025
grid_pts$y <- 0.57 - (grid_pts$row - 1) * 0.045
grid_pts$val <- c(
  0.25, 0.55, 0.35, 0.76, 0.45,
  0.65, 0.30, 0.78, 0.42, 0.26,
  0.44, 0.72, 0.36, 0.33, 0.68,
  0.22, 0.40, 0.58, 0.81, 0.51,
  0.78, 0.48, 0.29, 0.60, 0.35,
  0.38, 0.63, 0.70, 0.27, 0.54,
  0.55, 0.34, 0.49, 0.67, 0.30
)

grid_palette <- colorRampPalette(c("#DCEEF0", "#F4E7B5", "#C95A4A"))(100)
grid_pts$fill_code <- grid_palette[pmax(1, pmin(100, round(grid_pts$val * 99) + 1))]

perturb_df <- data.frame(
  x = c(0.11, 0.145, 0.18, 0.11, 0.18),
  y = c(0.64, 0.645, 0.63, 0.235, 0.245),
  xend = c(0.11, 0.145, 0.18, 0.11, 0.18),
  yend = c(0.69, 0.70, 0.68, 0.19, 0.20),
  color = c(pal$green, pal$gold, pal$red, pal$blue, pal$green)
)

arrow_df <- data.frame(
  x = c(0.235, 0.49, 0.73),
  y = c(0.48, 0.48, 0.48),
  xend = c(0.275, 0.525, 0.765),
  yend = c(0.48, 0.48, 0.48)
)

black_arrow <- data.frame(
  x = c(0.40, 0.40, 0.60),
  y = c(0.78, 0.785, 0.785),
  xend = c(0.40, 0.60, 0.60),
  yend = c(0.805, 0.785, 0.805)
)

# Small inferred directed network.
node_df <- data.frame(
  name = c("A", "B", "C", "D", "E", "F"),
  x = c(0.84, 0.90, 0.865, 0.815, 0.92, 0.855),
  y = c(0.58, 0.58, 0.48, 0.39, 0.36, 0.30),
  size = c(7.2, 4.6, 5.5, 3.8, 3.4, 3.2),
  fill = c(pal$red, pal$blue, pal$green, pal$gold, pal$teal, pal$grey)
)
edge_df <- data.frame(
  x = c(0.84, 0.84, 0.90, 0.865, 0.865, 0.815),
  y = c(0.58, 0.58, 0.58, 0.48, 0.48, 0.39),
  xend = c(0.90, 0.865, 0.865, 0.815, 0.92, 0.855),
  yend = c(0.58, 0.48, 0.48, 0.39, 0.36, 0.30),
  effect = c("+", "-", "+", "-", "+", "-")
)
edge_df$color <- ifelse(edge_df$effect == "+", pal$red, pal$blue)
edge_df$linetype <- ifelse(edge_df$effect == "+", "solid", "dashed")

black_nodes <- data.frame(
  x = c(0.40, 0.48, 0.56),
  y = c(0.845, 0.84, 0.845),
  label = c("?", "?", "?")
)

main_text <- data.frame(
  x = c(0.14, 0.38, 0.625, 0.865, 0.495),
  y = c(0.135, 0.135, 0.135, 0.135, 0.795),
  label = c(
    "WT plus targeted perturbations\nare observed after relaxation",
    "Steady-state equations expose\nnode-wise additive regulation",
    "Group-sparse regression selects\nsource-to-target effects",
    "Hubs, signs, and edge functions\nbecome testable hypotheses",
    ""
  )
)

p <- ggplot() +
  geom_rect(data = rect_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = rect_df$fill, color = rect_df$color, linewidth = 0.45) +
  geom_text(data = panel_labels,
            aes(x = x, y = y, label = label, color = color),
            size = 3.25, fontface = "bold", lineheight = 0.9) +
  scale_color_identity() +
  geom_segment(data = arrow_df,
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.025, "npc"), type = "closed"),
               color = pal$muted, linewidth = 0.6) +
  geom_segment(data = black_arrow,
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.018, "npc"), type = "closed"),
               color = pal$grey, linewidth = 0.45, linetype = "dashed") +
  annotate("text", x = 0.50, y = 0.825, label = "prediction possible, mechanism hidden",
           size = 2.55, color = pal$grey) +
  geom_point(data = black_nodes, aes(x = x, y = y),
             shape = 21, size = 8.0, fill = "#DDE1E3", color = pal$grey, stroke = 0.5) +
  geom_text(data = black_nodes, aes(x = x, y = y, label = label),
            size = 4.0, color = pal$grey, fontface = "bold") +
  annotate("text", x = 0.14, y = 0.625, label = "u", size = 3.3,
           color = pal$muted, fontface = "bold") +
  geom_segment(data = perturb_df,
               aes(x = x, y = y, xend = xend, yend = yend, color = color),
               arrow = arrow(length = unit(0.014, "npc"), type = "closed"),
               linewidth = 0.7) +
  geom_tile(data = grid_pts, aes(x = x, y = y, fill = fill_code),
            width = 0.019, height = 0.034, color = "white", linewidth = 0.25) +
  scale_fill_identity() +
  annotate("text", x = 0.14, y = 0.245, label = "steady-state\nprofile matrix X",
           size = 2.55, color = pal$muted, lineheight = 0.95) +
  annotate("text", x = 0.38, y = 0.55,
           label = "dx_j/dt = mu_j +\nsum_i f_ji(x_i) + u_j",
           size = 2.75, lineheight = 0.95, color = pal$ink) +
  annotate("text", x = 0.38, y = 0.445,
           label = "steady state:\n0 = mu_j + sum_i f_ji(x_i) + u_j",
           size = 2.55, lineheight = 0.95, color = pal$muted) +
  annotate("text", x = 0.38, y = 0.335,
           label = "basis functions keep f_ji(0) = 0",
           size = 2.65, color = pal$muted) +
  annotate("rect", xmin = 0.565, xmax = 0.685, ymin = 0.51, ymax = 0.61,
           fill = "#ECF5EF", color = "#BAD8C6", linewidth = 0.35) +
  annotate("rect", xmin = 0.565, xmax = 0.685, ymin = 0.39, ymax = 0.49,
           fill = "#F7EFE0", color = "#DDC48A", linewidth = 0.35) +
  annotate("text", x = 0.625, y = 0.56, label = "between-node\ngroup sparsity",
           size = 2.55, color = pal$green, fontface = "bold", lineheight = 0.95) +
  annotate("text", x = 0.625, y = 0.44, label = "within-edge\nfunction selection",
           size = 2.55, color = pal$gold, fontface = "bold", lineheight = 0.95) +
  annotate("text", x = 0.625, y = 0.30, label = "not only prediction:\nrecover directed coupling",
           size = 2.45, color = pal$muted, lineheight = 0.95) +
  geom_curve(data = edge_df,
             aes(x = x, y = y, xend = xend, yend = yend,
                 color = color, linetype = linetype),
             curvature = 0.18,
             arrow = arrow(length = unit(0.014, "npc"), type = "closed"),
             linewidth = 0.75, inherit.aes = FALSE) +
  scale_linetype_identity() +
  geom_point(data = node_df, aes(x = x, y = y, size = size, fill = fill),
             shape = 21, color = "white", stroke = 0.8) +
  scale_size_identity() +
  geom_text(data = node_df, aes(x = x, y = y, label = name),
            size = 2.7, color = "white", fontface = "bold") +
  annotate("text", x = 0.84, y = 0.655, label = "core regulator",
           size = 2.45, color = pal$red, fontface = "bold") +
  annotate("segment", x = 0.84, xend = 0.84, y = 0.635, yend = 0.595,
           arrow = arrow(length = unit(0.01, "npc"), type = "closed"),
           color = pal$red, linewidth = 0.4) +
  annotate("text", x = 0.895, y = 0.235,
           label = "solid red = positive local effect\ndashed blue = negative local effect",
           size = 2.35, color = pal$muted, lineheight = 0.95) +
  geom_text(data = main_text, aes(x = x, y = y, label = label),
            size = 2.35, color = pal$muted, lineheight = 0.95) +
  annotate("text", x = 0.50, y = 0.965,
           label = "PSS-Net learns interpretable coupling networks from perturbed steady states",
           size = 4.5, color = pal$ink, fontface = "bold") +
  annotate("text", x = 0.50, y = 0.935,
           label = "Perturbations turn equilibrium observations into constraints on a sparse additive dynamical system",
           size = 2.85, color = pal$muted) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0.08, 0.985), expand = FALSE, clip = "off") +
  theme_void(base_size = 10) +
  theme(
    plot.background = element_rect(fill = pal$light, color = NA),
    panel.background = element_rect(fill = pal$light, color = NA),
    plot.margin = margin(12, 12, 12, 12)
  )

pdf_file <- file.path(out_dir, "pss_concept_diagram.pdf")
png_file <- file.path(out_dir, "pss_concept_diagram.png")

ggsave(pdf_file, p, width = 10.8, height = 6.4, units = "in", device = "pdf")
ggsave(png_file, p, width = 10.8, height = 6.4, units = "in", dpi = 320, device = ragg::agg_png)

message("Wrote: ", pdf_file)
message("Wrote: ", png_file)

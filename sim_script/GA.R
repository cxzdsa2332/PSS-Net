################################################################################
# GA.R -- graphical abstract for PSS-Net
#
# Purpose:
#   Draw a manuscript graphical abstract that communicates the core PSS-Net
#   workflow in a domain-general way:
#     complex system -> perturbations -> steady states -> sparse nonlinear
#     inference -> interpretable directed coupling network.
#
# Inputs:
#   none
#
# Outputs:
#   results/figure/GA.pdf
#   results/figure/GA.png
#
# Notes:
#   This script intentionally uses only base R graphics/grid so that the figure
#   is reproducible without extra plotting dependencies.
################################################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(grid)
})

out_dir <- "results/figure"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_pdf <- file.path(out_dir, "GA.pdf")
out_png <- file.path(out_dir, "GA.png")

## ------------------------------------------------------------------- palette ----
pal <- list(
  bg       = "#F7FAF7",
  ink      = "#26343D",
  muted    = "#66727C",
  line     = "#CFD8D8",
  panel    = "#FFFFFF",
  panel2   = "#FDFBF5",
  teal     = "#2B8C9A",
  teal2    = "#8FC9C5",
  blue     = "#446AAE",
  blue2    = "#DDE7F7",
  green    = "#3A9A68",
  green2   = "#E7F2EA",
  gold     = "#E2A529",
  gold2    = "#FBF1D8",
  red      = "#C24A4D",
  red2     = "#F7DEDC",
  purple   = "#7B5CA7",
  purple2  = "#ECE4F4",
  greybox  = "#F1F4F3"
)

## ------------------------------------------------------------------ helpers ----
gpar_text <- function(col = pal$ink, fontsize = 10, fontface = "plain",
                      lineheight = 0.95) {
  gpar(col = col, fontsize = fontsize, fontface = fontface,
       lineheight = lineheight, fontfamily = "sans")
}

box <- function(x, y, w, h, fill = pal$panel, col = pal$line, lwd = 1.2,
                r = 0.012) {
  grid.roundrect(
    x = x, y = y, width = w, height = h,
    r = unit(r, "npc"),
    gp = gpar(fill = fill, col = col, lwd = lwd)
  )
}

label <- function(txt, x, y, size = 10, col = pal$ink, face = "plain",
                  just = "centre", lineheight = 0.95, rot = 0) {
  grid.text(
    txt, x = x, y = y, just = just, rot = rot,
    gp = gpar_text(col = col, fontsize = size, fontface = face,
                   lineheight = lineheight)
  )
}

arr <- function(x1, y1, x2, y2, col = pal$muted, lwd = 2.2,
                alpha = 1, curvature = 0, lty = 1,
                arrow_len = 0.018) {
  gp <- gpar(col = adjustcolor(col, alpha.f = alpha), lwd = lwd,
             lineend = "round", linejoin = "round", lty = lty)
  if (abs(curvature) < 1e-8) {
    grid.lines(
      x = unit(c(x1, x2), "npc"), y = unit(c(y1, y2), "npc"),
      arrow = arrow(type = "closed", length = unit(arrow_len, "npc")),
      gp = gp
    )
  } else {
    grid.curve(
      x1 = x1, y1 = y1, x2 = x2, y2 = y2, default.units = "npc",
      curvature = curvature, angle = 90,
      arrow = arrow(type = "closed", length = unit(arrow_len, "npc")),
      gp = gp
    )
  }
}

node <- function(x, y, r = 0.014, fill = pal$teal, col = "white",
                 txt = NULL, txt_col = "white", size = 8) {
  grid.circle(x = x, y = y, r = r,
              gp = gpar(fill = fill, col = col, lwd = 1.2))
  if (!is.null(txt)) label(txt, x, y, size = size, col = txt_col, face = "bold")
}

mini_axis <- function(x0, y0, w, h, col = "#D8DEDE") {
  grid.lines(unit(c(x0, x0), "npc"), unit(c(y0, y0 + h), "npc"),
             gp = gpar(col = col, lwd = 0.8))
  grid.lines(unit(c(x0, x0 + w), "npc"), unit(c(y0, y0), "npc"),
             gp = gpar(col = col, lwd = 0.8))
}

draw_network <- function(cx, cy, scale = 1, alpha = 1, labels = TRUE) {
  pts <- data.frame(
    id = c("x1", "x2", "x3", "x4", "x5", "x6"),
    x = cx + scale * c(-0.045, 0.040, -0.006, -0.052, 0.055, 0.015),
    y = cy + scale * c( 0.055, 0.047, -0.005, -0.055, -0.050, -0.085),
    fill = c(pal$blue, pal$red, pal$green, pal$gold, pal$teal, "#8C9AA3")
  )
  edge <- function(a, b, col, lty = 1, curv = 0.15, lwd = 2.2) {
    pa <- pts[pts$id == a, ]; pb <- pts[pts$id == b, ]
    arr(pa$x, pa$y, pb$x, pb$y, col = col, lwd = lwd * scale,
        alpha = alpha, curvature = curv, lty = lty)
  }
  edge("x1", "x2", pal$red, 1, 0.12)
  edge("x2", "x3", pal$green, 1, -0.10)
  edge("x3", "x5", pal$red, 1, 0.12)
  edge("x3", "x4", pal$blue, 1, -0.18)
  edge("x4", "x6", pal$blue, 1, 0.18)
  edge("x2", "x1", pal$blue, 2, -0.12, 1.8)
  for (i in seq_len(nrow(pts))) {
    node(pts$x[i], pts$y[i], r = 0.014 * scale, fill = pts$fill[i],
         txt = if (labels) pts$id[i] else NULL, size = 7 * scale)
  }
}

draw_recovered_network <- function(cx, cy, scale = 1) {
  pts <- data.frame(
    id = c("A", "B", "C", "D", "E", "F"),
    x = cx + scale * c(-0.052, 0.045, -0.006, -0.060, 0.062, 0.030),
    y = cy + scale * c( 0.060, 0.052, -0.005, -0.055, -0.048, -0.092),
    fill = c(pal$red, pal$blue, pal$green, pal$gold, pal$teal, "#8C9AA3")
  )
  edge <- function(a, b, col, curv = 0.10, lty = 1, lwd = 2.4) {
    pa <- pts[pts$id == a, ]; pb <- pts[pts$id == b, ]
    dx <- pb$x - pa$x
    dy <- pb$y - pa$y
    len <- sqrt(dx^2 + dy^2)
    shrink <- 0.018 * scale
    sx <- pa$x + dx / len * shrink
    sy <- pa$y + dy / len * shrink
    ex <- pb$x - dx / len * shrink
    ey <- pb$y - dy / len * shrink
    arr(sx, sy, ex, ey, col = col, lwd = lwd * scale,
        curvature = curv, lty = lty, arrow_len = 0.022)
  }
  edge("A", "B", pal$red, 0.08, 1, 2.5)
  edge("B", "C", pal$green, 0.12, 1, 2.4)
  edge("C", "E", pal$red, -0.10, 1, 2.4)
  edge("C", "D", pal$blue, 0.12, 2, 2.3)
  edge("D", "F", pal$blue, -0.18, 2, 2.1)
  for (i in seq_len(nrow(pts))) {
    node(pts$x[i], pts$y[i], r = ifelse(pts$id[i] == "A", 0.019, 0.015) * scale,
         fill = pts$fill[i], txt = pts$id[i], size = 7.2 * scale)
  }
}

draw_heatmap <- function(x0, y0, nrow = 6, ncol = 6, cell = 0.016,
                         gap = 0.003) {
  vals <- matrix(c(
    0.15, 0.65, 0.42, 0.27, 0.84, 0.34,
    0.28, 0.39, 0.76, 0.51, 0.31, 0.68,
    0.55, 0.20, 0.35, 0.71, 0.43, 0.22,
    0.80, 0.48, 0.23, 0.39, 0.62, 0.29,
    0.33, 0.73, 0.58, 0.18, 0.44, 0.77,
    0.22, 0.36, 0.64, 0.50, 0.28, 0.46
  ), nrow = nrow, byrow = TRUE)
  ramp <- colorRampPalette(c("#E9EED7", "#F1DCA3", "#DD8E64", "#B95B5D"))
  cols <- ramp(100)
  for (r in seq_len(nrow)) for (c in seq_len(ncol)) {
    grid.rect(
      x = x0 + (c - 1) * (cell + gap),
      y = y0 + (nrow - r) * (cell + gap),
      width = cell, height = cell,
      just = c("left", "bottom"),
      gp = gpar(fill = cols[pmax(1, pmin(100, round(vals[r, c] * 100)))],
                col = "white", lwd = 0.5)
    )
  }
}

draw_bars <- function(x0, y0, w, h) {
  axis_col <- "#D6DDDD"
  grid.lines(unit(c(x0, x0 + w), "npc"), unit(c(y0, y0), "npc"),
             gp = gpar(col = axis_col, lwd = 0.8))
  groups <- 5
  heights <- matrix(c(
    0.45, 0.82, 0.32,
    0.20, 0.77, 0.36,
    0.42, 0.31, 0.34,
    0.88, 0.71, 0.35,
    0.23, 0.39, 0.31
  ), nrow = groups, byrow = TRUE)
  cols <- c(pal$blue, pal$red, pal$green)
  bw <- w / (groups * 4.4)
  for (g in seq_len(groups)) for (j in 1:3) {
    xx <- x0 + (g - 1) * w / groups + (j - 1) * bw * 1.15 + bw * 0.7
    grid.rect(
      x = xx, y = y0, width = bw, height = h * heights[g, j],
      just = c("left", "bottom"),
      gp = gpar(fill = adjustcolor(cols[j], alpha.f = 0.82),
                col = NA)
    )
  }
  label("conditions", x0 + w / 2, y0 - 0.025, size = 6.6, col = pal$muted)
}

draw_traj <- function(x0, y0, w, h) {
  mini_axis(x0, y0, w, h)
  xs <- seq(0, 1, length.out = 60)
  curves <- list(
    list(y = 0.20 + 0.62 * (1 - exp(-4 * xs)), col = pal$red),
    list(y = 0.75 - 0.43 * (1 - exp(-5 * xs)), col = pal$blue),
    list(y = 0.35 + 0.10 * sin(2.8 * pi * xs) * exp(-3 * xs) + 0.18 * xs,
         col = pal$green)
  )
  for (cc in curves) {
    grid.lines(unit(x0 + xs * w, "npc"), unit(y0 + cc$y * h, "npc"),
               gp = gpar(col = cc$col, lwd = 1.8, lineend = "round"))
    grid.lines(unit(c(x0 + 0.78 * w, x0 + w), "npc"),
               unit(rep(y0 + tail(cc$y, 1) * h, 2), "npc"),
               gp = gpar(col = cc$col, lwd = 1.8, lty = 2))
  }
}

draw_basis_curves <- function(x0, y0, w, h) {
  mini_axis(x0, y0, w, h)
  xs <- seq(0, 1, length.out = 60)
  ys1 <- 0.15 + 0.68 * xs
  ys2 <- 0.08 + 0.76 * xs^2
  ys3 <- 0.12 + 0.65 * (1 - exp(-3.5 * xs))
  grid.lines(unit(x0 + xs * w, "npc"), unit(y0 + ys1 * h, "npc"),
             gp = gpar(col = pal$blue, lwd = 1.5))
  grid.lines(unit(x0 + xs * w, "npc"), unit(y0 + ys2 * h, "npc"),
             gp = gpar(col = pal$red, lwd = 1.5))
  grid.lines(unit(x0 + xs * w, "npc"), unit(y0 + ys3 * h, "npc"),
             gp = gpar(col = pal$green, lwd = 1.5))
  label("basis library", x0 + w / 2, y0 - 0.020, size = 6.8, col = pal$muted)
}

draw_function_inset <- function(x0, y0, w, h, col = pal$red,
                                title = "f_ji(x)") {
  box(x0 + w / 2, y0 + h / 2, w, h, fill = "#FFFFFF", col = "#D8E0E0", lwd = 0.8)
  mini_axis(x0 + 0.18 * w, y0 + 0.22 * h, 0.65 * w, 0.55 * h)
  xs <- seq(0, 1, length.out = 60)
  ys <- 0.15 + 0.70 * (xs - 0.25)^2 + 0.15 * xs
  ys <- (ys - min(ys)) / diff(range(ys)) * 0.55 + 0.20
  grid.lines(
    unit(x0 + 0.18 * w + xs * 0.65 * w, "npc"),
    unit(y0 + 0.22 * h + ys * 0.55 * h, "npc"),
    gp = gpar(col = col, lwd = 1.6)
  )
  label(title, x0 + 0.52 * w, y0 + 0.83 * h, size = 6.5, col = col,
        face = "bold")
}

draw_card <- function(x, y, w, h, title, subtitle, col, fill) {
  box(x, y, w, h, fill = fill, col = adjustcolor(col, alpha.f = 0.40), lwd = 1.1)
  label(title, x, y + h * 0.16, size = 8.8, col = col, face = "bold")
  label(subtitle, x, y - h * 0.13, size = 6.8, col = pal$muted,
        lineheight = 0.86)
}

draw_ga <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = pal$bg, col = NA))

  ## Title
  label("PSS-Net learns interpretable coupling networks from perturbed steady states",
        0.50, 0.955, size = 19, col = pal$ink, face = "bold")
  label("External inputs convert equilibrium measurements into constraints on sparse nonlinear edge functions",
        0.50, 0.918, size = 10.5, col = pal$muted)

  ## Panel coordinates
  px <- c(0.045, 0.285, 0.525, 0.765)
  py <- 0.155
  pw <- 0.195
  ph <- 0.690
  pcent <- px + pw / 2

  for (i in seq_along(px)) {
    box(px[i] + pw / 2, py + ph / 2, pw, ph,
        fill = if (i %% 2 == 1) pal$panel else pal$panel2,
        col = "#D5DEDE", lwd = 1.2)
  }

  ## Flow arrows between panels
  arr(px[1] + pw + 0.010, 0.505, px[2] - 0.012, 0.505, col = "#6B777E", lwd = 2.6)
  arr(px[2] + pw + 0.010, 0.505, px[3] - 0.012, 0.505, col = "#6B777E", lwd = 2.6)
  arr(px[3] + pw + 0.010, 0.505, px[4] - 0.012, 0.505, col = "#6B777E", lwd = 2.6)

  ## Panel 1
  label("1  perturb", pcent[1], py + ph - 0.040, size = 12.5, col = pal$teal,
        face = "bold")
  label("a complex system", pcent[1], py + ph - 0.073, size = 8.5,
        col = pal$muted)
  draw_network(pcent[1], py + 0.460, scale = 1.05, alpha = 0.85, labels = TRUE)
  label("hidden directed coupling", pcent[1], py + 0.590, size = 7.5,
        col = pal$muted)
  ## Perturbation arrows: a separate input symbol below the unknown network.
  ## They intentionally do not touch nodes, avoiding the impression that these
  ## arrows are inferred network edges.
  arrow_cols <- c(pal$blue, pal$gold, pal$green)
  arrow_x <- pcent[1] + c(-0.055, 0.000, 0.055)
  for (i in 1:3) {
    arr(arrow_x[i], py + 0.255, arrow_x[i], py + 0.350,
        col = arrow_cols[i], lwd = 2.6, alpha = 0.95, arrow_len = 0.015)
  }
  box(pcent[1], py + 0.160, 0.118, 0.050, fill = pal$greybox, col = "#DEE5E5")
  label("known inputs  u", pcent[1], py + 0.160, size = 8.2,
        col = pal$ink, face = "bold")
  label("not labels, not trajectories:\ncontrolled experimental pushes",
        pcent[1], py + 0.065, size = 7.0, col = pal$muted, lineheight = 0.86)

  ## Panel 2
  label("2  measure", pcent[2], py + ph - 0.040, size = 12.5, col = pal$blue,
        face = "bold")
  label("post-perturbation steady states", pcent[2], py + ph - 0.075, size = 8.2,
        col = pal$muted)
  draw_traj(px[2] + 0.036, py + 0.430, 0.124, 0.140)
  draw_bars(px[2] + 0.030, py + 0.270, 0.138, 0.125)
  label("profile matrix  X*", pcent[2], py + 0.230, size = 8.2, col = pal$ink,
        face = "bold")
  draw_heatmap(px[2] + 0.047, py + 0.075, nrow = 6, ncol = 6, cell = 0.014,
               gap = 0.0028)
  label("rows = conditions\ncolumns = system variables",
        pcent[2], py + 0.040, size = 6.7, col = pal$muted, lineheight = 0.86)

  ## Panel 3
  label("3  infer", pcent[3], py + ph - 0.040, size = 12.5, col = pal$green,
        face = "bold")
  label("sparse nonlinear steady-state effects", pcent[3], py + ph - 0.075,
        size = 8.2, col = pal$muted)
  box(pcent[3], py + 0.535, 0.155, 0.102, fill = "#F6FBF8", col = "#BBDAC9")
  label("0 = mu_j + f_jj(x_j*)\n    + sum_i f_ji(x_i*) + u_j",
        pcent[3], py + 0.548, size = 8.2, col = pal$ink, face = "bold",
        lineheight = 0.90)
  label("steady-state equation", pcent[3], py + 0.490, size = 7.0,
        col = pal$muted)
  draw_basis_curves(px[3] + 0.043, py + 0.315, 0.115, 0.105)
  draw_card(pcent[3], py + 0.225, 0.150, 0.070,
            "source groups", "which variables influence node j", pal$green,
            pal$green2)
  draw_card(pcent[3], py + 0.130, 0.150, 0.070,
            "basis terms", "linear / curved / saturating effects", pal$gold,
            pal$gold2)
  label("ADSIHT selects a compact model", pcent[3], py + 0.057,
        size = 7.2, col = pal$muted)

  ## Panel 4
  label("4  recover", pcent[4], py + ph - 0.040, size = 12.5, col = pal$red,
        face = "bold")
  label("directed map + model-based functions", pcent[4], py + ph - 0.075, size = 8.2,
        col = pal$muted)
  draw_recovered_network(pcent[4], py + 0.430, scale = 1.18)
  label("directed edges + local signs", pcent[4], py + 0.575, size = 7.6,
        col = pal$muted)
  draw_function_inset(px[4] + 0.028, py + 0.178, 0.070, 0.084, col = pal$red,
                      title = "f_ji(x)")
  draw_function_inset(px[4] + 0.107, py + 0.178, 0.070, 0.084, col = pal$blue,
                      title = "J_ji")
  label("f_ji(x) is drawn after choosing\na steady-state model / basis\n(e.g. gLV or additive library)",
        pcent[4], py + 0.070, size = 6.6, col = pal$muted, lineheight = 0.86)

  ## Bottom take-home ribbon
  box(0.50, 0.065, 0.78, 0.070, fill = "#FFFFFF", col = "#DDE5E5", lwd = 1.0,
      r = 0.020)
  label("Core idea: perturbations produce steady states; PSS-Net recovers a directed signed map; specified models add edge functions.",
        0.50, 0.067, size = 9.2, col = pal$ink, face = "bold")

  ## Tiny legend
  grid.lines(unit(c(0.812, 0.852), "npc"), unit(c(0.086, 0.086), "npc"),
             gp = gpar(col = pal$red, lwd = 2.4))
  label("positive", 0.874, 0.086, size = 6.8, col = pal$muted, just = "left")
  grid.lines(unit(c(0.812, 0.852), "npc"), unit(c(0.060, 0.060), "npc"),
             gp = gpar(col = pal$blue, lwd = 2.4, lty = 2))
  label("negative", 0.874, 0.060, size = 6.8, col = pal$muted, just = "left")
}

## --------------------------------------------------------------------- export ----
pdf(out_pdf, width = 14.5, height = 8.0, useDingbats = FALSE)
draw_ga()
dev.off()

png(out_png, width = 2600, height = 1430, res = 220)
draw_ga()
dev.off()

cat("Saved graphical abstract:\n")
cat(" -", out_pdf, "\n")
cat(" -", out_png, "\n")

rm(list = ls())

################################################################################
# Fig3c_structure_dependence.R — PSS-Net on scale-free vs ER (Erdos-Renyi)
#
# Purpose: compare PSS-Net network recovery across two contrasting topologies,
#          using the SAME estimator, basis and sample-budget grid, so the only
#          thing that changes is the degree structure:
#            scalefree — preferential attachment, power-law out-degree, a few hub
#                        sources influencing many targets, variable in-degree;
#            er        — homogeneous Erdos-Renyi-like graph, fixed in-degree,
#                        roughly Poisson out-degree, no hubs.
#          Both topologies share the edge magnitude range, the linear steady
#          state, noise level and the N / (avg_in * log p) budget axis, and both
#          use the expected in-degree avg_in so edge density is matched.
#
#          For each topology two PSS-Net variants are run:
#            nodewise — p independent ADSIHT/DSIC solves (default);
#            joint    — one block-diagonal solve, group = rep(1:(p*p), each = M)
#                       (CLAUDE.md rule). Global DSIC pooling is expected to help
#                       on heterogeneous (hub) structure.
#          Metrics: edge recovery (Pr / Re / MCC) and hub identification
#          (estimated vs true out-degree Spearman rho, top-k hub hit rate).
#
# Reuse: the scale-free generator and the inference / metric helpers are taken
#        verbatim from sim_script/03_robustness_benchmarks/pss_net_scalefree.R;
#        only make_er and the topology loop are new.
#
# Input:   none
# Output:  results/sim_results/Fig3c_structure_dependence.csv
################################################################################

suppressMessages(library(ADSIHT))
set.seed(1)
M_ord <- 2L

# scale-free: preferential attachment gives a power-law out-degree; in-degree
# is variable (1 + Poisson). pa_power sharpens the attachment (source chosen
# with probability proportional to out-degree^pa_power): pa_power = 2 produces
# pronounced hubs (max out-degree ~20 at p = 50), making the network strongly
# heterogeneous. This is the regime where the joint block-diagonal solve has the
# most shared-source structure to pool across targets, so its advantage over
# node-wise (edge MCC and especially hub recovery) is clearest; with the plain
# pa_power = 1 attachment the two schemes are close and noisy.
make_scalefree <- function(p, avg_in = 2, seed = 1, pa_power = 2) {
  set.seed(seed)
  A <- matrix(0, p, p)
  outdeg <- rep(1, p)                                  # +1 smoothing, pref. weight
  for (j in seq_len(p)) {
    k_j <- min(1 + rpois(1, avg_in - 1), p - 1)        # variable in-degree
    chosen <- integer(0)
    for (t in seq_len(k_j)) {
      w <- outdeg^pa_power; w[c(j, chosen)] <- 0       # no self-loop / repeat
      i <- sample.int(p, 1, prob = w)                  # P(source) prop. to out-degree^power
      chosen <- c(chosen, i)
      A[j, i] <- runif(1, 0.15, 0.35) * sample(c(-1, 1), 1)
      outdeg[i] <- outdeg[i] + 1
    }
  }
  gamma <- rowSums(abs(A)) + runif(p, 1.0, 1.5)
  r <- runif(p, 0.8, 1.5)
  list(p = p, A = A, gamma = gamma, r = r, adj = (A != 0) * 1,
       outdeg_true = colSums(A != 0))
}

# ER / homogeneous control: every target has exactly avg_in incoming edges drawn
# from uniformly random sources (no preferential attachment), so the out-degree
# is roughly Poisson and there are no hubs. Edge magnitude / gamma / r are kept
# identical to make_scalefree so only the degree structure differs.
make_er <- function(p, avg_in = 2, seed = 1) {
  set.seed(seed)
  A <- matrix(0, p, p)
  for (j in seq_len(p)) {
    src <- sample(setdiff(seq_len(p), j), avg_in)
    A[j, src] <- runif(avg_in, 0.15, 0.35) * sample(c(-1, 1), avg_in, TRUE)
  }
  gamma <- rowSums(abs(A)) + runif(p, 1.0, 1.5)
  r <- runif(p, 0.8, 1.5)
  list(p = p, A = A, gamma = gamma, r = r, adj = (A != 0) * 1,
       outdeg_true = colSums(A != 0))
}

make_system <- function(topology, p, avg_in, seed) {
  if (topology == "scalefree") make_scalefree(p, avg_in, seed)
  else make_er(p, avg_in, seed)
}

ss_lin <- function(sys, U) t(apply(U, 1, function(u)
  as.numeric(solve(diag(sys$gamma) - sys$A, sys$r + u))))

psi_row <- function(xv) as.vector(sapply(xv, function(x) x^(seq_len(M_ord))))

build_design <- function(U, X) {
  p <- ncol(X)
  Psi <- t(apply(X, 1, psi_row))
  Psi_c <- sweep(Psi, 2, colMeans(Psi))
  sdv <- apply(Psi_c, 2, sd); sdv[sdv < 1e-10] <- 1e-10
  Psi_cs <- sweep(Psi_c, 2, sdv, FUN = "/")
  rhs <- sapply(seq_len(p), function(j) -(U[, j] - mean(U[, j])))
  list(Psi_cs = Psi_cs, rhs = rhs, p = p)
}

adj_from_theta <- function(th, p) {
  gn <- sapply(seq_len(p), function(i)
    sqrt(sum(th[((i - 1) * M_ord + 1):(i * M_ord)]^2)))
  gn[!is.finite(gn)] <- 0                              # guard NaN coefficients
  as.integer(gn >= 1e-8)
}

infer_nodewise <- function(d) {
  p <- d$p; group <- rep(seq_len(p), each = M_ord)
  adj <- matrix(0, p, p)
  for (j in seq_len(p)) {
    fit <- tryCatch(ADSIHT(d$Psi_cs, matrix(d$rhs[, j]), group,
                           ic.type = "dsic"), error = function(e) NULL)
    if (is.null(fit)) next
    adj[j, ] <- adj_from_theta(fit$beta[, which.min(fit$ic)], p)
  }
  diag(adj) <- 0; adj
}

infer_joint <- function(d) {
  p <- d$p
  Xbig <- kronecker(diag(p), d$Psi_cs)
  Ybig <- as.vector(d$rhs)
  group <- rep(seq_len(p * p), each = M_ord)          # v0.1.txt: p*p groups
  fit <- tryCatch(ADSIHT(Xbig, matrix(Ybig), group, ic.type = "dsic"),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  beta <- fit$beta[, which.min(fit$ic)]; pM <- p * M_ord
  adj <- matrix(0, p, p)
  for (j in seq_len(p))
    adj[j, ] <- adj_from_theta(beta[((j - 1) * pM + 1):(j * pM)], p)
  diag(adj) <- 0; adj
}

mcc_of <- function(est, true) {
  off <- which(row(true) != col(true))
  e <- est[off]; t <- true[off]
  TP <- sum(e == 1 & t == 1); FP <- sum(e == 1 & t == 0)
  TN <- sum(e == 0 & t == 0); FN <- sum(e == 0 & t == 1)
  pr <- ifelse(TP + FP == 0, 0, TP / (TP + FP))
  re <- ifelse(TP + FN == 0, 0, TP / (TP + FN))
  den <- sqrt(as.numeric(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  c(Pr = pr, Re = re, MCC = ifelse(den == 0, 0, (TP * TN - FP * FN) / den))
}

# hub identification: out-degree Spearman + top-k hit (k = top 10% true out-deg)
hub_of <- function(adj_est, sys) {
  out_est <- colSums(adj_est); out_true <- sys$outdeg_true
  k <- max(3, round(0.1 * sys$p))
  true_top <- order(out_true, decreasing = TRUE)[seq_len(k)]
  est_top <- order(out_est, decreasing = TRUE)[seq_len(k)]
  rho <- suppressWarnings(cor(out_est, out_true, method = "spearman"))
  c(rho = ifelse(is.na(rho), 0, rho),
    topk_hit = length(intersect(true_top, est_top)) / k)
}

## ----------------------------------------------------------- main loop ----
p <- 50L; avg_in <- 2
k_grid <- c(2, 3, 4, 5)
R <- 10
sigma <- 0.03
topologies <- c("scalefree", "er")

rows <- list()
base <- avg_in * log(p)
N_grid <- unique(ceiling(k_grid * base))
for (topology in topologies) {
  for (seed in seq_len(R)) {
    sys <- make_system(topology, p, avg_in = avg_in, seed = 1000 + seed)
    max_out <- max(sys$outdeg_true); n_edge <- sum(sys$adj)
    for (N in N_grid) {
      U <- matrix(runif(N * p, -0.3, 0.5), N, p); U[1, ] <- 0
      X <- ss_lin(sys, U)
      ok <- apply(X, 1, function(rr) all(is.finite(rr)) && all(rr > 0))
      Xo <- (X + matrix(rnorm(length(X), sd = sigma), nrow(X)))[ok, , drop = FALSE]
      d <- build_design(U[ok, , drop = FALSE], Xo)
      adj_nw <- infer_nodewise(d); adj_jt <- infer_joint(d)
      m_nw <- mcc_of(adj_nw, sys$adj)
      m_jt <- if (is.null(adj_jt)) c(Pr = NA, Re = NA, MCC = NA) else mcc_of(adj_jt, sys$adj)
      h_nw <- hub_of(adj_nw, sys)
      h_jt <- if (is.null(adj_jt)) c(rho = NA, topk_hit = NA) else hub_of(adj_jt, sys)
      rows[[length(rows) + 1]] <- data.frame(
        topology = topology, p = p, N = N, seed = seed,
        n_edge = n_edge, max_outdeg = max_out,
        MCC_nodewise = m_nw["MCC"], MCC_joint = m_jt["MCC"],
        rho_nodewise = h_nw["rho"], rho_joint = h_jt["rho"],
        topk_nodewise = h_nw["topk_hit"], topk_joint = h_jt["topk_hit"])
      cat(sprintf("%-9s seed=%d N=%2d | maxOut=%2d | MCC nw=%.3f jt=%.3f | hubRho nw=%.2f jt=%.2f\n",
                  topology, seed, N, max_out, m_nw["MCC"], m_jt["MCC"],
                  h_nw["rho"], h_jt["rho"]))
    }
  }
}
df <- do.call(rbind, rows); rownames(df) <- NULL

dir.create("results/sim_results", showWarnings = FALSE, recursive = TRUE)
write.csv(df, "results/sim_results/Fig3c_structure_dependence.csv", row.names = FALSE)

agg <- aggregate(cbind(MCC_nodewise, MCC_joint, rho_nodewise, rho_joint,
                       topk_nodewise, topk_joint) ~ topology + N, df, mean)
cat("\n===== PSS-Net structure dependence (p=", p, ", avg_in=", avg_in, ", ",
    R, " seeds) =====\n", sep = "")
print(agg, row.names = FALSE)
cat("\nSaved: results/sim_results/Fig3c_structure_dependence.csv\n")

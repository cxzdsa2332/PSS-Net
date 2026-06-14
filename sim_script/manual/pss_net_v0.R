rm(list = ls())

################################################################################
# pss_net_v0.R  —  PSS-Net: 100-OTU microbial community, 30 PSS conditions
#
# Scenario: realistic underdetermined regime (N=30, p=100, M=1)
#   - Sparse random GLV, strong self-regulation for ODE stability
#   - 30 perturbation conditions simulated to steady state
#   - Inference: ADSIHT (group lasso omitted at this scale)
#   - No ground-truth evaluation — treated as real-data analysis
#
# Model:  dx_j/dt = μ_j + Σ_i A_{ji}·x_i − γ_j·x_j + u_j   (true GLV)
# Basis:  ψ(x) = x  (M=1, linear)  →  Ψ ∈ R^{N×p}
# At SS:  Ψ_c · θ_j = −u_{c,j}
################################################################################

library(deSolve)
library(ADSIHT)
library(ggplot2)
library(igraph)

set.seed(2025)

# ── Parameters ────────────────────────────────────────────────────────────────
n_sp  <- 100
N_cond <- 30
M      <- 2
sigma  <- 0.03

otu_names <- paste0("OTU", sprintf("%03d", seq_len(n_sp)))

cat("================================================================\n")
cat("  PSS-Net v0: 100-OTU microbial community\n")
cat(sprintf("  %d OTUs | %d PSS conditions | M=%d | σ=%.2f\n",
            n_sp, N_cond, M, sigma))
cat("================================================================\n\n")

# ── Sparse random GLV ────────────────────────────────────────────────────────
# Mean in-degree k_in=3; weak off-diagonal weights; strong self-regulation
# guarantees diagonal dominance → unique stable equilibrium per condition.
cat("Generating sparse random GLV network...\n")

k_in   <- 3
p_edge <- k_in / (n_sp - 1)

r_true   <- runif(n_sp, 0.5, 1.5)
gam_true <- runif(n_sp, 4.0, 6.0)   # strong self-regulation

A_true <- matrix(0, n_sp, n_sp)
for (j in seq_len(n_sp)) {
  srcs <- setdiff(which(runif(n_sp) < p_edge), j)
  for (i in srcs)
    A_true[j, i] <- sample(c(-1,1), 1) * runif(1, 0.05, 0.20)
}
n_edges <- sum(A_true != 0)
cat(sprintf("  Edges: %d  (mean in-degree %.1f, density %.2f%%)\n",
            n_edges, n_edges/n_sp, 100*n_edges/(n_sp*(n_sp-1))))

# Verify diagonal dominance: max|row sum off-diag| < min(gam)
max_offdiag <- max(rowSums(abs(A_true)))
cat(sprintf("  max|Σ_i |A_{ji}|| = %.3f  vs  min(gam) = %.3f  → %s\n",
            max_offdiag, min(gam_true),
            if (max_offdiag < min(gam_true)) "STABLE" else "WARNING"))

# ── ODE ───────────────────────────────────────────────────────────────────────
ode_func <- function(t, state, parms) {
  x <- pmax(state, 0)
  list(r_true + as.numeric(A_true %*% x) - gam_true * x + parms$u)
}

# ── Step 1: PSS simulation ────────────────────────────────────────────────────
cat("\nStep 1: Simulating PSS data...\n")

out_wt <- ode(rep(0.3, n_sp), c(0, 1e4), ode_func, list(u=rep(0,n_sp)),
              method="lsoda", rtol=1e-12, atol=1e-14)
x_wt   <- as.numeric(out_wt[nrow(out_wt), 2:(n_sp+1)])
cat(sprintf("  x_wt: range=[%.3f, %.3f]  mean=%.3f\n",
            min(x_wt), max(x_wt), mean(x_wt)))

# Perturbation magnitude bounded so system stays positive at SS
u_scale <- min(r_true) * 0.4
U_mat   <- matrix(runif(N_cond * n_sp, -u_scale, u_scale*2), N_cond, n_sp)
U_mat[1, ] <- 0

ss_mat <- matrix(NA, N_cond, n_sp)
for (k in seq_len(N_cond)) {
  tryCatch({
    out <- ode(x_wt, c(0, 3000), ode_func, list(u=U_mat[k,]),
               method="lsoda", rtol=1e-9, atol=1e-11)
    ss_mat[k, ] <- out[nrow(out), 2:(n_sp+1)]
  }, error=function(e) NULL)
}

ok    <- apply(ss_mat, 1, function(r) all(is.finite(r) & r > 0))
X_obs <- ss_mat[ok, ] + matrix(rnorm(sum(ok)*n_sp, 0, sigma), sum(ok), n_sp)
X_obs <- pmax(X_obs, 1e-6)
U_obs <- U_mat[ok, ]
N_ok  <- nrow(X_obs)
cat(sprintf("  Valid conditions: %d / %d\n\n", N_ok, N_cond))


colnames(X_obs) = otu_names
colnames(U_obs) = otu_names

#以上是生成一个模拟数据


# ── Step 2: Linear basis + centering + standardization ───────────────────────
cat(sprintf("Step 2: Design matrix (%d × %d)...\n", N_ok, n_sp*M))

# M=1: Ψ_{k,i} = x_i^(k)  →  same as X_obs
Psi <- X_obs   # N × p  (M=1 means one basis per OTU)

group   <- seq_len(n_sp)          # each OTU is its own group (size M=1)
Psi_bar <- colMeans(Psi)
Psi_c   <- sweep(Psi, 2, Psi_bar)
Psi_sd  <- pmax(apply(Psi_c, 2, sd), 1e-10)
Psi_cs  <- sweep(Psi_c, 2, Psi_sd, "/")
U_bar   <- colMeans(U_obs)
U_c     <- sweep(U_obs, 2, U_bar)
cat(sprintf("  kappa(Ψ_cs) = %.1f\n\n", kappa(Psi_cs, exact=FALSE)))

# ── Step 3: ADSIHT ───────────────────────────────────────────────────────────
cat("Step 3: ADSIHT inference (100 OTUs)...\n")

ALPHA_ads <- matrix(0, n_sp, n_sp)   # M=1: theta_{ji} is scalar
MU_ads    <- numeric(n_sp)

pb_step <- max(1L, n_sp %/% 5L)
for (j in seq_len(n_sp)) {
  if (j %% pb_step == 0) cat(sprintf("  OTU %3d / %d\n", j, n_sp))
  rhs_c <- -U_c[, j]
  fit   <- tryCatch(
    ADSIHT(Psi_cs, matrix(rhs_c), group, ic.type="dsic"),
    error = function(e) NULL)
  if (is.null(fit) || ncol(fit$beta) == 0) next
  best <- which.min(fit$ic)
  if (length(best) == 0 || !is.finite(fit$ic[best])) next
  ALPHA_ads[j, ] <- fit$beta[, best] / Psi_sd
  MU_ads[j]      <- -U_bar[j] - sum(Psi_bar * ALPHA_ads[j, ])
}
cat("  Done.\n\n")

# ── Step 4: Adjacency + Jacobian ─────────────────────────────────────────────
# M=1: group norm = |theta_{ji}|; Jacobian = theta_{ji} (linear)
tau_ads  <- 1e-10
adj_hat  <- (abs(ALPHA_ads) >= tau_ads) * 1L
diag(adj_hat) <- 0

A_hat <- ALPHA_ads   # M=1: Jacobian = coefficient directly
diag(A_hat) <- 0
g_hat <- -diag(ALPHA_ads)

n_inf  <- sum(adj_hat)
in_deg  <- colSums(adj_hat)
out_deg <- rowSums(adj_hat)

el <- which(adj_hat != 0, arr.ind=TRUE)
n_pos <- if (nrow(el)>0) sum(A_hat[el] > 0) else 0
n_neg <- if (nrow(el)>0) sum(A_hat[el] < 0) else 0

# ── Summary ───────────────────────────────────────────────────────────────────
cat("================================================================\n")
cat(sprintf("  PSS-Net v0: 100-OTU inferred network\n"))
cat(sprintf("  N=%d valid conditions | M=%d | τ=%.0e\n", N_ok, M, tau_ads))
cat(sprintf("  Inferred edges: %d  (+)=%d  (−)=%d\n", n_inf, n_pos, n_neg))
cat(sprintf("  In-degree:  min=%d max=%d mean=%.2f median=%.1f\n",
            min(in_deg), max(in_deg), mean(in_deg), median(in_deg)))
cat(sprintf("  Out-degree: min=%d max=%d mean=%.2f median=%.1f\n",
            min(out_deg), max(out_deg), mean(out_deg), median(out_deg)))
cat(sprintf("  μ̂ range: [%.3f, %.3f]\n", min(MU_ads), max(MU_ads)))
cat(sprintf("  γ̂ range: [%.3f, %.3f]\n", min(g_hat),  max(g_hat)))

if (nrow(el) > 0) {
  ew  <- A_hat[el]
  ord <- order(abs(ew), decreasing=TRUE)
  cat("\n  Top 15 inferred interactions (|A_hat_{ji}|):\n")
  cat(sprintf("  %-10s %-10s  %8s\n", "Target j","Source i","A_hat"))
  for (k in seq_len(min(15, nrow(el)))) {
    jj <- el[ord[k],1]; ii <- el[ord[k],2]
    cat(sprintf("  %-10s %-10s  %+8.4f  %s\n",
                otu_names[jj], otu_names[ii], A_hat[jj,ii],
                ifelse(A_hat[jj,ii]>0,"(+)","(−)")))
  }
}
cat("================================================================\n\n")

# ── Plot 1: Degree distribution ───────────────────────────────────────────────
deg_df <- rbind(
  data.frame(degree=in_deg,  type="In-degree  (regulators)"),
  data.frame(degree=out_deg, type="Out-degree (targets)")
)
p_deg <- ggplot(deg_df, aes(x=degree, fill=type)) +
  geom_histogram(binwidth=1, color="white", linewidth=0.25, alpha=0.85) +
  facet_wrap(~type, ncol=2, scales="free_y") +
  scale_fill_manual(values=c("In-degree  (regulators)"="steelblue3",
                             "Out-degree (targets)"="tomato3")) +
  labs(title="PSS-Net v0: Inferred degree distribution (100 OTUs)",
       subtitle=sprintf("%d edges | N=%d | M=%d", n_inf, N_ok, M),
       x="Degree", y="Count") +
  theme_bw(base_size=10) +
  theme(legend.position="none",
        strip.background=element_rect(fill="grey92"),
        strip.text=element_text(face="bold", size=9),
        plot.title=element_text(size=11, face="bold"),
        plot.subtitle=element_text(size=9, color="grey45"))
print(p_deg)

# ── Plot 2: Interaction strength distribution ─────────────────────────────────
if (nrow(el) > 0) {
  strength_df <- data.frame(
    A = A_hat[el],
    sign = ifelse(A_hat[el] > 0, "Promotion (+)", "Inhibition (−)")
  )
  p_str <- ggplot(strength_df, aes(x=A, fill=sign)) +
    geom_histogram(bins=25, color="white", linewidth=0.2, alpha=0.85) +
    geom_vline(xintercept=0, linetype="dashed", color="grey40", linewidth=0.5) +
    scale_fill_manual(values=c("Promotion (+)"="tomato3",
                               "Inhibition (−)"="steelblue3")) +
    labs(title="PSS-Net v0: Inferred interaction strength A_hat",
         subtitle=sprintf("n(+)=%d  n(−)=%d  |  Jacobian at x_wt", n_pos, n_neg),
         x=expression(hat(A)[ji]), y="Count", fill=NULL) +
    theme_bw(base_size=10) +
    theme(legend.position="bottom",
          plot.title=element_text(size=11, face="bold"),
          plot.subtitle=element_text(size=9, color="grey45"))
  print(p_str)
}

# ── Plot 3: Adjacency heatmap (reordered by in-degree) ────────────────────────
ord_nodes  <- order(in_deg, decreasing=TRUE)
Aheat_ord  <- A_hat[ord_nodes, ord_nodes] * adj_hat[ord_nodes, ord_nodes]
ticks <- seq(1, n_sp, by=10)

heat_df <- data.frame(
  target = rep(seq_len(n_sp), times=n_sp),
  source = rep(seq_len(n_sp), each=n_sp),
  value  = as.vector(Aheat_ord)
)
heat_df <- heat_df[heat_df$value != 0, ]

p_heat <- ggplot(heat_df, aes(x=source, y=target, fill=value)) +
  geom_tile() +
  scale_fill_gradient2(low="steelblue3", mid="white", high="tomato3",
                       midpoint=0, name=expression(hat(A)[ji])) +
  scale_x_continuous(breaks=ticks, labels=otu_names[ord_nodes][ticks], expand=c(0,0)) +
  scale_y_continuous(breaks=ticks, labels=otu_names[ord_nodes][ticks], expand=c(0,0)) +
  labs(title="PSS-Net v0: Inferred interaction matrix",
       subtitle="Rows/cols ordered by in-degree (high → low)  |  blank = no edge",
       x="Source OTU", y="Target OTU") +
  theme_bw(base_size=8) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=5.5),
        axis.text.y=element_text(size=5.5),
        plot.title=element_text(size=10, face="bold"),
        plot.subtitle=element_text(size=8, color="grey45"),
        legend.key.height=unit(0.9,"cm"))
print(p_heat)

# ── Plot 4: igraph network (force-directed, all n_sp nodes) ──────────────────
if (n_inf > 0) {
  # Build graph with all n_sp nodes so vertex indices match OTU indices exactly
  g_net <- make_empty_graph(n_sp, directed=TRUE)
  g_net <- add_edges(g_net, as.vector(t(cbind(el[,2], el[,1]))))
  V(g_net)$name   <- otu_names
  E(g_net)$weight <- A_hat[el]
  
  vsize  <- 3 + 8*(in_deg - min(in_deg)) / (max(in_deg) - min(in_deg) + 1e-9)
  vlabel <- otu_names
  ecol   <- ifelse(E(g_net)$weight > 0, "tomato3", "steelblue3")
  
  set.seed(42)
  lay <- layout_with_fr(g_net, niter = 800, weights = NA)
  
  par(mar=c(1,1,3.5,1), bg="white")
  plot(g_net,
       layout=lay,
       vertex.color="grey88", vertex.frame.color="grey60",
       vertex.size=vsize, vertex.label=vlabel,
       vertex.label.cex=0.40, vertex.label.color="grey10",
       edge.color=ecol, edge.width=0.9,
       edge.arrow.size=0.25, edge.curved=0.2,
       main=sprintf(
         "PSS-Net v0 — Inferred network\n%d OTUs | %d edges | N=%d conditions | M=%d",
         n_sp, n_inf, N_ok, M))
  legend("bottomleft", bty="n", cex=0.8,
         legend=c("Promotion (+)", "Inhibition (-)"),
         col=c("tomato3","steelblue3"), lwd=2)
  par(mar=c(5,4,4,2))
}

# ── Plot 5: μ̂ vs γ̂ scatter per OTU ────────────────────────────────────────
param_df <- data.frame(
  mu  = MU_ads,
  gam = g_hat,
  otu = otu_names,
  deg = in_deg
)
p_par <- ggplot(param_df, aes(x=mu, y=gam, color=deg)) +
  geom_point(size=2, alpha=0.8) +
  scale_color_gradient(low="grey80", high="tomato3", name="In-degree") +
  labs(title="PSS-Net v0: Inferred μ̂ vs γ̂ per OTU",
       subtitle="Color = in-degree (darker = more regulators)",
       x=expression(hat(mu)[j]~"(intrinsic growth)"),
       y=expression(hat(gamma)[j]~"(self-regulation)")) +
  theme_bw(base_size=10) +
  theme(plot.title=element_text(size=11, face="bold"),
        plot.subtitle=element_text(size=9, color="grey45"))
print(p_par)

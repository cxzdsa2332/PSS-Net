################################################################################
# pss_net_compare.R  —  PSS-Net: 10-run MCC benchmark: v3 (no smooth) vs v1 (nonlinear)
#
# Input:   none (generates its own simulated data per seed)
# Output:  results/mcc_comparison.csv   — per-run metrics
#          results/mcc_comparison.txt   — summary table
################################################################################

rm(list = ls())

library(deSolve)
library(ADSIHT)
library(grpreg)
library(splines)

# ── Fixed true parameters (shared across all runs) ────────────────────────────
n_sp     <- 8
sp_names <- paste0("Sp", seq_len(n_sp))
M        <- 2

r_true   <- c(0.8, 1.2, 0.6, 1.0, 0.7, 1.1, 0.5, 0.9)
gam_true <- c(1.5, 1.8, 1.2, 1.6, 1.4, 1.5, 1.0, 1.6)

A_true <- matrix(0, n_sp, n_sp)
A_true[1,3] <-  0.40; A_true[1,5] <- -0.30
A_true[2,1] <-  0.30; A_true[2,4] <- -0.40; A_true[2,7] <-  0.20
A_true[3,2] <- -0.30; A_true[3,6] <-  0.30
A_true[4,1] <-  0.40; A_true[4,3] <- -0.20; A_true[4,8] <-  0.30
A_true[5,2] <-  0.30; A_true[5,6] <- -0.20
A_true[6,4] <-  0.20; A_true[6,5] <-  0.30
A_true[7,3] <- -0.30; A_true[7,8] <-  0.40
A_true[8,1] <-  0.20; A_true[8,6] <- -0.20

adj_true <- (A_true != 0) * 1L
n_edges  <- sum(A_true != 0)

alpha_true <- matrix(0, n_sp, n_sp * M)
for (j in seq_len(n_sp))
  for (i in seq_len(n_sp))
    alpha_true[j, (i - 1) * M + 1] <- A_true[j, i]

ode_func <- function(t, state, parms) {
  x <- pmax(state, 0)
  list(r_true + as.numeric(A_true %*% x) - gam_true * x + parms$u)
}

# ── Core inference function ───────────────────────────────────────────────────
# smooth = FALSE → v3 (no pre-smoothing)
# smooth = TRUE  → v1 (B-spline pre-smoothing)
run_pss_net <- function(seed, smooth = FALSE) {
  set.seed(seed)

  # ── Simulate PSS data ──────────────────────────────────────────────────────
  N_cond <- 300
  U_mat  <- matrix(runif(N_cond * n_sp, -0.4, 0.8), N_cond, n_sp)
  U_mat[1, ] <- 0

  ss_mat <- matrix(NA, N_cond, n_sp)
  for (k in seq_len(N_cond)) {
    tryCatch({
      out <- ode(runif(n_sp, 0.5, 2.0), c(0, 5000), ode_func,
                 list(u = U_mat[k, ]), method = "lsoda",
                 rtol = 1e-9, atol = 1e-11)
      ss_mat[k, ] <- out[nrow(out), 2:(n_sp + 1)]
    }, error = function(e) NULL)
  }

  ok    <- apply(ss_mat, 1, function(r) all(is.finite(r) & r > 0))
  X_obs <- ss_mat[ok, ] + matrix(rnorm(sum(ok)*n_sp, 0, 0.03), sum(ok), n_sp)
  X_obs <- pmax(X_obs, 1e-6)
  U_obs <- U_mat[ok, ]
  N     <- nrow(X_obs)

  # Noise-free WT steady state
  out_wt <- ode(rep(1, n_sp), c(0, 1e4), ode_func, list(u = rep(0, n_sp)),
                method = "lsoda", rtol = 1e-12, atol = 1e-14)
  x_wt   <- as.numeric(out_wt[nrow(out_wt), 2:(n_sp + 1)])

  # ── Pre-smoothing (v1 only) ────────────────────────────────────────────────
  X_basis <- X_obs
  if (smooth) {
    for (j in seq_len(n_sp)) {
      fit_sm        <- lm(X_obs[, j] ~ bs(U_obs[, j], df = 6))
      X_basis[, j]  <- pmax(fitted(fit_sm), 1e-6)
    }
  }

  # ── Polynomial basis ───────────────────────────────────────────────────────
  Psi <- matrix(0, N, n_sp * M)
  for (i in seq_len(n_sp))
    for (m in seq_len(M))
      Psi[, (i - 1) * M + m] <- X_basis[, i]^m

  group   <- rep(seq_len(n_sp), each = M)
  Psi_bar <- colMeans(Psi)
  Psi_c   <- sweep(Psi, 2, Psi_bar)
  Psi_sd  <- pmax(apply(Psi_c, 2, sd), 1e-10)
  Psi_cs  <- sweep(Psi_c, 2, Psi_sd, "/")
  U_bar   <- colMeans(U_obs)
  U_c     <- sweep(U_obs, 2, U_bar)

  # ── Sparse regression ─────────────────────────────────────────────────────
  ALPHA_ads <- matrix(0, n_sp, n_sp * M)
  ALPHA_gl  <- matrix(0, n_sp, n_sp * M)
  MU_ads    <- numeric(n_sp)
  MU_gl     <- numeric(n_sp)

  for (j in seq_len(n_sp)) {
    rhs_c <- -U_c[, j]

    fit_a         <- ADSIHT(Psi_cs, matrix(rhs_c), group, ic.type = "dsic")
    best_a        <- which.min(fit_a$ic)
    ALPHA_ads[j,] <- fit_a$beta[, best_a] / Psi_sd
    MU_ads[j]     <- -U_bar[j] - sum(Psi_bar * ALPHA_ads[j, ])

    cv            <- cv.grpreg(Psi_cs, rhs_c, group = group,
                               penalty = "grLasso", nfolds = 5)
    ALPHA_gl[j, ] <- coef(cv)[-1] / Psi_sd
    MU_gl[j]      <- -U_bar[j] - sum(Psi_bar * ALPHA_gl[j, ])
  }

  # ── Jacobian ───────────────────────────────────────────────────────────────
  get_jacobian <- function(alpha_mat) {
    j_mat <- matrix(0, n_sp, n_sp)
    for (j in seq_len(n_sp))
      for (i in seq_len(n_sp)) {
        cols        <- (i - 1) * M + seq_len(M)
        dpsi_i      <- seq_len(M) * x_wt[i]^pmax(seq_len(M) - 1, 0)
        j_mat[j, i] <- sum(alpha_mat[j, cols] * dpsi_i)
      }
    j_mat
  }

  get_A   <- function(J) { A <- J; diag(A) <- 0; A }
  get_gam <- function(J) -diag(J)

  J_ads <- get_jacobian(ALPHA_ads)
  J_gl  <- get_jacobian(ALPHA_gl)

  # ── Edge detection ─────────────────────────────────────────────────────────
  group_norm_mat <- function(ALPHA) {
    mat <- matrix(0, n_sp, n_sp)
    for (j in seq_len(n_sp))
      for (i in seq_len(n_sp)) {
        cols      <- (i - 1) * M + seq_len(M)
        mat[j, i] <- sqrt(sum(ALPHA[j, cols]^2))
      }
    mat
  }

  norm_ads   <- group_norm_mat(ALPHA_ads)
  norm_gl    <- group_norm_mat(ALPHA_gl)
  tau_gl_vec <- apply(norm_gl, 1, max) * 0.01

  adj_ads <- { adj <- norm_ads >= 1e-10; diag(adj) <- FALSE; adj * 1L }
  adj_gl  <- {
    adj <- matrix(FALSE, n_sp, n_sp)
    for (j in seq_len(n_sp)) adj[j, ] <- norm_gl[j, ] >= tau_gl_vec[j]
    diag(adj) <- FALSE; adj * 1L
  }

  # ── Metrics ────────────────────────────────────────────────────────────────
  calc_mets <- function(adj_h, ALPHA_h) {
    TP <- sum(adj_true == 1 & adj_h == 1)
    FP <- sum(adj_true == 0 & adj_h == 1)
    FN <- sum(adj_true == 1 & adj_h == 0)
    TN <- sum(adj_true == 0 & adj_h == 0)
    pr    <- TP / (TP + FP + 1e-9)
    re    <- TP / (TP + FN + 1e-9)
    f1    <- 2 * pr * re / (pr + re + 1e-9)
    denom <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    mcc   <- if (denom < 1e-9) 0 else (TP * TN - FP * FN) / denom
    coef_l2  <- mean(sqrt(rowSums((ALPHA_h - alpha_true)^2)))
    jac_rmse <- sqrt(mean((get_A(get_jacobian(ALPHA_h)) - A_true)^2))
    c(TP = TP, FP = FP, FN = FN, TN = TN,
      Pr = pr, Re = re, F1 = f1, MCC = mcc,
      CoefL2 = coef_l2, JacRMSE = jac_rmse)
  }

  list(
    seed = seed, N = N,
    ads = calc_mets(adj_ads, ALPHA_ads),
    gl  = calc_mets(adj_gl,  ALPHA_gl)
  )
}

# ── Run 10 replicates ─────────────────────────────────────────────────────────
seeds <- 1:10
n_runs <- length(seeds)

cat(sprintf("Running %d replicates × 2 versions × 2 methods...\n\n", n_runs))

results_v3 <- vector("list", n_runs)
results_v1 <- vector("list", n_runs)

for (s in seq_along(seeds)) {
  cat(sprintf("[Seed %2d] ", seeds[s]))
  cat("v3 (no smooth)... ")
  results_v3[[s]] <- run_pss_net(seeds[s], smooth = FALSE)
  cat(sprintf("ADSIHT MCC=%.3f  grLasso MCC=%.3f  | ",
              results_v3[[s]]$ads["MCC"], results_v3[[s]]$gl["MCC"]))
  cat("v1 (pre-smooth)... ")
  results_v1[[s]] <- run_pss_net(seeds[s], smooth = TRUE)
  cat(sprintf("ADSIHT MCC=%.3f  grLasso MCC=%.3f\n",
              results_v1[[s]]$ads["MCC"], results_v1[[s]]$gl["MCC"]))
}

# ── Build results data frame ──────────────────────────────────────────────────
build_df <- function(res_list, version) {
  do.call(rbind, lapply(res_list, function(r) {
    rbind(
      data.frame(version = version, method = "ADSIHT",
                 seed = r$seed, N = r$N,
                 t(r$ads), stringsAsFactors = FALSE),
      data.frame(version = version, method = "grLasso",
                 seed = r$seed, N = r$N,
                 t(r$gl),  stringsAsFactors = FALSE)
    )
  }))
}

df_v3 <- build_df(results_v3, "v3_no_smooth")
df_v1 <- build_df(results_v1, "v1_pre_smooth")
df_all <- rbind(df_v3, df_v1)
rownames(df_all) <- NULL

# ── Save CSV ──────────────────────────────────────────────────────────────────
dir.create("results", showWarnings = FALSE)
write.csv(df_all, "results/mcc_comparison.csv", row.names = FALSE)
cat("\nSaved: results/mcc_comparison.csv\n")

# ── Summary statistics ────────────────────────────────────────────────────────
summary_tab <- function(df, ver, meth) {
  sub <- df[df$version == ver & df$method == meth, ]
  c(version = ver, method = meth,
    MCC_mean  = round(mean(sub$MCC),  3),
    MCC_sd    = round(sd(sub$MCC),    3),
    MCC_min   = round(min(sub$MCC),   3),
    MCC_max   = round(max(sub$MCC),   3),
    F1_mean   = round(mean(sub$F1),   3),
    Pr_mean   = round(mean(sub$Pr),   3),
    Re_mean   = round(mean(sub$Re),   3),
    TP_mean   = round(mean(sub$TP),   2),
    FP_mean   = round(mean(sub$FP),   2),
    FN_mean   = round(mean(sub$FN),   2),
    CoefL2_mean  = round(mean(sub$CoefL2),   4),
    JacRMSE_mean = round(mean(sub$JacRMSE),  4))
}

combos <- list(
  c("v3_no_smooth",  "ADSIHT"),
  c("v3_no_smooth",  "grLasso"),
  c("v1_pre_smooth", "ADSIHT"),
  c("v1_pre_smooth", "grLasso")
)

sum_df <- do.call(rbind, lapply(combos, function(x)
  as.data.frame(t(summary_tab(df_all, x[1], x[2])),
                stringsAsFactors = FALSE)))

# ── Print summary ─────────────────────────────────────────────────────────────
cat("\n")
cat("================================================================\n")
cat(sprintf("  PSS-Net: v3 (no smooth) vs v1 (pre-smooth)  |  %d runs\n", n_runs))
cat(sprintf("  %d species | %d true edges | M=%d | N_cond=300 | σ=0.03\n",
            n_sp, n_edges, M))
cat("================================================================\n")
cat(sprintf("  %-18s %-8s  %6s±%-5s  [%5s,%5s]  F1=%5s  Pr=%5s  Re=%5s  TP=%4s FP=%4s FN=%4s\n",
            "Version", "Method",
            "MCC", "sd", "min", "max", "mean", "mean", "mean", "avg", "avg", "avg"))
cat("  ", strrep("-", 100), "\n", sep = "")

for (i in seq_len(nrow(sum_df))) {
  r <- sum_df[i, ]
  # separator between versions
  if (i == 3) cat("  ", strrep("·", 100), "\n", sep = "")
  cat(sprintf("  %-18s %-8s  %6s±%-5s  [%5s,%5s]  F1=%5s  Pr=%5s  Re=%5s  TP=%4s FP=%4s FN=%4s\n",
              r$version, r$method,
              r$MCC_mean, r$MCC_sd,
              r$MCC_min, r$MCC_max,
              r$F1_mean, r$Pr_mean, r$Re_mean,
              r$TP_mean, r$FP_mean, r$FN_mean))
}
cat("================================================================\n")

# ── Per-run MCC table ─────────────────────────────────────────────────────────
cat("\n── Per-run MCC ────────────────────────────────────────────────\n")
cat(sprintf("  %5s  %10s  %10s  %10s  %10s  | Δ_ADS  Δ_GL\n",
            "Seed",
            "v3-ADSIHT", "v1-ADSIHT",
            "v3-grLasso", "v1-grLasso"))
cat("  ", strrep("-", 72), "\n", sep = "")

for (s in seq_along(seeds)) {
  mcc_v3_a <- results_v3[[s]]$ads["MCC"]
  mcc_v1_a <- results_v1[[s]]$ads["MCC"]
  mcc_v3_g <- results_v3[[s]]$gl["MCC"]
  mcc_v1_g <- results_v1[[s]]$gl["MCC"]
  cat(sprintf("  %5d  %10.3f  %10.3f  %10.3f  %10.3f  | %+.3f  %+.3f\n",
              seeds[s],
              mcc_v3_a, mcc_v1_a,
              mcc_v3_g, mcc_v1_g,
              mcc_v1_a - mcc_v3_a,
              mcc_v1_g - mcc_v3_g))
}
cat("  ", strrep("-", 72), "\n", sep = "")
cat(sprintf("  %5s  %10.3f  %10.3f  %10.3f  %10.3f  | %+.3f  %+.3f\n",
            "Mean",
            mean(sapply(results_v3, function(r) r$ads["MCC"])),
            mean(sapply(results_v1, function(r) r$ads["MCC"])),
            mean(sapply(results_v3, function(r) r$gl["MCC"])),
            mean(sapply(results_v1, function(r) r$gl["MCC"])),
            mean(sapply(results_v1, function(r) r$ads["MCC"])) -
              mean(sapply(results_v3, function(r) r$ads["MCC"])),
            mean(sapply(results_v1, function(r) r$gl["MCC"])) -
              mean(sapply(results_v3, function(r) r$gl["MCC"]))))
cat("\n")

# ── Save text summary ─────────────────────────────────────────────────────────
sink("results/mcc_comparison.txt")

cat("================================================================\n")
cat(sprintf("  PSS-Net: v3 (no smooth) vs v1 (pre-smooth)  |  %d runs\n", n_runs))
cat(sprintf("  %d species | %d true edges | M=%d | N_cond=300 | σ=0.03\n",
            n_sp, n_edges, M))
cat(sprintf("  Pre-smoothing: lm(x_j ~ bs(u_j, df=6)) per species\n"))
cat("================================================================\n")
cat(sprintf("  %-18s %-8s  %6s±%-5s  [%5s,%5s]  F1=%5s  Pr=%5s  Re=%5s  TP=%4s FP=%4s FN=%4s\n",
            "Version", "Method",
            "MCC", "sd", "min", "max", "mean", "mean", "mean", "avg", "avg", "avg"))
cat("  ", strrep("-", 100), "\n", sep = "")

for (i in seq_len(nrow(sum_df))) {
  r <- sum_df[i, ]
  if (i == 3) cat("  ", strrep("·", 100), "\n", sep = "")
  cat(sprintf("  %-18s %-8s  %6s±%-5s  [%5s,%5s]  F1=%5s  Pr=%5s  Re=%5s  TP=%4s FP=%4s FN=%4s\n",
              r$version, r$method,
              r$MCC_mean, r$MCC_sd,
              r$MCC_min, r$MCC_max,
              r$F1_mean, r$Pr_mean, r$Re_mean,
              r$TP_mean, r$FP_mean, r$FN_mean))
}
cat("================================================================\n\n")

cat("── Per-run MCC ────────────────────────────────────────────────\n")
cat(sprintf("  %5s  %10s  %10s  %10s  %10s  | Δ_ADS  Δ_GL\n",
            "Seed", "v3-ADSIHT", "v1-ADSIHT", "v3-grLasso", "v1-grLasso"))
cat("  ", strrep("-", 72), "\n", sep = "")

for (s in seq_along(seeds)) {
  mcc_v3_a <- results_v3[[s]]$ads["MCC"]
  mcc_v1_a <- results_v1[[s]]$ads["MCC"]
  mcc_v3_g <- results_v3[[s]]$gl["MCC"]
  mcc_v1_g <- results_v1[[s]]$gl["MCC"]
  cat(sprintf("  %5d  %10.3f  %10.3f  %10.3f  %10.3f  | %+.3f  %+.3f\n",
              seeds[s], mcc_v3_a, mcc_v1_a, mcc_v3_g, mcc_v1_g,
              mcc_v1_a - mcc_v3_a, mcc_v1_g - mcc_v3_g))
}
cat("  ", strrep("-", 72), "\n", sep = "")
cat(sprintf("  %5s  %10.3f  %10.3f  %10.3f  %10.3f  | %+.3f  %+.3f\n",
            "Mean",
            mean(sapply(results_v3, function(r) r$ads["MCC"])),
            mean(sapply(results_v1, function(r) r$ads["MCC"])),
            mean(sapply(results_v3, function(r) r$gl["MCC"])),
            mean(sapply(results_v1, function(r) r$gl["MCC"])),
            mean(sapply(results_v1, function(r) r$ads["MCC"])) -
              mean(sapply(results_v3, function(r) r$ads["MCC"])),
            mean(sapply(results_v1, function(r) r$gl["MCC"])) -
              mean(sapply(results_v3, function(r) r$gl["MCC"]))))

sink()
cat("Saved: results/mcc_comparison.txt\n")

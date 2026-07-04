#!/usr/bin/env Rscript
# ============================================================================
# P2 (scoped): a joint influence-function (delta-method/sandwich) SE for the
# D-corrected win-MR estimand tau_hat = gamma_hat / D_hat, in the ONE-SAMPLE
# AR setting (mrwin_ar_onesample), where the needed per-subject pieces are
# already native to the package's .mrwin_ar1_blocks().
#
# MOTIVATION (from probe #2, p2-tau-portable-winmr-twosample.R): naively
# rescaling gamma_hat's CI by a FIXED D_hat (tau_ci <- ar$ci / D_hat) improves
# coverage of the portable target tau* over the raw CI at heavy censoring
# (0.875 vs 0.662 at cens=5.0) but is not fully calibrated -- because it
# ignores D_hat's OWN sampling variance and its covariance with gamma_hat
# (both are computed from the SAME sample). This is exactly a job for a
# proper multivariate delta method / EIF-style combination.
#
# DERIVATION. mrwin_ar_onesample's AR estimator solves (approximately, via the
# Hajek/first-order projection) the estimating equation
#   E[ gh_i(beta) - beta*... ] ... more precisely: gh_i(beta) := gh_i - beta*gx_i
#   (gh_i, gx_i as in .mrwin_ar1_blocks; c0=Var(gh), c1=Cov(gh,gx), c2=Var(gx),
#   a=mean(gh_i(0))-ish population moment, b = the "denominator" slope), with
#   gamma_hat = a/b solving mean(gh_i - gamma*gx_i) = 0 in the linearized sense.
# By M-estimation: IF_gamma(o) = (gh(o;gamma_hat)) / b_hat   [standard sandwich:
#   -1/(d/dbeta E[gh_i(beta)]) * gh_i(gamma_hat), and d/dbeta E[gh_i(beta)]=-b]
#
# D_hat = mean pairwise "decided" rate is a degree-2 U-statistic; its Hajek
# projection gives IF_D(o) = 2*(d(o) - D_hat), d(o) = per-subject decided rate
# (already computed as an intermediate inside .mrwin_ar1_blocks, here recomputed
# directly via the same subj() kernel).
#
# Delta method for tau=gamma/D: IF_tau(o) = IF_gamma(o)/D - gamma*IF_D(o)/D^2
#   Var_hat(tau_hat) = var(IF_tau_i) / n   (empirical sandwich variance, no bootstrap)
#
# TEST: across replications, does mean(SE_hat) calibrate against the empirical
# SD of tau_hat (ratio ~ 1, CLAUDE.md's calibration check), and does the Wald
# CI tau_hat +/- 1.96*SE_hat achieve closer-to-nominal coverage of tau* than
# the naive ar$ci/D_hat plug-in, especially at heavy censoring? Public sim +
# package (mrwin_config/mrwin_simulate-style one-sample DGP) only.
# ============================================================================
suppressMessages(library(mrwin))
ncores <- max(1L, min(4L, parallel::detectCores()))
subj <- getFromNamespace("mrwin_subject_win_loss_cpp", "mrwin")
blocks <- getFromNamespace(".mrwin_ar1_blocks", "mrwin")

cfg0 <- mrwin_config(m_snps = 1L, n_outcome = 1200L)
set.seed(3); maf <- runif(1, cfg0$maf_low, cfg0$maf_high); beta1 <- rnorm(1, 0, cfg0$sigma_beta)

gen_one <- function(n, seed, alpha_s, cens) {
  set.seed(seed)
  g <- scale(rbinom(n, 2, maf)); storage.mode(g) <- "double"
  Z <- as.numeric(g) * beta1
  u <- rnorm(n); w <- rgamma(n, 1 / cfg0$theta_f, scale = cfg0$theta_f)
  cn <- pmin(rexp(n, cens), cfg0$max_follow_up); uf <- matrix(runif(n * 3), n, 3)
  X <- alpha_s * Z + cfg0$alpha_u * u + rnorm(n)
  et <- matrix(NA_real_, n, 3); lw <- log(w)
  for (j in 1:3) { lp <- cfg0$alpha_x[j] * X + cfg0$nu_u[j] * u + lw
    et[, j] <- (-log(uf[, j]) / (cfg0$baseline_haz[j] * exp(lp)))^(1 / cfg0$shape_weibull) }
  td <- et[, 1]; th <- et[, 2]; tr <- et[, 3]
  list(time = cbind(pmin(td, cn), pmin(th, cn, td), pmin(tr, cn, td)),
       status = cbind(as.integer(td <= cn), as.integer(th <= cn & th <= td), as.integer(tr <= cn & tr <= td)),
       X = X, Z = Z)
}
one_rep <- function(n, seed, alpha_s, cens) {
  d <- gen_one(n, seed, alpha_s, cens)
  N <- length(d$X)
  m <- blocks(d$time, d$status, d$X, d$Z)
  gamma_hat <- m$a / m$b

  # per-subject gh_i, gx_i (recomputed directly; gh_i already appears as `m$s_i`-derived inside blocks()
  # but the per-subject VECTOR isn't returned by the exported list, so recompute identically here)
  cw <- subj(d$time, d$status, d$time, d$status, rep(1, N)); s_i <- cw[, 1] - cw[, 2]
  cz <- subj(d$time, d$status, d$time, d$status, d$Z); r_i <- cz[, 1] - cz[, 2]
  gh_i <- (d$Z * s_i - r_i) / (N - 1)
  SX <- sum(d$X); SZ <- sum(d$Z); SXZ <- sum(d$X * d$Z)
  gx_i <- ((N - 1) * d$X * d$Z - d$X * (SZ - d$Z) - d$Z * (SX - d$X) + (SXZ - d$X * d$Z)) / (N - 1)
  d_i <- (cw[, 1] + cw[, 2]) / (N - 1)
  D_hat <- mean(d_i)

  IF_gamma <- (gh_i - gamma_hat * gx_i) / m$b
  IF_D <- 2 * (d_i - D_hat)
  tau_hat <- gamma_hat / D_hat
  IF_tau <- IF_gamma / D_hat - gamma_hat * IF_D / D_hat^2
  se_tau <- sqrt(stats::var(IF_tau) / N)

  # naive comparator: AR's own quadratic-inversion CI, plug-in rescaled by (fixed) D_hat
  ar <- mrwin_ar_onesample(d$time, d$status, d$X, d$Z, degeneracy_check = FALSE)
  naive_ci <- ar$ci / D_hat

  c(gamma_hat = gamma_hat, D_hat = D_hat, tau_hat = tau_hat, se_tau = se_tau,
    naive_lo = naive_ci[1], naive_hi = naive_ci[2],
    delta_lo = tau_hat - 1.96 * se_tau, delta_hi = tau_hat + 1.96 * se_tau)
}

cat("== P2: delta-method SE for tau_hat=gamma_hat/D_hat vs naive ar$ci/D_hat, one-sample AR ==\n")
N <- 3000L; R <- 150L; alpha_s <- 0.15
cens_grid <- c(0.02, 0.2, 0.8, 2.5, 6.0)
cat(sprintf("N=%d, R=%d, alpha_s=%.2f (one-sample AR, mrwin_ar_onesample)\n\n", N, R, alpha_s))
cat(sprintf("%-8s %-12s %-14s %-14s %-16s %-16s\n",
            "cens", "mean(tau_hat)", "empSD(tau_hat)", "mean(SE_delta)", "calib.ratio", "D_hat(mean)"))
allres <- list()
for (cens in cens_grid) {
  rr <- do.call(rbind, parallel::mclapply(seq_len(R), function(r) {
    tryCatch(one_rep(N, 50000L + round(cens * 1000) + r, alpha_s, cens), error = function(e) NULL)
  }, mc.cores = ncores))
  rr <- rr[stats::complete.cases(rr), , drop = FALSE]
  allres[[as.character(cens)]] <- rr
  emp_sd <- sd(rr[, "tau_hat"]); mean_se <- mean(rr[, "se_tau"])
  cat(sprintf("%-8.2f %-12.4f %-14.4f %-14.4f %-16.3f %-16.4f\n",
              cens, mean(rr[, "tau_hat"]), emp_sd, mean_se, mean_se / emp_sd, mean(rr[, "D_hat"])))
}

cat("\n== Reference tau* (large N, near-uncensored one-sample DGP) ==\n")
ref <- one_rep(60000L, 999L, alpha_s, cens_grid[1])
tau_star <- ref["tau_hat"]
cat(sprintf("tau* = %.4f (large-N gamma_hat=%.4f / D_hat=%.4f)\n\n", tau_star, ref["gamma_hat"], ref["D_hat"]))

cat("== Coverage of tau* : naive (ar$ci/D_hat) vs delta-method Wald CI ==\n")
cat(sprintf("%-8s %-16s %-16s\n", "cens", "naive-cov(tau*)", "delta-cov(tau*)"))
for (cens in cens_grid) {
  rr <- allres[[as.character(cens)]]
  naive_cov <- mean(tau_star >= rr[, "naive_lo"] & tau_star <= rr[, "naive_hi"])
  delta_cov <- mean(tau_star >= rr[, "delta_lo"] & tau_star <= rr[, "delta_hi"])
  cat(sprintf("%-8.2f %-16.3f %-16.3f\n", cens, naive_cov, delta_cov))
}
cat("DONE-P2-DELTA\n")

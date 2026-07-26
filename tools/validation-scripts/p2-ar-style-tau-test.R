#!/usr/bin/env Rscript
# ============================================================================
# P2, the "correct fix" flagged as open in dev/findings-transport-collapsibility.md
# §11: both a naive plug-in CI (ar$ci/D_hat, treats D_hat as fixed) and a
# delta-method Wald CI on tau_hat=gamma_hat/D_hat (treats tau_hat as a ratio
# point estimate) were shown to be inadequate under weak instruments -- the
# delta-method one badly so (SE miscalibrated 2-27x), because BOTH approaches
# eventually divide by a possibly-small/noisy quantity, which is exactly the
# irregularity Anderson-Rubin-type tests were invented to avoid for gamma_hat
# itself. This script builds the analogous AR-type test DIRECTLY for tau,
# never forming tau_hat as a free-standing ratio.
#
# DERIVATION. tau solves the population moment  M(tau) := a - b*tau*D = 0
# (since gamma=a/b=tau*D by definition of tau=gamma/D). Its per-subject
# (Hajek/linearized) representation, treating tau as FIXED (not estimated):
#   IF_M_i(tau) = gh_i - tau*(D_hat*gx_i + b_hat*IF_D_i),   IF_D_i = 2*(d_i-D_hat)
# (product rule for tau*b*D: d(bD) = D*db + b*dD). Writing u_i := D_hat*gx_i +
# b_hat*IF_D_i (fixed per-subject vector, no tau dependence), this is LINEAR in
# tau: IF_M_i(tau) = gh_i - tau*u_i, so its variance is an EXACT quadratic form
#   Var(IF_M(tau)) = c0' - 2*tau*c1' + tau^2*c2',  c0'=Var(gh), c1'=Cov(gh,u),
#   c2'=Var(u)   [same shape as .mrwin_ar1_blocks's c0,c1,c2 for gamma, with
#   gx replaced by u -- so the whole quadratic-inversion machinery of
#   mrwin_ar_onesample carries over exactly, just relabeling b->b*D, gx->u].
#
# AR_tau(tau) := n*(a - b*D*tau)^2 / (4*(c0'-2*tau*c1'+tau^2*c2'))
# Confidence set: {tau : AR_tau(tau) <= qchisq(1-alpha,1)}, solved by the SAME
# quadratic A*tau^2+B*tau+C<=0 closed form as mrwin_ar_onesample (A,B,C below).
#
# KEY DIFFERENCE from both prior attempts: this never treats D_hat as a fixed
# constant (its own sampling noise enters c1',c2' through u_i) AND never forms
# tau_hat as a free ratio requiring delta-method regularity -- it is a direct
# generalization of the AR philosophy (test candidate values via a moment,
# invert) from gamma to tau.
#
# TEST: coverage of tau* across the SAME censoring ladder used in
# p2-delta-method-tau-se.R, compared against both prior (inadequate) methods.
# One-sample AR DGP; public sim + package only.
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

# AR-style quadratic-inversion CI for tau, at level alpha (mirrors
# mrwin_ar_onesample's A,B,C construction exactly, with b->b*D, gx->u)
ar_tau_ci <- function(d, alpha = 0.05) {
  N <- length(d$X)
  m <- blocks(d$time, d$status, d$X, d$Z)
  cw <- subj(d$time, d$status, d$time, d$status, rep(1, N)); s_i <- cw[, 1] - cw[, 2]
  cz <- subj(d$time, d$status, d$time, d$status, d$Z); r_i <- cz[, 1] - cz[, 2]
  gh_i <- (d$Z * s_i - r_i) / (N - 1)
  SX <- sum(d$X); SZ <- sum(d$Z); SXZ <- sum(d$X * d$Z)
  gx_i <- ((N - 1) * d$X * d$Z - d$X * (SZ - d$Z) - d$Z * (SX - d$X) + (SXZ - d$X * d$Z)) / (N - 1)
  d_i <- (cw[, 1] + cw[, 2]) / (N - 1); D_hat <- mean(d_i)
  IF_D_i <- 2 * (d_i - D_hat)

  bD <- m$b * D_hat
  u_i <- D_hat * gx_i + m$b * IF_D_i
  c0p <- stats::var(gh_i); c1p <- stats::cov(gh_i, u_i); c2p <- stats::var(u_i)

  q <- stats::qchisq(1 - alpha, 1)
  A <- N * bD^2 - 4 * q * c2p
  B <- -2 * N * m$a * bD + 8 * q * c1p
  C <- N * m$a^2 - 4 * q * c0p
  disc <- B^2 - 4 * A * C
  tau_point <- m$a / bD  # argmin (still a ratio numerically, but the CI below does not rely on its regularity)
  bounded <- FALSE; empty <- FALSE; ci <- c(-Inf, Inf)
  if (abs(A) < .Machine$double.eps^0.5) {
    if (abs(B) < .Machine$double.eps^0.5) { empty <- C > 0 } else {
      root <- -C / B; ci <- if (B > 0) c(-Inf, root) else c(root, Inf)
    }
  } else if (disc < 0) {
    if (A > 0) empty <- TRUE else ci <- c(-Inf, Inf)
  } else {
    r1 <- (-B - sqrt(disc)) / (2 * A); r2 <- (-B + sqrt(disc)) / (2 * A)
    lo <- min(r1, r2); hi <- max(r1, r2)
    if (A > 0) { ci <- c(lo, hi); bounded <- TRUE } else { ci <- c(-Inf, Inf) }
  }
  list(tau_point = tau_point, ci = ci, bounded = bounded, empty = empty, D_hat = D_hat, gamma_hat = m$a / m$b)
}

# comparators from the earlier (inadequate) attempts, for a head-to-head
naive_and_delta <- function(d) {
  N <- length(d$X)
  m <- blocks(d$time, d$status, d$X, d$Z)
  cw <- subj(d$time, d$status, d$time, d$status, rep(1, N)); s_i <- cw[, 1] - cw[, 2]
  cz <- subj(d$time, d$status, d$time, d$status, d$Z); r_i <- cz[, 1] - cz[, 2]
  gh_i <- (d$Z * s_i - r_i) / (N - 1)
  SX <- sum(d$X); SZ <- sum(d$Z); SXZ <- sum(d$X * d$Z)
  gx_i <- ((N - 1) * d$X * d$Z - d$X * (SZ - d$Z) - d$Z * (SX - d$X) + (SXZ - d$X * d$Z)) / (N - 1)
  d_i <- (cw[, 1] + cw[, 2]) / (N - 1); D_hat <- mean(d_i)
  gamma_hat <- m$a / m$b; tau_hat <- gamma_hat / D_hat
  IF_gamma <- (gh_i - gamma_hat * gx_i) / m$b
  IF_D <- 2 * (d_i - D_hat)
  IF_tau <- IF_gamma / D_hat - gamma_hat * IF_D / D_hat^2
  se_tau <- sqrt(stats::var(IF_tau) / N)
  ar <- mrwin_ar_onesample(d$time, d$status, d$X, d$Z, degeneracy_check = FALSE)
  list(naive_ci = ar$ci / D_hat, delta_ci = c(tau_hat - 1.96 * se_tau, tau_hat + 1.96 * se_tau))
}

cat("== P2 fix: AR-style quadratic-inversion test DIRECTLY for tau vs. naive plug-in and delta-method ==\n")
N <- 3000L; R <- 150L; alpha_s <- 0.15
cens_grid <- c(0.02, 0.2, 0.8, 2.5, 6.0)
cat(sprintf("N=%d, R=%d, alpha_s=%.2f (one-sample AR)\n\n", N, R, alpha_s))

ref <- ar_tau_ci(gen_one(60000L, 999L, alpha_s, cens_grid[1]))
tau_star <- ref$tau_point
cat(sprintf("tau* (large-N reference) = %.4f (gamma_hat=%.4f, D_hat=%.4f)\n\n", tau_star, ref$gamma_hat, ref$D_hat))

cat(sprintf("%-8s %-14s %-16s %-16s %-16s\n", "cens", "AR-tau bounded%", "AR-tau cov(tau*)", "naive cov(tau*)", "delta cov(tau*)"))
for (cens in cens_grid) {
  rr <- do.call(rbind, parallel::mclapply(seq_len(R), function(r) {
    d <- gen_one(N, 50000L + round(cens * 1000) + r, alpha_s, cens)
    art <- tryCatch(ar_tau_ci(d), error = function(e) NULL)
    if (is.null(art)) return(NULL)
    nd <- tryCatch(naive_and_delta(d), error = function(e) NULL)
    if (is.null(nd)) return(NULL)
    c(art_bounded = as.numeric(art$bounded), art_empty = as.numeric(art$empty),
      art_cov = as.numeric(!art$empty && tau_star >= art$ci[1] && tau_star <= art$ci[2]),
      naive_cov = as.numeric(tau_star >= nd$naive_ci[1] && tau_star <= nd$naive_ci[2]),
      delta_cov = as.numeric(tau_star >= nd$delta_ci[1] && tau_star <= nd$delta_ci[2]),
      art_width = if (art$bounded) diff(art$ci) else NA_real_)
  }, mc.cores = ncores))
  rr <- rr[stats::complete.cases(rr[, 1:5]), , drop = FALSE]
  cat(sprintf("%-8.2f %-14.1f %-16.3f %-16.3f %-16.3f\n",
              cens, 100 * mean(rr[, "art_bounded"]), mean(rr[, "art_cov"]),
              mean(rr[, "naive_cov"]), mean(rr[, "delta_cov"])))
}
cat("\n(art_bounded% = how often the AR-tau confidence set was a finite interval,\n")
cat(" not (-Inf,Inf) or empty -- an honest weak-identification diagnostic in its own right.)\n")
cat("DONE-P2-ARTAU\n")

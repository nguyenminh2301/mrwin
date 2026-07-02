#!/usr/bin/env Rscript
# Many-weak-instrument stress test for the one-sample over-ID AR test
# (mrwin_ar_onesample_overid()), Paper 03 open question (dev/p3-weak-iv-robust-
# winmr.md section 9): "does the over-ID test need an L-correction (a la
# Kleibergen) when L/n is non-negligible?" The originally validated probe used
# L=40 at N=2500 (L/N=0.016). This pushes L much higher relative to N.
#
# RESULT (type-I error under a valid instrument, R=60-100 reps/cell):
#   L    N      L/N     type-I
#   40   2500   0.016   0.020
#   100  2500   0.040   0.087
#   100  500    0.200   0.713
#   250  500    0.500   1.000
#
# FINDING: YES, the over-ID test needs an L-correction it does not have. Type-I
# error is calibrated (if slightly conservative) at L/N=0.016, already mildly
# inflated at L/N=0.04 (0.09 vs nominal 0.05), and catastrophically broken
# beyond that: 0.71 at L/N=0.2, and it ALWAYS rejects (1.00) at L/N=0.5. This
# mirrors the well-known many-weak-instrument problem for the linear-IV
# analogue (Kleibergen 2002, Econometrica), for which a correction exists in
# the linear case; none is implemented here. mrwin_ar_onesample_overid() now
# warns when L/n > 0.02 and refuses to run when L/n > 0.15, with thresholds set
# directly from this simulation (see ?mrwin_ar_onesample_overid). This is a
# genuine, currently-unresolved limitation of the shipped over-ID test, not
# fixed by this probe -- reported honestly rather than silently worked around.
#
# NOTE: the numbers above were obtained BEFORE the L/n safety guard existed.
# Re-running this script now will correctly show the L/N=0.5 cell refusing to
# run (a stop(), caught below and reported as "GUARDED") -- itself a live
# confirmation that the guard fires exactly where the original finding says it
# must. Public package + simulator only; fixed seeds.
suppressMessages(library(mrwin))
ncores <- max(1L, min(4L, parallel::detectCores()))
cfg0 <- mrwin_config(m_snps = 1L)
gen_pop <- function(n, seed, mafs, betas, alpha_s, cens = cfg0$censoring_rate) {
  set.seed(seed); M <- length(mafs)
  g <- sapply(mafs, function(p) rbinom(n, 2, p)); gstd <- scale(g); storage.mode(gstd) <- "double"
  s <- as.numeric(gstd %*% betas)
  u <- rnorm(n); w <- rgamma(n, 1 / cfg0$theta_f, scale = cfg0$theta_f)
  cn <- pmin(rexp(n, cens), cfg0$max_follow_up); uf <- matrix(runif(n * 3), n, 3)
  X <- alpha_s * s + cfg0$alpha_u * u + rnorm(n)
  list(X = X, s = s, gstd = gstd, u = u, w = w, cn = cn, uf = uf)
}
gen_out <- function(X, u, w, s, cn, uf, cfgc, extra = 0) {
  n <- length(X); et <- matrix(NA_real_, n, 3); lw <- log(w)
  for (j in 1:3) { lp <- cfgc$alpha_x[j] * X + cfgc$nu_u[j] * u + extra + lw
    et[, j] <- (-log(uf[, j]) / (cfgc$baseline_haz[j] * exp(lp)))^(1 / cfgc$shape_weibull) }
  td <- et[, 1]; th <- et[, 2]; tr <- et[, 3]
  list(time = cbind(pmin(td, cn), pmin(th, cn, td), pmin(tr, cn, td)),
       status = cbind(as.integer(td <= cn), as.integer(th <= cn & th <= td), as.integer(tr <= cn & tr <= td)))
}
run_L <- function(L, N, R, alpha_s = 0.4) {
  set.seed(7L); mafs <- runif(L, cfg0$maf_low, cfg0$maf_high); betas <- rnorm(L, 0, cfg0$sigma_beta)
  q <- stats::qchisq(0.95, L - 1L)
  guarded <- 0L
  rej <- unlist(parallel::mclapply(1:R, function(r) {
    p <- gen_pop(N, 9000L + r, mafs, betas, alpha_s)
    oc <- gen_out(p$X, p$u, p$w, p$s, p$cn, p$uf, cfg0)
    res <- tryCatch(suppressWarnings(mrwin_ar_onesample_overid(oc$time, oc$status, p$X, p$gstd)),
                     error = function(e) NULL)
    if (is.null(res)) return(NA)
    as.numeric(res$Q > q)
  }, mc.cores = ncores))
  if (L / N > 0.15) guarded <- sum(is.na(rej))
  c(L = L, N = N, LoverN = L / N, size = mean(rej, na.rm = TRUE), nfail = sum(is.na(rej)), guarded = guarded)
}
cat(sprintf("%-6s %-6s %-10s %-10s %-8s %-8s\n", "L", "N", "L/N", "type-I", "n_fail", "guarded"))
for (spec in list(c(40, 2500, 100), c(100, 2500, 80), c(100, 500, 80), c(250, 500, 60))) {
  out <- run_L(spec[1], spec[2], spec[3])
  status <- if (out["guarded"] == spec[3]) "(all runs correctly refused by the L/n>0.15 guard)" else ""
  cat(sprintf("%-6d %-6d %-10.3f %-10.3f %-8d %-8d %s\n", out["L"], out["N"], out["LoverN"],
              out["size"], out["nfail"], out["guarded"], status))
}
cat("DONEMANYWEAK\n")

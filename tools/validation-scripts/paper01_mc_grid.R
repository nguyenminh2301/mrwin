#!/usr/bin/env Rscript
# Paper-01 full Monte-Carlo validation grid: type-I, power, weak-instrument, and
# pleiotropy (SDPD). Complements the WP19 calibration cells (which covered type-I
# and coverage only). The primary reported CI is Fieller throughout (the
# calibrated interval; see inst/spec/validation-findings.md). Inference is the
# analytic influence-function path (validated to match the bootstrap, sigma_beta=0
# pure closed-form). Public package + simulator only; fixed seeds; parallel reps.
#
# Reads: nothing. Writes: a results table to stdout (rate +/- Monte-Carlo SE).
suppressMessages(library(mrwin))
ncores <- max(1L, min(4L, parallel::detectCores()))
cat(sprintf("paper01 MC grid | cores=%d | %s\n", ncores, "analytic + Fieller"))

# one simulated cohort -> the quantities every grid needs
rep_fit <- function(seed, N, m_snps, ax, alpha_s, gd, sigma_beta, sdpd, infer = "analytic") {
  cfg <- mrwin_config(n_outcome = N, m_snps = m_snps, alpha_s = alpha_s,
                      alpha_x = rep(ax, 3), gamma_direct = rep(gd, 3), seed = seed)
  dat <- mrwin_simulate(cfg, seed = seed)
  ep <- mrwin_endpoint(dat$time, dat$status, colnames(dat$time))
  gw <- mrwin_gwas(dat$true_betas, rep(sigma_beta, length(dat$true_betas)))
  be <- if (infer == "analytic") "dense" else "fast"
  ctrl <- mrwin_controls(n_strata = 5L, inference = infer, backend = be,
                         run_sdpd = sdpd, bootstrap = 150L, seed = seed)
  fit <- tryCatch(mrwin(endpoint = ep, genotype = dat$G, exposure = dat$X,
                        gwas = gw, controls = ctrl), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  ci <- fit$inference$ci95_delta_fieller
  bounded <- is.finite(ci[[1]]) && is.finite(ci[[2]])
  list(bounded = bounded,
       reject0 = bounded && ((ci[[1]] > 0) || (ci[[2]] < 0)),
       lo = ci[[1]], hi = ci[[2]],
       delta = fit$point$delta_gls,
       weak = any(vapply(fit$warnings, function(w) w$code == "weak_instrument", logical(1))),
       sdpd_rej = if (sdpd) isTRUE(fit$sdpd$rejected) else NA,
       underpow = if (sdpd) isTRUE(fit$sdpd$underpowered) else NA)
}

# run M reps in parallel, return the list of non-null fits
run_cell <- function(M, base, ...) {
  res <- parallel::mclapply(seq_len(M), function(m) rep_fit(base + m, ...),
                            mc.cores = ncores)
  res[!vapply(res, is.null, logical(1))]
}
mcse <- function(p, n) sqrt(p * (1 - p) / max(n, 1))

# truth (population DS-CWR delta) for coverage: mean delta_gls at large N
truth_of <- function(ax, alpha_s, base = 70000L, reps = 5L, N = 10000L) {
  d <- parallel::mclapply(seq_len(reps), function(i) {
    f <- rep_fit(base + i, N, 100L, ax, alpha_s, 0, 0, FALSE)
    if (is.null(f)) NA_real_ else f$delta
  }, mc.cores = ncores)
  mean(unlist(d), na.rm = TRUE)
}

flush_cat <- function(...) { cat(...); flush.console() }

## ====================== GRID A — TYPE-I (null) ======================
flush_cat("\n== A. TYPE-I (Fieller reject H0:delta=0 | ax=0, gd=0), target 0.05 ==\n")
for (N in c(2000L, 4000L)) {
  M <- 500L
  f <- run_cell(M, base = 1000L + N, N = N, m_snps = 100L, ax = 0, alpha_s = 0.4,
                gd = 0, sigma_beta = 0, sdpd = FALSE)
  r <- mean(vapply(f, `[[`, logical(1), "reject0"))
  flush_cat(sprintf("  N=%5d  type-I=%.3f (SE %.3f, n=%d)\n", N, r, mcse(r, length(f)), length(f)))
}

## ====================== GRID B — POWER (valid IV) ======================
flush_cat("\n== B. POWER (Fieller reject H0:delta=0 | true effect) ==\n")
flush_cat("  B1 by N (ax=-0.4, alpha_s=0.4):\n")
for (N in c(2000L, 4000L)) {
  M <- 400L
  f <- run_cell(M, base = 11000L + N, N = N, m_snps = 100L, ax = -0.4, alpha_s = 0.4,
                gd = 0, sigma_beta = 0, sdpd = FALSE)
  r <- mean(vapply(f, `[[`, logical(1), "reject0"))
  flush_cat(sprintf("    N=%5d  power=%.3f (SE %.3f, n=%d)\n", N, r, mcse(r, length(f)), length(f)))
}
flush_cat("  B2 by effect size (N=4000):\n")
for (ax in c(-0.1, -0.2, -0.4)) {
  M <- 400L
  f <- run_cell(M, base = as.integer(13000 + abs(ax) * 1000), N = 4000L, m_snps = 100L,
                ax = ax, alpha_s = 0.4, gd = 0, sigma_beta = 0, sdpd = FALSE)
  r <- mean(vapply(f, `[[`, logical(1), "reject0"))
  flush_cat(sprintf("    ax=%.2f  power=%.3f (SE %.3f, n=%d)\n", ax, r, mcse(r, length(f)), length(f)))
}

## ====================== GRID C — WEAK INSTRUMENT ======================
flush_cat("\n== C. WEAK-INSTRUMENT sweep (vary alpha_s, N=4000, ax=-0.4) ==\n")
for (as_ in c(0.1, 0.2, 0.4)) {
  M <- 300L
  tr <- truth_of(-0.4, as_, base = 80000L + as.integer(as_ * 100))
  f <- run_cell(M, base = 21000L + as.integer(as_ * 1000), N = 4000L, m_snps = 100L,
                ax = -0.4, alpha_s = as_, gd = 0, sigma_beta = 0, sdpd = FALSE)
  weakr <- mean(vapply(f, `[[`, logical(1), "weak"))
  unb <- mean(!vapply(f, `[[`, logical(1), "bounded"))
  bd <- vapply(f, `[[`, logical(1), "bounded")
  cov <- mean(mapply(function(g) g$bounded && g$lo <= tr && tr <= g$hi, f))
  bias <- mean(vapply(f, `[[`, numeric(1), "delta"), na.rm = TRUE) - tr
  flush_cat(sprintf("  alpha_s=%.1f truth=%.3f  weak_flag=%.2f  unbounded=%.2f  cover=%.3f  bias=%+.3f (n=%d)\n",
                    as_, tr, weakr, unb, cov, bias, length(f)))
}

## ====================== GRID D — PLEIOTROPY (SDPD) ======================
flush_cat("\n== D. PLEIOTROPY / SDPD rejection (N=4000, ax=-0.4) | gd=0 is type-I ==\n")
for (gd in c(0.0, 0.15, 0.3, 0.5)) {
  M <- 300L
  f <- run_cell(M, base = 31000L + as.integer(gd * 1000), N = 4000L, m_snps = 100L,
                ax = -0.4, alpha_s = 0.4, gd = gd, sigma_beta = 0, sdpd = TRUE)
  rej <- mean(vapply(f, `[[`, logical(1), "sdpd_rej"), na.rm = TRUE)
  up <- mean(vapply(f, `[[`, logical(1), "underpow"), na.rm = TRUE)
  lbl <- if (gd == 0) "type-I" else "power "
  flush_cat(sprintf("  gamma_direct=%.2f  SDPD reject (%s)=%.3f (SE %.3f)  underpowered=%.2f (n=%d)\n",
                    gd, lbl, rej, mcse(rej, length(f)), up, length(f)))
}
flush_cat("\nDONE\n")

#!/usr/bin/env Rscript
# Weak-instrument behaviour of the DS-CWR Fieller interval, done correctly. The
# main grid (paper01_mc_grid.R) used an estimator-based truth that is itself
# unstable at weak instruments (the ratio estimand log-theta/Delta-X has its
# denominator -> 0), and counted only BOUNDED Fieller intervals as covering. Both
# undersell the method. Here:
#   - truth = the package estimand at a STRONG instrument + large N (well
#     identified, consistent; the estimand is instrument-standardized so this is
#     the right instrument-independent target);
#   - coverage counts an UNBOUNDED Fieller interval as covering, which is the
#     honest convention: an unbounded CI is the method saying "I cannot bound
#     this" (it does not exclude the truth).
# We sweep instrument strength alpha_s and report the weak-instrument flag rate,
# the unbounded rate, and both coverage conventions. Public package + simulator.
suppressMessages(library(mrwin))
ncores <- max(1L, min(4L, parallel::detectCores()))

fit_cell <- function(seed, N, alpha_s, infer = "analytic", run_strong = FALSE) {
  cfg <- mrwin_config(n_outcome = N, m_snps = 100L, alpha_s = alpha_s,
                      alpha_x = rep(-0.4, 3), gamma_direct = c(0, 0, 0), seed = seed)
  dat <- mrwin_simulate(cfg, seed = seed)
  ep <- mrwin_endpoint(dat$time, dat$status, colnames(dat$time))
  gw <- mrwin_gwas(dat$true_betas, rep(0, length(dat$true_betas)))
  be <- if (run_strong) "fast" else if (infer == "analytic") "dense" else "fast"
  ctrl <- mrwin_controls(n_strata = 5L, inference = infer, backend = be,
                         run_sdpd = FALSE, bootstrap = 100L, seed = seed)
  fit <- tryCatch(mrwin(endpoint = ep, genotype = dat$G, exposure = dat$X,
                        gwas = gw, controls = ctrl), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  ci <- fit$inference$ci95_delta_fieller
  list(delta = fit$point$delta_gls, lo = ci[[1]], hi = ci[[2]],
       bounded = is.finite(ci[[1]]) && is.finite(ci[[2]]),
       weak = any(vapply(fit$warnings, function(w) w$code == "weak_instrument", logical(1))))
}

# truth: strong instrument (alpha_s=0.4) + large N, fast backend (point estimate)
tr <- parallel::mclapply(1:5, function(i) {
  f <- fit_cell(90000L + i, 20000L, 0.4, infer = "bootstrap", run_strong = TRUE)
  if (is.null(f)) NA_real_ else f$delta
}, mc.cores = ncores)
truth <- mean(unlist(tr), na.rm = TRUE)
cat(sprintf("truth (alpha_s=0.4, N=20000) = %.3f (rep sd %.3f)\n\n", truth, sd(unlist(tr), na.rm = TRUE)))

cat("weak-instrument sweep (N=4000, ax=-0.4), Fieller CI:\n")
for (as_ in c(0.10, 0.15, 0.20, 0.30, 0.40)) {
  M <- 300L
  res <- parallel::mclapply(seq_len(M), function(m) fit_cell(60000L + as.integer(as_*1000) + m, 4000L, as_),
                            mc.cores = ncores)
  res <- res[!vapply(res, is.null, logical(1))]
  bd <- vapply(res, `[[`, logical(1), "bounded")
  weak <- mean(vapply(res, `[[`, logical(1), "weak"))
  unb <- mean(!bd)
  cov_b <- mean(vapply(res, function(g) g$bounded && g$lo <= truth && truth <= g$hi, logical(1)))
  cov_t <- mean(vapply(res, function(g) (!g$bounded) || (g$lo <= truth && truth <= g$hi), logical(1)))
  bias <- median(vapply(res[bd], `[[`, numeric(1), "delta")) - truth   # bounded only
  cat(sprintf("  alpha_s=%.2f  weak_flag=%.2f  unbounded=%.2f  cover(bounded)=%.3f  cover(incl.unbounded)=%.3f  median_bias(bounded)=%+.3f (n=%d)\n",
              as_, weak, unb, cov_b, cov_t, bias, length(res)))
  flush.console()
}
cat("DONE\n")

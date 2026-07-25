#!/usr/bin/env Rscript
# ============================================================================
# PROBE #3, attempt 3 -- REPLACE the bootstrap, don't screen it.
#
# Attempts 1 (D_hat) and 2 (n*zeta1_min) both tried to PREDICT the bootstrap
# flag from a cheap proxy, and both failed at a useful false-safe rate. The
# framing was wrong: the flag is `c2_boot_ci[1] < 1`, i.e. a LOWER CONFIDENCE
# BOUND on c^2 = n*zeta1(gamma). So the right object to compute cheaply is not
# a correlate of the flag but the CONFIDENCE BOUND ITSELF.
#
# KEY OBSERVATION (from reading R/onesample_ar.R): the bootstrap holds `gamma`
# FIXED at the original point estimate inside the resampling loop --
#     max(mb$n * (mb$c0 - 2*gamma*mb$c1 + gamma^2*mb$c2), 0)
# -- so the bootstrapped quantity is exactly zeta1(gamma) = Var_i(g_i) at fixed
# gamma, with g_i := gh_i - gamma*gx_i. No delta-method term for gamma_hat's own
# randomness is needed. That makes an analytic SE available in closed form.
#
# DERIVATION. zeta1_hat = sample variance of g_1..g_n. The influence function of
# a variance functional is IF(g) = (g - mu)^2 - sigma^2, so
#     Var(zeta1_hat) = (m4 - m2^2)/n,   m2 := E[(g-mu)^2], m4 := E[(g-mu)^4]
# and since c^2 = n*zeta1_hat,
#     SE(c^2) = n * sqrt((m4 - m2^2)/n) = sqrt(n*(m4 - m2^2)).
# Two candidate lower bounds at one-sided level alpha1:
#   (i)  normal:    c2_hat - z*SE
#   (ii) lognormal: c2_hat * exp(-z*SE/c2_hat)   [variances are right-skewed,
#                   so a log-scale interval should track the bootstrap's own
#                   skew better -- tested, not assumed]
#
# CAVEAT THIS SCRIPT IS DESIGNED TO TEST: g_i are NOT iid draws -- gh_i and gx_i
# are U-statistic (leave-one-out-style) projections recomputed over the whole
# sample, which the bootstrap captures fully and the iid influence function
# above does not. Whether that gap matters is an empirical question, and is
# exactly what the comparison below measures. The decision-relevant metric is
# FLAG AGREEMENT (does the analytic bound produce the same degenerate/not call
# as the bootstrap), not correlation of the bounds.
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

# analytic c^2 lower bounds (normal + lognormal), reusing ONLY quantities the
# point estimate already needs -- no resampling.
analytic_c2_bounds <- function(d, alpha1 = 0.025) {
  N <- length(d$X)
  m <- blocks(d$time, d$status, d$X, d$Z)
  gamma <- m$a / m$b
  cw <- subj(d$time, d$status, d$time, d$status, rep(1, N)); s_i <- cw[, 1] - cw[, 2]
  cz <- subj(d$time, d$status, d$time, d$status, d$Z); r_i <- cz[, 1] - cz[, 2]
  gh_i <- (d$Z * s_i - r_i) / (N - 1)
  SX <- sum(d$X); SZ <- sum(d$Z); SXZ <- sum(d$X * d$Z)
  gx_i <- ((N - 1) * d$X * d$Z - d$X * (SZ - d$Z) - d$Z * (SX - d$X) + (SXZ - d$X * d$Z)) / (N - 1)
  g_i <- gh_i - gamma * gx_i
  gc <- g_i - mean(g_i)
  m2 <- mean(gc^2); m4 <- mean(gc^4)
  c2_hat <- max(N * stats::var(g_i), 0)
  se_c2 <- sqrt(max(N * (m4 - m2^2), 0))
  z <- stats::qnorm(1 - alpha1)
  lo_norm <- c2_hat - z * se_c2
  lo_lnorm <- if (c2_hat > 0) c2_hat * exp(-z * se_c2 / c2_hat) else 0
  list(c2_hat = c2_hat, se_c2 = se_c2, lo_norm = lo_norm, lo_lnorm = lo_lnorm, gamma = gamma)
}

one_rep <- function(n, seed, alpha_s, cens) {
  d <- gen_one(n, seed, alpha_s, cens)
  t0 <- proc.time()[["elapsed"]]
  an <- analytic_c2_bounds(d)
  t_an <- proc.time()[["elapsed"]] - t0
  t1 <- proc.time()[["elapsed"]]
  full <- mrwin_ar_onesample(d$time, d$status, d$X, d$Z, degeneracy_check = TRUE,
                              boot_reps = 200L, seed = seed)
  t_boot <- proc.time()[["elapsed"]] - t1
  c(c2_hat = an$c2_hat, se_c2 = an$se_c2,
    lo_norm = an$lo_norm, lo_lnorm = an$lo_lnorm,
    lo_boot = full$c2_boot_ci[1],
    flag_boot = as.numeric(full$degenerate),
    flag_norm = as.numeric(an$lo_norm < 1),
    flag_lnorm = as.numeric(an$lo_lnorm < 1),
    t_an = t_an, t_boot = t_boot)
}

cens_grid <- c(0.01, 0.05, 0.15, 0.40, 1.00, 2.50, 6.00)
N_grid <- c(800L, 2500L)
R <- 30L
cat("== PROBE 3 attempt 3: ANALYTIC c^2 lower bound vs the 200-rep bootstrap bound ==\n")
cat("(target = replace the bootstrap outright with the SAME decision rule, not a proxy screen)\n\n")
all_rows <- list()
for (alpha_s in c(0.05, 0.35)) {
  cat(sprintf("-- alpha_s=%.2f --\n", alpha_s))
  cat(sprintf("%-6s %-6s %-12s %-12s %-12s %-10s %-10s\n",
              "N", "cens", "lo_boot", "lo_norm", "lo_lnorm", "P(dgn.bt)", "speedup"))
  for (n in N_grid) for (cens in cens_grid) {
    rr <- do.call(rbind, parallel::mclapply(seq_len(R), function(r) {
      tryCatch(one_rep(n, 30000L * n + 7L * round(cens * 100) + r + round(alpha_s * 1000), alpha_s, cens),
               error = function(e) NULL)
    }, mc.cores = ncores))
    rr <- rr[stats::complete.cases(rr), , drop = FALSE]
    all_rows[[length(all_rows) + 1]] <- data.frame(N = n, cens = cens, alpha_s = alpha_s, rr)
    m <- colMeans(rr)
    cat(sprintf("%-6d %-6.2f %-12.4f %-12.4f %-12.4f %-10.3f %-10.1f\n",
                n, cens, m["lo_boot"], m["lo_norm"], m["lo_lnorm"], m["flag_boot"],
                m["t_boot"] / m["t_an"]))
  }
  cat("\n")
}
dat <- do.call(rbind, all_rows)
cat(sprintf("Total rows: %d\n\n", nrow(dat)))

cat("== Bound agreement (the analytic bound should TRACK the bootstrap bound) ==\n")
cat(sprintf("Pearson  cor(lo_boot, lo_norm)  = %.4f\n", cor(dat$lo_boot, dat$lo_norm)))
cat(sprintf("Pearson  cor(lo_boot, lo_lnorm) = %.4f\n", cor(dat$lo_boot, dat$lo_lnorm)))
cat(sprintf("Spearman cor(lo_boot, lo_lnorm) = %.4f\n", cor(dat$lo_boot, dat$lo_lnorm, method = "spearman")))
cat(sprintf("median(lo_norm / lo_boot)  = %.4f   median(lo_lnorm / lo_boot) = %.4f  (1.0 = unbiased)\n",
            median(dat$lo_norm / pmax(dat$lo_boot, 1e-12)),
            median(dat$lo_lnorm / pmax(dat$lo_boot, 1e-12))))

cat("\n== DECISION-RELEVANT metric: flag agreement with the bootstrap ==\n")
for (nm in c("flag_norm", "flag_lnorm")) {
  agree <- mean(dat[[nm]] == dat$flag_boot)
  # false-safe = analytic says NOT degenerate but bootstrap says degenerate (the dangerous error)
  fs <- sum(dat[[nm]] == 0 & dat$flag_boot == 1)
  fs_rate <- if (sum(dat[[nm]] == 0) > 0) fs / sum(dat[[nm]] == 0) else NA
  fc <- sum(dat[[nm]] == 1 & dat$flag_boot == 0)
  cat(sprintf("%-11s agreement=%.4f   false-SAFE=%d (%.1f%% of its 'safe' calls)   false-ALARM=%d\n",
              nm, agree, fs, 100 * fs_rate, fc))
}
cat(sprintf("\nMean speedup (analytic vs 200-rep bootstrap): %.1fx\n", mean(dat$t_boot / dat$t_an)))
cat("DONE-P3-ANALYTIC\n")

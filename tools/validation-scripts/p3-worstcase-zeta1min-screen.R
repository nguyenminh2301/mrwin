#!/usr/bin/env Rscript
# ============================================================================
# P3, a second attempt at a cheap degeneracy pre-screen, after the D_hat screen
# (p3-dthreshold-degeneracy-screen.R) failed: D_hat alone cannot see the
# weak-instrument point-estimate wandering that dominates zeta1(beta_hat).
#
# KEY IDEA: zeta1(beta) = c0 - 2*beta*c1 + beta^2*c2 is an upward parabola in
# beta (c2=Var(gx_i)>=0). Its MINIMUM over ALL beta (not just the estimated
# beta_hat) is
#   zeta1_min := c0 - c1^2/c2   (at beta = c1/c2)
# and by construction zeta1(beta_hat) >= zeta1_min for ANY realized beta_hat,
# NO MATTER how far weak-instrument noise pushes beta_hat from the true value.
# So a threshold on n*zeta1_min gives a DETERMINISTIC (not merely empirical)
# lower bound on the actual bootstrapped quantity c^2=n*zeta1(beta_hat) --
# it should be immune to the exact failure mode that broke the D_hat screen.
# zeta1_min costs NOTHING extra: c0,c1,c2 are already computed by
# .mrwin_ar1_blocks() for the point estimate itself.
#
# TEST: same DGP/grid as p3-dthreshold-degeneracy-screen.R (weak + strong IV,
# N in {800,2500}, 7 censoring rates), ground truth = mrwin_ar_onesample's own
# bootstrap `degenerate` flag (boot_reps=200). Does n*zeta1_min > threshold
# reliably rule out degenerate==TRUE, INCLUDING in the weak-instrument cells
# where D_hat failed?
# ============================================================================
suppressMessages(library(mrwin))
ncores <- max(1L, min(4L, parallel::detectCores()))
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
  m <- blocks(d$time, d$status, d$X, d$Z)
  zeta1_min <- m$c0 - m$c1^2 / m$c2   # cheap, deterministic lower bound on zeta1(beta_hat), any beta_hat
  full <- mrwin_ar_onesample(d$time, d$status, d$X, d$Z, degeneracy_check = TRUE,
                              boot_reps = 200L, seed = seed)
  c(zeta1_min = zeta1_min, n_zeta1_min = n * zeta1_min, c2_hat_full = full$c2_hat,
    degenerate = as.numeric(full$degenerate))
}

cens_grid <- c(0.01, 0.05, 0.15, 0.40, 1.00, 2.50, 6.00)
N_grid <- c(800L, 2500L)
R <- 30L
cat("== n*zeta1_min (cheap, deterministic lower bound) vs bootstrap degenerate flag ==\n")
all_rows <- list()
for (alpha_s in c(0.05, 0.35)) {
  cat(sprintf("\n-- alpha_s=%.2f --\n", alpha_s))
  cat(sprintf("%-6s %-6s %-14s %-12s\n", "N", "cens", "n*zeta1_min", "P(degen)"))
  for (n in N_grid) for (cens in cens_grid) {
    rr <- do.call(rbind, parallel::mclapply(seq_len(R), function(r) {
      tryCatch(one_rep(n, 20000L * n + 7L * round(cens * 100) + r + round(alpha_s * 1000), alpha_s, cens),
               error = function(e) NULL)
    }, mc.cores = ncores))
    rr <- rr[stats::complete.cases(rr), , drop = FALSE]
    all_rows[[length(all_rows) + 1]] <- data.frame(N = n, cens = cens, alpha_s = alpha_s, rr)
    m <- colMeans(rr)
    cat(sprintf("%-6d %-6.2f %-14.4f %-12.3f\n", n, cens, m["n_zeta1_min"], m["degenerate"]))
  }
}
dat <- do.call(rbind, all_rows)
cat(sprintf("\nTotal rows: %d\n", nrow(dat)))
cat(sprintf("Spearman(n*zeta1_min, degenerate) = %.4f\n", cor(dat$n_zeta1_min, dat$degenerate, method = "spearman")))

cat("\n== Screen calibration: does n*zeta1_min > tau ever co-occur with degenerate==1? ==\n")
cat(sprintf("%-10s %-14s %-16s %-16s\n", "tau", "n_above_tau", "n_degen_above", "coverage(skip%)"))
for (tau in c(1, 2, 5, 10, 20, 50, 100, 200)) {
  above <- dat$n_zeta1_min > tau
  n_above <- sum(above); n_bad <- sum(above & dat$degenerate == 1)
  cat(sprintf("%-10.1f %-14d %-16d %-16.1f\n", tau, n_above, n_bad, 100 * n_above / nrow(dat)))
}
cat("\nMinimum n*zeta1_min observed among degenerate==1 rows (should be small if the screen is valid):\n")
deg1 <- dat[dat$degenerate == 1, ]
cat(sprintf("  min=%.4f  median=%.4f  max=%.4f  (n=%d degenerate rows)\n",
            min(deg1$n_zeta1_min), median(deg1$n_zeta1_min), max(deg1$n_zeta1_min), nrow(deg1)))
cat("Maximum n*zeta1_min observed among degenerate==0 rows (context for choosing tau):\n")
deg0 <- dat[dat$degenerate == 0, ]
cat(sprintf("  min=%.4f  median=%.4f  max=%.4f  (n=%d non-degenerate rows)\n",
            min(deg0$n_zeta1_min), median(deg0$n_zeta1_min), max(deg0$n_zeta1_min), nrow(deg0)))
cat("DONE-ZETA1MIN\n")

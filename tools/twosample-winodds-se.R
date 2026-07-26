#!/usr/bin/env Rscript
# Closed-form SE for the per-SNP win-odds coefficient, and its payoff for
# two-sample win-ratio MR (real inverse-variance IVW weights + valid CIs).
# -----------------------------------------------------------------------------
# mrwin_win_snp() estimates the per-SNP log-win-odds slope on dosage with the
# influence-function variance of log(W/L) (the estimator used by
# mrwin_analytic_covariance), evaluated subquadratically: the point estimate is
# on the fast kernel, the variance from per-subject win/loss counts on a capped
# subject subsample (mrwin_subject_win_loss_cpp).
#
# Two validations against ground truth:
#   1. per-SNP: closed-form se_log_theta vs the replication SD of log-theta.
#   2. IVW: mrwin_twosample_ivw() SE vs the empirical SD of gamma_hat, and 95% CI
#      coverage of the interventional oracle gamma*.
# Public package + simulator only.
suppressMessages(library(mrwin))
cfg <- mrwin_config(m_snps = 40L)

gen_out <- function(X, u, w, s, cn, uf, cfg) {
  n <- length(X); et <- matrix(NA_real_, n, 3); lw <- log(w)
  for (j in 1:3) {
    lp <- cfg$alpha_x[j]*X + cfg$nu_u[j]*u + cfg$gamma_direct[j]*s + lw
    et[, j] <- (-log(uf[, j]) / (cfg$baseline_haz[j]*exp(lp)))^(1/cfg$shape_weibull)
  }
  td <- et[,1]; th <- et[,2]; tr <- et[,3]
  list(time = cbind(pmin(td,cn), pmin(th,cn,td), pmin(tr,cn,td)),
       status = cbind(as.integer(td<=cn), as.integer(th<=cn & th<=td), as.integer(tr<=cn & tr<=td)))
}
gen_pop <- function(n, seed, mafs, betas) {
  set.seed(seed)
  graw <- sapply(mafs, function(p) rbinom(n, 2, p)); g <- scale(graw); storage.mode(g) <- "double"
  s <- drop(g %*% betas); u <- rnorm(n); w <- rgamma(n, 1/cfg$theta_f, scale = cfg$theta_f)
  cn <- pmin(rexp(n, cfg$censoring_rate), cfg$max_follow_up); uf <- matrix(runif(n*3), n, 3)
  X <- cfg$alpha_s*s + cfg$alpha_u*u + rnorm(n); oc <- gen_out(X, u, w, s, cn, uf, cfg)
  list(time = oc$time, status = oc$status, graw = graw, X = X, s = s, u = u, w = w, cn = cn, uf = uf)
}
WLr <- function(Ht, Hs, Lt, Ls) {
  v <- mrwin_fast_pair_win_loss(Ht, Hs, Lt, Ls)
  log((as.numeric(v[["wins"]]) + .5) / (as.numeric(v[["losses"]]) + .5))
}
set.seed(7); mafs <- runif(40, cfg$maf_low, cfg$maf_high); betas <- rnorm(40, 0, cfg$sigma_beta)

## ---- 1. per-SNP closed-form SE vs replication SD of log-theta ----
k <- which.min(abs(mafs - 0.3))
p0 <- gen_pop(4000L, 101L, mafs, betas)
r0 <- mrwin_win_snp(p0$time, p0$status, p0$graw[, k], max_subjects = 10000L)  # exact (no subsample)
R <- 300; lt <- numeric(R)
for (r in 1:R) {
  pr <- gen_pop(4000L, 20000L + r, mafs, betas)
  rr <- mrwin_win_snp(pr$time, pr$status, pr$graw[, k], max_subjects = 10000L)
  lt[r] <- rr$log_theta
}
cat(sprintf("per-SNP (SNP %d):  closed-form se_log_theta = %.4f   replication SD = %.4f   ratio = %.2f\n\n",
            k, r0$se_log_theta, sd(lt), r0$se_log_theta / sd(lt)))

## ---- 2. two-sample IVW: SE calibration + CI coverage of oracle gamma* ----
po <- gen_pop(150000L, 99L, mafs, betas); eps <- 0.25
Bp <- gen_out(po$X+eps, po$u,po$w,po$s,po$cn,po$uf, cfg)
Bm <- gen_out(po$X-eps, po$u,po$w,po$s,po$cn,po$uf, cfg)
gstar <- WLr(Bp$time, Bp$status, Bm$time, Bm$status) / (2*eps)
pe <- gen_pop(150000L, 555L, mafs, betas)
bGX <- apply(pe$graw, 2, function(d) cov(pe$X, d) / var(d))   # fixed exposure GWAS

Rc <- 80; Nout <- 10000L; cap <- 700L
G <- S <- CL <- CU <- numeric(Rc)
for (r in 1:Rc) {
  po <- gen_pop(Nout, 30000L + r, mafs, betas)
  wg <- mrwin_win_gwas(po$time, po$status, po$graw, max_subjects = cap, seed = r)
  iv <- mrwin_twosample_ivw(bGX, wg$delta, wg$se)
  G[r] <- iv$gamma; S[r] <- iv$se; CL[r] <- iv$ci95[1]; CU[r] <- iv$ci95[2]
}
cat(sprintf("two-sample IVW (R=%d, Nout=%d, oracle gamma*=%.3f):\n", Rc, Nout, gstar))
cat(sprintf("  mean gamma_hat = %.3f (bias %+.3f, weak-IV; -> 0 as N grows)\n", mean(G), mean(G)-gstar))
cat(sprintf("  mean estimated SE = %.3f   empirical SD = %.3f   ratio = %.2f\n", mean(S), sd(G), mean(S)/sd(G)))
cat(sprintf("  95%% CI coverage of gamma* = %.1f%%\n", 100*mean(CL <= gstar & gstar <= CU)))

#!/usr/bin/env Rscript
# Feasibility probe: TWO-SAMPLE / summary-data win-ratio MR.
# -----------------------------------------------------------------------------
# Question: can a hierarchical-endpoint win-ratio MR effect be estimated from
# *summary statistics* (no individual-level data shared across cohorts), the way
# ordinary two-sample MR pools per-SNP (beta_GX, beta_GY) by IVW?
#
# Insight: the outcome-side summary statistic is a per-SNP WIN-ODDS coefficient
# delta_G,win (a "win-odds GWAS"): the log-win-odds slope of the priority-ranked
# outcome on allele dosage, computed in the outcome cohort and shared as a
# summary. Under a valid instrument and local linearity, delta_G,m ~= gamma *
# beta_GX,m, so fixed-effect IVW of the per-SNP coefficients recovers the causal
# win-odds gradient gamma. The interventional oracle (paper 01) supplies gamma*.
#
# Finding (see dev/findings-two-sample-winratio.md): FEASIBLE — IVW from two
# independent cohorts recovers the right sign and magnitude, with a modest
# downward attenuation (~15-30%) whose source (weak instrument vs win-odds
# non-collapsibility) is the open question for a full method.
#
# Public package + simulator only; no confidential data.
suppressMessages(library(mrwin))

gen_out <- function(X, u, w, s, cn, uf, cfg) {
  n <- length(X); et <- matrix(NA_real_, n, 3); lw <- log(w)
  for (j in 1:3) {
    lp <- cfg$alpha_x[j]*X + cfg$nu_u[j]*u + cfg$gamma_direct[j]*s + lw
    et[, j] <- (-log(uf[, j]) / (cfg$baseline_haz[j]*exp(lp)))^(1/cfg$shape_weibull)
  }
  td <- et[,1]; th <- et[,2]; tr <- et[,3]
  time <- cbind(pmin(td,cn), pmin(th,cn,td), pmin(tr,cn,td))
  st <- cbind(as.integer(td<=cn), as.integer(th<=cn & th<=td), as.integer(tr<=cn & tr<=td))
  colnames(time) <- colnames(st) <- c("death","hf","renal"); list(time=time, status=st)
}
# Fixed genetic architecture (mafs, betas) shared across samples; only the
# individuals differ between cohorts -- exactly the two-sample MR setting.
gen_pop <- function(n, cfg, seed, mafs, betas) {
  set.seed(seed)
  g <- scale(sapply(mafs, function(p) rbinom(n, 2, p))); storage.mode(g) <- "double"
  s <- drop(g %*% betas); u <- rnorm(n); w <- rgamma(n, 1/cfg$theta_f, scale=cfg$theta_f)
  cn <- pmin(rexp(n, cfg$censoring_rate), cfg$max_follow_up); uf <- matrix(runif(n*3), n, 3)
  X <- cfg$alpha_s*s + cfg$alpha_u*u + rnorm(n); oc <- gen_out(X, u, w, s, cn, uf, cfg)
  list(time=oc$time, status=oc$status, G=g, X=X, s=s, u=u, w=w, cn=cn, uf=uf)
}
WL <- function(th,sh,tl,sl) {
  v <- mrwin_fast_pair_win_loss(th,sh,tl,sl)
  log(as.numeric(v[["wins"]]) / as.numeric(v[["losses"]]))
}
beta_GX  <- function(p) apply(p$G, 2, function(gk) cov(p$X, gk) / var(gk))   # SNP -> exposure
delta_Gwin <- function(p) apply(p$G, 2, function(gk) {                       # SNP -> win-odds
  hi <- which(gk > median(gk)); lo <- which(gk <= median(gk))
  if (length(hi) < 20 || length(lo) < 20) return(NA_real_)
  WL(p$time[hi,,drop=FALSE], p$status[hi,,drop=FALSE],
     p$time[lo,,drop=FALSE], p$status[lo,,drop=FALSE]) / (mean(gk[hi]) - mean(gk[lo]))
})
ivw <- function(bx, dy) { ok <- is.finite(bx) & is.finite(dy); sum(bx[ok]*dy[ok]) / sum(bx[ok]^2) }

cfg <- mrwin_config(m_snps = 40L)
set.seed(7); mafs <- runif(40, cfg$maf_low, cfg$maf_high); betas <- rnorm(40, 0, cfg$sigma_beta)

## oracle gamma* (interventional, secant over the PRS-stratum exposure span)
po <- gen_pop(120000L, cfg, 99L, mafs, betas); ref <- gen_pop(20000L, cfg, 1L, mafs, betas)
sb <- findInterval(ref$s, quantile(ref$s, seq(0,1,length.out=6)), rightmost.closed=TRUE, all.inside=TRUE)
xlo <- mean(ref$X[sb==1]); xhi <- mean(ref$X[sb==5])
A <- gen_out(rep(xlo,120000L), po$u,po$w,po$s,po$cn,po$uf,cfg)
B <- gen_out(rep(xhi,120000L), po$u,po$w,po$s,po$cn,po$uf,cfg)
gstar <- WL(B$time,B$status,A$time,A$status) / (xhi - xlo)
cat(sprintf("oracle gamma* = %.3f\n\n", gstar))

## two independent samples (exposure-GWAS cohort A; win-GWAS cohort B), bias vs N
cat("two-sample IVW of per-SNP win-odds, bias vs N (R=6):\n")
for (N in c(8000L, 16000L, 32000L)) {
  e <- numeric(6)
  for (r in 1:6) {
    pa <- gen_pop(N, cfg, 1000L + N + r, mafs, betas)   # exposure GWAS
    pb <- gen_pop(N, cfg, 9000L + N + r, mafs, betas)   # win-odds GWAS (independent)
    e[r] <- ivw(beta_GX(pa), delta_Gwin(pb))
  }
  cat(sprintf("  N=%6d  gamma_hat=%.3f (se %.3f)  bias=%+.3f\n",
              N, mean(e), sd(e)/sqrt(6), mean(e) - gstar))
}

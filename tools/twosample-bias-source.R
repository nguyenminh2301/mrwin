#!/usr/bin/env Rscript
# Source of the two-sample win-ratio MR attenuation: weak instrument vs
# win-odds non-collapsibility. (Follow-up to tools/twosample-feasibility.R.)
# -----------------------------------------------------------------------------
# The feasibility probe showed IVW of per-SNP win-odds coefficients recovers the
# causal gradient but with a downward attenuation at small N. Two competing
# explanations:
#   H_weak  : finite-sample weak-instrument bias  -> vanishes as N grows.
#   H_ncoll : win-odds non-collapsibility          -> structural, persists at N=inf
#             and scales with residual outcome heterogeneity.
#
# Three decisive tests, all against the interventional oracle gamma* (do(X+/-eps)):
#   1. convergence in N (two independent samples): does bias -> 0 ?
#   2. per-SNP linearity: is delta_G,m ~ gamma * beta_GX,m (the IVW assumption)?
#   3. attenuation vs residual heterogeneity: does the ratio IVW/gamma* depend on
#      the residual outcome spread (the signature of non-collapsibility)?
#
# RESULT (see dev/findings-two-sample-winratio.md): H_weak. The estimator is
# CONSISTENT -- bias is ordinary weak-instrument finite-sample bias (toward null
# at small N), the per-SNP coefficient is linear in beta_GX with slope gamma*,
# and the IVW/gamma* ratio is ~1 regardless of heterogeneity. Non-collapsibility
# is ruled out as the mechanism.
suppressMessages(library(mrwin))
cfg <- mrwin_config(m_snps = 40L)

gen_out <- function(X, u, w, s, cn, uf, cfg, nu_sc = 1, gd_sc = 1) {
  n <- length(X); et <- matrix(NA_real_, n, 3); lw <- log(w)
  for (j in 1:3) {
    lp <- cfg$alpha_x[j]*X + nu_sc*cfg$nu_u[j]*u + gd_sc*cfg$gamma_direct[j]*s + lw
    et[, j] <- (-log(uf[, j]) / (cfg$baseline_haz[j]*exp(lp)))^(1/cfg$shape_weibull)
  }
  td <- et[,1]; th <- et[,2]; tr <- et[,3]
  time <- cbind(pmin(td,cn), pmin(th,cn,td), pmin(tr,cn,td))
  st <- cbind(as.integer(td<=cn), as.integer(th<=cn & th<=td), as.integer(tr<=cn & tr<=td))
  colnames(time) <- colnames(st) <- c("death","hf","renal"); list(time=time, status=st)
}
# nu_sc/gd_sc/frail_sc scale residual outcome heterogeneity for test 3.
gen_pop <- function(n, seed, mafs, betas, nu_sc = 1, gd_sc = 1, frail_sc = 1) {
  set.seed(seed)
  graw <- sapply(mafs, function(p) rbinom(n, 2, p)); g <- scale(graw); storage.mode(g) <- "double"
  s <- drop(g %*% betas); u <- rnorm(n); thf <- cfg$theta_f*frail_sc; w <- rgamma(n, 1/thf, scale=thf)
  cn <- pmin(rexp(n, cfg$censoring_rate), cfg$max_follow_up); uf <- matrix(runif(n*3), n, 3)
  X <- cfg$alpha_s*s + cfg$alpha_u*u + rnorm(n); oc <- gen_out(X, u, w, s, cn, uf, cfg, nu_sc, gd_sc)
  list(time=oc$time, status=oc$status, graw=graw, X=X, s=s, u=u, w=w, cn=cn, uf=uf)
}
WLr <- function(th, sh, tl, sl) {                      # +0.5 continuity guards log(0)
  v <- mrwin_fast_pair_win_loss(th, sh, tl, sl)
  log((as.numeric(v[["wins"]]) + 0.5) / (as.numeric(v[["losses"]]) + 0.5))
}
oracle <- function(N, seed, mafs, betas, nu_sc = 1, gd_sc = 1, frail_sc = 1, eps = 0.25) {
  p <- gen_pop(N, seed, mafs, betas, nu_sc, gd_sc, frail_sc)
  Bp <- gen_out(p$X+eps, p$u,p$w,p$s,p$cn,p$uf, cfg, nu_sc, gd_sc)
  Bm <- gen_out(p$X-eps, p$u,p$w,p$s,p$cn,p$uf, cfg, nu_sc, gd_sc)
  WLr(Bp$time, Bp$status, Bm$time, Bm$status) / (2*eps)
}
sub   <- function(p, idx) list(time=p$time[idx,,drop=FALSE], status=p$status[idx,,drop=FALSE])
bGXf  <- function(p) apply(p$graw, 2, function(d) cov(p$X, d) / var(d))          # per raw allele
dMedf <- function(p) apply(p$graw, 2, function(d) {                              # win-odds slope
  hi <- which(d > median(d)); lo <- which(d <= median(d))
  if (length(hi) < 50 || length(lo) < 50) return(NA_real_)
  H <- sub(p, hi); L <- sub(p, lo)
  WLr(H$time, H$status, L$time, L$status) / (mean(d[hi]) - mean(d[lo]))
})
ivw <- function(b, d) { ok <- is.finite(b) & is.finite(d); sum(b[ok]*d[ok]) / sum(b[ok]^2) }

set.seed(7); mafs <- runif(40, cfg$maf_low, cfg$maf_high); betas <- rnorm(40, 0, cfg$sigma_beta)
g0 <- oracle(200000L, 99L, mafs, betas)
cat(sprintf("oracle gamma*_local = %.3f\n\n", g0))

cat("test 1 -- two independent samples, convergence in N (R=4):\n")
for (N in c(16000L, 50000L, 150000L)) {
  e <- numeric(4)
  for (r in 1:4) {
    pa <- gen_pop(N, 2000L+N+r, mafs, betas); pb <- gen_pop(N, 7000L+N+r, mafs, betas)
    e[r] <- ivw(bGXf(pa), dMedf(pb))
  }
  cat(sprintf("  N=%6d  IVW=%.3f (se %.3f)  bias=%+.3f\n", N, mean(e), sd(e)/sqrt(4), mean(e)-g0))
}

cat("\ntest 2 -- per-SNP linearity delta_med ~ beta_GX (N=150k, one sample):\n")
p <- gen_pop(150000L, 123L, mafs, betas); b <- bGXf(p); d <- dMedf(p); ok <- is.finite(b) & is.finite(d)
cat(sprintf("  slope(thru origin)=%.3f  cor=%.3f  n=%d\n",
            sum(b[ok]*d[ok])/sum(b[ok]^2), cor(b[ok], d[ok]), sum(ok)))

cat("\ntest 3 -- attenuation vs residual heterogeneity (N=150k):\n")
for (h in c(1.0, 0.4, 0.1)) {
  gh <- oracle(150000L, 99L, mafs, betas, nu_sc=h, gd_sc=h, frail_sc=h)
  p  <- gen_pop(150000L, 8000L, mafs, betas, nu_sc=h, gd_sc=h, frail_sc=h)
  ga <- ivw(bGXf(p), dMedf(p))
  cat(sprintf("  het_scale=%.1f  gamma*=%.3f  IVW=%.3f  ratio=%.3f\n", h, gh, ga, ga/gh))
}

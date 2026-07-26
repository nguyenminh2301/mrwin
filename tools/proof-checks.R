#!/usr/bin/env Rscript
# Numerical checks for the two non-trivial premises of Proposition 1 (Identification).
#
#  V1  Covariance structure of the adjacent instrument-standardized gradients
#      (ISGs) delta_1..delta_{D-1}. Lemma "Two-sample U-statistic limits"
#      predicts a block-tridiagonal covariance (adjacent pairs share one stratum;
#      non-adjacent pairs share none -> independent). We confirm the off-diagonal
#      couplings are negligible, so GLS is in practice inverse-variance weighting.
#
#  V2  Local linearity of the causal win-log-odds Lambda^c(x', x) -- the exact
#      content of assumption (M). We evaluate the interventional oracle on an
#      exposure grid and confirm Lambda^c is linear in the exposure contrast
#      through the origin (R^2 ~ 1), with slope = gamma*.
#
# Public package + simulator only; no confidential data.   Rscript tools/proof-checks.R
suppressMessages(library(mrwin))

gen_outcomes <- function(X, u, w, s, cn, uf, cfg) {
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
gen_pop <- function(n, cfg, seed) {
  set.seed(seed); m <- cfg$m_snps
  mafs <- runif(m, cfg$maf_low, cfg$maf_high)
  g <- scale(sapply(mafs, function(p) rbinom(n,2,p))); storage.mode(g) <- "double"
  b <- rnorm(m, 0, cfg$sigma_beta); s <- drop(g %*% b)
  u <- rnorm(n); w <- rgamma(n, 1/cfg$theta_f, scale=cfg$theta_f)
  cn <- pmin(rexp(n, cfg$censoring_rate), cfg$max_follow_up); uf <- matrix(runif(n*3), n, 3)
  X <- cfg$alpha_s*s + cfg$alpha_u*u + rnorm(n)
  oc <- gen_outcomes(X, u, w, s, cn, uf, cfg)
  list(time=oc$time, status=oc$status, X=X, s=s, u=u, w=w, cn=cn, uf=uf)
}
WL <- function(th,sh,tl,sl) {
  v <- mrwin_fast_pair_win_loss(th,sh,tl,sl)
  c(w = as.numeric(v[["wins"]]), l = as.numeric(v[["losses"]]))
}
qs <- function(z, D) findInterval(z, quantile(z, seq(0,1,length.out=D+1)),
                                  rightmost.closed=TRUE, all.inside=TRUE)
cfg <- mrwin_config(m_snps = 40L); D <- 5L

## V1 -- empirical covariance of adjacent ISGs ---------------------------------
isg_vec <- function(p, D) {
  s <- qs(p$s, D); v <- numeric(D-1)
  for (d in 1:(D-1)) {
    hi <- which(s == d+1); lo <- which(s == d)
    wl <- WL(p$time[hi,,drop=FALSE], p$status[hi,,drop=FALSE],
             p$time[lo,,drop=FALSE], p$status[lo,,drop=FALSE])
    v[d] <- log(wl["w"]/wl["l"]) / (mean(p$X[hi]) - mean(p$X[lo]))
  }
  v
}
R <- 300L; M <- matrix(NA_real_, R, D-1)
for (r in 1:R) M[r, ] <- isg_vec(gen_pop(8000L, cfg, 5000L + r), D)
M <- M[is.finite(rowSums(M)), ]
C <- cor(M); off <- abs(row(C) - col(C))
cat(sprintf("V1: correlation matrix of (delta_1..delta_%d), R=%d cohorts\n", D-1, nrow(M)))
print(round(C, 3))
cat(sprintf("  adjacent mean|corr| = %.3f | gap>=2 mean|corr| = %.3f | MC s.e. ~ %.3f\n\n",
            mean(abs(C[off==1])), mean(abs(C[off>=2])), 1/sqrt(nrow(M))))

## V2 -- linearity of the causal win-log-odds (assumption M) -------------------
Lc <- function(x1, x0, n, seed) {
  p <- gen_pop(n, cfg, seed)
  A <- gen_outcomes(rep(x0,n), p$u,p$w,p$s,p$cn,p$uf,cfg)
  B <- gen_outcomes(rep(x1,n), p$u,p$w,p$s,p$cn,p$uf,cfg)
  wl <- WL(B$time,B$status,A$time,A$status); log(wl["w"]/wl["l"])
}
grid <- seq(-0.18, 0.18, by = 0.06); x0 <- -0.18
lam <- sapply(grid, function(x) Lc(x, x0, 120000L, 99L))
fit <- lm(lam ~ I(grid - x0))
cat("V2: Lambda^c(x, x0=-0.18) vs exposure contrast (n=120k paired)\n")
print(round(rbind(contrast = grid - x0, Lambda_c = as.numeric(lam)), 3))
cat(sprintf("  linear fit: slope = %.4f, intercept = %.4f, R^2 = %.5f\n",
            coef(fit)[2], coef(fit)[1], summary(fit)$r.squared))

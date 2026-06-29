#!/usr/bin/env Rscript
# BREAKTHROUGH PROBE (Paper 04 candidate): within-family / sibling win-ratio MR.
#
# Key structural identity: the win statistic is a PAIRWISE comparison h(O_i,O_j), and
# the sibling design is a PAIRWISE genetic contrast (Z_i - Z_j between sibs). They are
# the SAME object. Restricting the Paper-3 pairwise IV moment to sibling pairs makes
# the contrast obey Mendelian segregation, so dZ is INDEPENDENT of population
# stratification / assortative mating / dynastic (parental) effects -- the confounders
# population (all-individual) win-MR cannot control.
#
# Falsifiable test under population stratification (stratum raises trait-RAISING
# alleles -> systematic PRS shift, AND worsens the outcome -> a pure Z<->stratum->Y
# backdoor; stratum does NOT affect the exposure):
#   - population win-MR  beta = Cov(Z, w(O)/(n-1)) / Cov(Z, X)         [should be biased]
#   - within-family      beta = sum_sib h*dZ / sum_sib dX*dZ           [should be unbiased]
#   - oracle             interventional net-benefit gradient do(X +/- eps).
#
# RESULT (R=16, F=4000): oracle 0.099; population -0.794 (bias -0.893, sign flipped by
# the stratification confounding); within-family 0.098 (bias -0.001 -- recovers the
# causal effect). The pairwise win moment restricted to sibling pairs is robust to
# population stratification; population win-MR is catastrophically biased.
suppressMessages(library(mrwin))
ncores <- max(1L, min(4L, parallel::detectCores())); cfg <- mrwin_config()
subj <- getFromNamespace("mrwin_subject_win_loss_cpp", "mrwin")

# hierarchical win h(a,b) for matched sibling vectors (vectorized over pairs), K=3
hpair <- function(t1, s1, t2, s2) {
  F <- nrow(t1); res <- integer(F)
  for (k in 1:3) { dec <- res != 0
    iwin <- (s2[, k] == 1 & t1[, k] > t2[, k]); ilos <- (s1[, k] == 1 & t2[, k] > t1[, k])
    res[!dec & iwin] <- 1L; res[!dec & ilos & !iwin] <- -1L }
  res
}
outcome <- function(X, u, w, strat, cn, uf, gstrat) {
  F <- length(X); et <- matrix(NA_real_, F, 3); lw <- log(w)
  for (j in 1:3) { lp <- cfg$alpha_x[j]*X + cfg$nu_u[j]*u + gstrat*strat + lw
    et[, j] <- (-log(uf[, j]) / (cfg$baseline_haz[j]*exp(lp)))^(1/cfg$shape_weibull) }
  td <- et[,1]; th <- et[,2]; tr <- et[,3]
  list(time = cbind(pmin(td,cn), pmin(th,cn,td), pmin(tr,cn,td)),
       status = cbind(as.integer(td<=cn), as.integer(th<=cn & th<=td), as.integer(tr<=cn & tr<=td)))
}
transmit <- function(G) ifelse(G == 2, 1, ifelse(G == 0, 0, rbinom(length(G), 1, 0.5)))
gen_fam <- function(Fn, seed, gstrat = 2.5, fgap = 0.18) {
  set.seed(seed); M <- 30; base <- runif(M, 0.15, 0.35); betas <- rnorm(M, 0, 0.15)
  strat <- rbinom(Fn, 1, 0.5)
  p <- matrix(base, Fn, M, byrow = TRUE) + outer(strat, fgap*sign(betas)); p <- pmin(pmax(p, 0.02), 0.95)
  Gf <- matrix(rbinom(Fn*M, 2, p), Fn, M); Gm <- matrix(rbinom(Fn*M, 2, p), Fn, M)
  sibg <- function() t(apply(Gf, 1, transmit)) + t(apply(Gm, 1, transmit))
  mk <- function(Gs) { Z <- as.numeric(Gs %*% betas)
    u <- rnorm(Fn); X <- cfg$alpha_s*Z + cfg$alpha_u*u + rnorm(Fn)        # stratum NOT in X
    w <- rgamma(Fn, 1/cfg$theta_f, scale = cfg$theta_f); cn <- pmin(rexp(Fn, cfg$censoring_rate), cfg$max_follow_up)
    uf <- matrix(runif(Fn*3), Fn, 3); oc <- outcome(X, u, w, strat, cn, uf, gstrat)
    list(Z = Z, X = X, u = u, w = w, cn = cn, uf = uf, time = oc$time, status = oc$status) }
  list(s1 = mk(sibg()), s2 = mk(sibg()), strat = strat, gstrat = gstrat)
}
pop_est <- function(d) {
  ti <- rbind(d$s1$time, d$s2$time); si <- rbind(d$s1$status, d$s2$status)
  Z <- c(d$s1$Z, d$s2$Z); X <- c(d$s1$X, d$s2$X); cc <- subj(ti, si, ti, si, rep(1, length(Z)))
  cov(Z, cc[,1] - cc[,2]) / cov(Z, X) / (length(Z) - 1)                   # /(n-1): net-win probability scale
}
wf_est <- function(d) {
  h <- hpair(d$s1$time, d$s1$status, d$s2$time, d$s2$status)
  dZ <- d$s1$Z - d$s2$Z; dX <- d$s1$X - d$s2$X; sum(h*dZ) / sum(dX*dZ)
}
oracle <- function(d, eps = 0.25) { A <- d$s1; npop <- length(A$X)
  Bp <- outcome(A$X+eps, A$u, A$w, d$strat, A$cn, A$uf, d$gstrat); Bm <- outcome(A$X-eps, A$u, A$w, d$strat, A$cn, A$uf, d$gstrat)
  v <- mrwin_fast_pair_win_loss(Bp$time, Bp$status, Bm$time, Bm$status)
  (as.numeric(v[["wins"]]) - as.numeric(v[["losses"]])) / (as.numeric(npop)^2) / (2*eps)
}
cat("within-family vs population win-MR under population stratification (R=16, F=4000)\n")
res <- do.call(rbind, parallel::mclapply(1:16, function(r) {
  d <- gen_fam(4000L, 1000L + r); c(orc = oracle(d), pop = pop_est(d), wf = wf_est(d)) }, mc.cores = ncores))
cat(sprintf("oracle causal gradient = %.3f (sd %.3f)\n", mean(res[,"orc"]), sd(res[,"orc"])))
cat(sprintf("population win-MR      = %.3f (sd %.3f)  bias = %+.3f\n", mean(res[,"pop"]), sd(res[,"pop"]), mean(res[,"pop"]-res[,"orc"])))
cat(sprintf("within-family win-MR   = %.3f (sd %.3f)  bias = %+.3f\n", mean(res[,"wf"]), sd(res[,"wf"]), mean(res[,"wf"]-res[,"orc"])))

#!/usr/bin/env Rscript
# ============================================================================
# P1(a) partial result: does the transport-velocity/decidability structure
# LIFT tier-by-tier to the K=3 lexicographic (hierarchical, cascading-censored,
# shared-frailty) composite the mrwin programme actually uses?
#
# EXACT STRUCTURAL IDENTITY (holds by construction, no assumptions):
#   h(A,B) = sum_{k=1}^K  1{R_k(A,B)} * sign_k(A,B),   R_k := {tiers 1..k-1 tied}
#   => NB(x+eps,x-eps) = sum_k E[ 1{R_k} * sign_k(A,B) ]  =: sum_k contrib_k
# This telescoping decomposition is a tautology (unpacks the cascade rule); the
# OPEN QUESTION is whether each term further factors as
#   contrib_k ~ (conditional decidability of tier k | reached tier k) x (tier-k
#   causal velocity), i.e. whether the K=1 result (scalar script: 4*INT f^2 v;
#   this session's new PH+censoring closed form: -alpha*D(p)) lifts tier-wise
#   once we condition on reaching tier k -- or whether cross-tier correlation
#   from the SHARED FRAILTY breaks the clean product form (a selection effect:
#   "reaching tier k" is not independent of the tier-k causal gradient, because
#   frailty correlates outcomes across tiers).
#
# TEST: track, for the do(X+-eps) oracle comparison, WHICH tier resolves each
# pair (exact O(n^2) tier cascade on a moderate subsample), and separately
# measure each tier's OWN conditional decidability at the BASE population.
# Compare: does contrib_k / [R_k * D_k^cond] recover a roughly K-INDEPENDENT
# constant across tiers and across a censoring sweep (the K=1 pattern), or
# does it drift with tier index / censoring (evidence of an entangled,
# frailty-driven correction term the K=1 theory does not capture)?
#
# Ground truth: mrwin's own frailty simulator (Gamma frailty, cascading
# censoring); public sim + package only.
# ============================================================================
suppressMessages(library(mrwin))
mk_cfg <- function(maxfu) { c <- mrwin_config(m_snps = 40L); c$max_follow_up <- maxfu; c }
gen_out <- function(X, u, w, cn, uf, cfg) {
  n <- length(X); et <- matrix(NA_real_, n, 3); lw <- log(w)
  for (j in 1:3) { lp <- cfg$alpha_x[j]*X + cfg$nu_u[j]*u + lw
    et[, j] <- (-log(uf[, j]) / (cfg$baseline_haz[j]*exp(lp)))^(1/cfg$shape_weibull) }
  td <- et[,1]; th <- et[,2]; tr <- et[,3]
  list(time = cbind(pmin(td,cn), pmin(th,cn,td), pmin(tr,cn,td)),
       status = cbind(as.integer(td<=cn), as.integer(th<=cn & th<=td), as.integer(tr<=cn & tr<=td)))
}
gen_pop <- function(n, seed, mafs, betas, alpha_s, cens, cfg) {
  set.seed(seed)
  g <- scale(sapply(mafs, function(p) rbinom(n, 2, p))); storage.mode(g) <- "double"
  s <- drop(g %*% betas); u <- rnorm(n); w <- rgamma(n, 1/cfg$theta_f, scale = cfg$theta_f)
  cn <- pmin(rexp(n, cens), cfg$max_follow_up); uf <- matrix(runif(n*3), n, 3)
  X <- alpha_s*s + cfg$alpha_u*u + rnorm(n); oc <- gen_out(X, u, w, cn, uf, cfg)
  list(time = oc$time, status = oc$status, X = X, u = u, w = w, cn = cn, uf = uf, cfg = cfg)
}

# exact O(n1*n2) tier cascade: returns, PER TIER, (n_win, n_lose, n_reached) counts
# over the full pairwise comparison of group 1 (n1) vs group 2 (n2).
tier_cascade <- function(time1, status1, time2, status2) {
  n1 <- nrow(time1); n2 <- nrow(time2); K <- ncol(time1)
  reached <- matrix(TRUE, n1, n2)  # tiers 1..k-1 all tied so far
  out <- vector("list", K)
  for (k in seq_len(K)) {
    t1k <- time1[, k]; s1k <- status1[, k]; t2k <- time2[, k]; s2k <- status2[, k]
    W <- outer(t1k, t2k, ">") & matrix(s2k == 1, n1, n2, byrow = TRUE) & reached
    L <- outer(t1k, t2k, "<") & matrix(s1k == 1, n1, n2) & reached & !W
    out[[k]] <- c(n_reached = sum(reached), n_win = sum(W), n_lose = sum(L))
    reached <- reached & !W & !L
  }
  out
}

cfgH <- mk_cfg(maxfu = 300)
set.seed(11); mafs <- runif(40, cfgH$maf_low, cfgH$maf_high); betas <- rnorm(40, 0, cfgH$sigma_beta)
N <- 1400L; alpha_s <- 0.4
cens_grid <- c(0.02, 0.15, 0.6, 2.0)
eps <- 0.25

cat("== P1(a) per-tier decomposition: contrib_k vs R_k * D_k^cond, across tiers and censoring ==\n")
cat(sprintf("N=%d (exact O(N^2) tier cascade), K=3, alpha_s=%.2f\n\n", N, alpha_s))
cat(sprintf("%-8s %-4s %-10s %-12s %-14s %-16s %-16s\n",
            "cens", "k", "R_k(reach)", "D_k^cond", "contrib_k", "contrib_k/(R_k*D_k)", "cum.contrib/g_total"))
allrows <- list()
for (cr in cens_grid) {
  d <- gen_pop(N, 5000L + round(cr*100), mafs, betas, alpha_s, cr, cfgH)
  # base-population per-tier conditional decidability (self x self cascade, remove diagonal effect by using n(n-1) denom approx via n^2 since n large)
  base_casc <- tier_cascade(d$time, d$status, d$time, d$status)
  # oracle do(X+-eps) per-tier cascade
  Bp <- gen_out(d$X + eps, d$u, d$w, d$cn, d$uf, cfgH); Bm <- gen_out(d$X - eps, d$u, d$w, d$cn, d$uf, cfgH)
  or_casc <- tier_cascade(Bp$time, Bp$status, Bm$time, Bm$status)
  tot <- as.numeric(N)^2
  g_total <- sum(sapply(or_casc, function(o) o["n_win"] - o["n_lose"])) / tot / (2 * eps)
  cum <- 0
  for (k in seq_len(3)) {
    Rk <- base_casc[[k]]["n_reached"] / tot
    Dk_cond <- (base_casc[[k]]["n_win"] + base_casc[[k]]["n_lose"]) / base_casc[[k]]["n_reached"]
    contrib_k <- (or_casc[[k]]["n_win"] - or_casc[[k]]["n_lose"]) / tot / (2 * eps)
    cum <- cum + contrib_k
    ratio <- contrib_k / (Rk * Dk_cond)
    cat(sprintf("%-8.2f %-4d %-10.4f %-12.4f %-14.5f %-16.5f %-16.4f\n",
                cr, k, Rk, Dk_cond, contrib_k, ratio, cum / g_total))
    allrows[[length(allrows) + 1]] <- c(cens = cr, k = k, Rk = unname(Rk), Dk = unname(Dk_cond),
                                         contrib = unname(contrib_k), ratio = unname(ratio))
  }
  cat(sprintf("         total g (all tiers) = %.5f  (check: sum of contrib_k above)\n\n", g_total))
}
dat <- do.call(rbind, allrows)
cat("Per-tier ratio contrib_k/(R_k*D_k^cond) summary across the whole grid (tiers x censoring):\n")
for (k in 1:3) {
  sub <- dat[dat[, "k"] == k, "ratio"]
  cat(sprintf("  tier %d: mean=%.4f  sd=%.4f  CV=%.3f  (n=%d)\n", k, mean(sub), sd(sub), sd(sub)/mean(sub), length(sub)))
}
cat(sprintf("\nAcross ALL tiers pooled: CV=%.3f (compare: aggregate whole-composite CV from the hierarchical\n", sd(dat[,"ratio"])/mean(dat[,"ratio"])))
cat("script's g/D was already shown far more stable than raw g -- this checks if the SAME holds tier-by-tier)\n")
cat("DONE-P1A-TIER\n")

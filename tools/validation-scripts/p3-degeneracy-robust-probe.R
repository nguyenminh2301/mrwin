#!/usr/bin/env Rscript
# P3 spec section 6 — the degeneracy-robust AR statistic, isolated on a synthetic
# degree-2 U-statistic with EXACTLY controllable first-projection variance, to
# separate two questions: (a) is the unified limit-law correct? (b) can the
# boundary nuisance c^2 = n*zeta1 be estimated?
#
# Kernel  phi(Wi,Wj) = eps*(Ai+Aj) + Ai*Aj,  Ai ~ N(0,1)  =>  E[phi]=0 (H0 true),
# first projection g=eps*A (zeta1 = eps^2), degenerate part Ai*Aj (rank-1, operator
# eigenvalue 1). Unified limit:  n*U_n -> R(c) = 2c*Z0 + sum_k lam_k (Z_k^2 - 1),
# c^2 = n*zeta1, lam_k the second-projection eigenvalues.
#
# Three H0 tests (size target 0.05):
#   naive   : n*U_n^2 / (4*zeta1_hat)  vs  chi2_1            (standard AR)
#   oracle  : (n*U_n)^2 vs quantile(R(c_true, lam_true)^2)   (TRUE c, lam)
#   plug-in : same with c_hat (bias-corrected) and lam_hat (empirical spectrum)
#
# FINDING (dev/p3-weak-iv-robust-winmr.md section 6): the unified limit is CORRECT
# (oracle holds ~0.05 across the strong->degenerate continuum); naive AR over-rejects
# badly under degeneracy (size ~0.18-0.20); the plug-in is UNSTABLE because c^2=n*zeta1
# is not consistently estimable at the boundary (zeta1 and its estimation bias are
# both O(1/n)). The fix is robust-to-non-estimable-nuisance inference (Andrews-Cheng
# sup-over-CI / least-favorable, or a degenerate-U bootstrap), not plug-in calibration.
ncores <- max(1L, min(4L, parallel::detectCores()))

one <- function(A, eps, B = 6000, r = 6) {
  n <- length(A)
  P <- outer(A, A) + eps*(outer(A, rep(1, n)) + outer(rep(1, n), A)); diag(P) <- 0
  Un <- sum(P[upper.tri(P)]) / choose(n, 2)
  gi <- rowSums(P) / (n - 1); z1 <- var(gi)
  Pt <- P - outer(gi, rep(1, n)) - outer(rep(1, n), gi) + Un; diag(Pt) <- 0  # doubly-centered
  d2 <- sum(Pt^2) / (n*(n - 1)); z1bc <- max(z1 - d2/(n - 1), 0)             # bias-corrected zeta1
  ev <- eigen(Pt, symmetric = TRUE, only.values = TRUE)$values / n          # operator eigenvalues
  lam <- ev[order(abs(ev), decreasing = TRUE)][1:r]
  nUn2 <- (n*Un)^2
  simq <- function(cc, ll) { Z0 <- rnorm(B); Zk <- matrix(rnorm(B*length(ll)), B, length(ll))
    as.numeric(quantile((2*cc*Z0 + as.numeric((Zk^2 - 1) %*% ll))^2, 0.95)) }
  c(naive  = as.numeric(n*Un^2/(4*z1) > 3.841),
    oracle = as.numeric(nUn2 > simq(sqrt(n)*eps, 1)),
    plugin = as.numeric(nUn2 > simq(sqrt(n*z1bc), lam)))
}

cat("H0 size (target 0.05), n=300, R=1500\n")
cat(sprintf("%-7s %-9s %-9s %-9s %-9s\n", "eps", "zeta1", "naive", "oracle", "plug-in"))
for (eps in c(1.0, 0.3, 0.1, 0.05, 0.02, 0.0)) {
  res <- do.call(rbind, parallel::mclapply(1:1500, function(r) {
    set.seed(20000L + as.integer(eps*1000) + r); one(rnorm(300), eps) }, mc.cores = ncores))
  cat(sprintf("%-7.2f %-9.4f %-9.3f %-9.3f %-9.3f\n",
              eps, eps^2, mean(res[, "naive"]), mean(res[, "oracle"]), mean(res[, "plugin"])))
}

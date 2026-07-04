#!/usr/bin/env Rscript
# ============================================================================
# CROSS-CUTTING PROBE #4: closed-form Var(w) as a function of D, K=1, Type-I
# (fixed administrative) censoring, T ~ F continuous.
#
# DERIVATION (by hand; see dev/findings-transport-collapsibility.md for the
# write-up). Using the package's exact win rule (i beats j iff delta_j=1 &
# X_i>X_j; i loses iff delta_i=1 & X_j>X_i), with p := P(T<=c) = F(c):
#   D      = p(2-p)                    (decidability: P(a random pair decided))
#   Var(w) = p(p^2 - 3p + 3) / 3       (variance of the per-subject win score)
# Both are DISTRIBUTION-FREE (depend on F only through p, not on F's shape) --
# a consequence of w(o) and the decided-indicator being invertible functions of
# F(T) ~ Uniform(0,1) under H0 (no systematic effect assumed here; this is the
# BASE-population functional the degeneracy diagnostics condition on).
# KEY ASYMPTOTIC: as p->0 (heavy censoring, D->0), Var(w)/D -> 1/2 EXACTLY --
# a clean linear limit in the degenerate regime that matters most for Paper 03.
#
# VALIDATION: exact O(n log n) leave-one-out computation of the per-subject win
# score (no Monte-Carlo opponent-subsampling noise): for subject i,
#   L_i = #{j != i : delta_j=1, X_j < X_i}   (i beats j)
#   R_i = #{j != i : X_j > X_i}              (i loses to j, only if delta_i=1)
#   w_i = (L_i - delta_i * R_i) / (n - 1)
# L_i via event times being a.s. distinct (T continuous) -> findInterval on
# sorted event times, corrected for self-count when delta_i=1. R_i via
# rank(ties.method="max"). Public sim only (rexp), no package/confidential data.
# ============================================================================
verify_one <- function(p, n = 200000, seed = 1) {
  set.seed(seed)
  Tt <- rexp(n, 1)
  cc <- -log(1 - p)
  X <- pmin(Tt, cc); delta <- as.integer(Tt <= cc)

  eventX <- sort(X[delta == 1])
  Lraw <- findInterval(X, eventX)          # # events with X_j <= X_i
  L <- Lraw - delta                        # subtract self if delta_i==1 (event times a.s. distinct)
  rank_le <- rank(X, ties.method = "max")  # # j with X_j <= X_i (ties only among delta=0 at X=c)
  R <- n - rank_le                         # # j with X_j > X_i, self-excluded automatically

  w <- (L - delta * R) / (n - 1)
  a <- (L + delta * R) / (n - 1)           # per-subject "decided" rate

  c(D_hat = mean(a), Var_w_hat = var(w), Ew_hat = mean(w),
    D_theory = p * (2 - p), Var_w_theory = p * (p^2 - 3 * p + 3) / 3)
}
cat(sprintf("%-6s %-10s %-10s %-12s %-12s %-10s %-10s\n",
            "p", "D_hat", "D_theory", "Var(w)_hat", "Var(w)_th", "ratio_hD", "E[w]_hat"))
for (p in c(0.02, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.98)) {
  r <- verify_one(p)
  cat(sprintf("%-6.2f %-10.5f %-10.5f %-12.6f %-12.6f %-10.4f %-10.2e\n",
              p, r["D_hat"], r["D_theory"], r["Var_w_hat"], r["Var_w_theory"],
              r["Var_w_hat"]/r["D_hat"], r["Ew_hat"]))
}
cat("\nSmall-p asymptotic check: Var(w)/D -> 1/2 as D->0 ?\n")
for (p in c(0.001, 0.01, 0.05)) {
  r <- verify_one(p, n = 400000)
  cat(sprintf("p=%.4f  D_hat=%.5f D_theory=%.5f  Var(w)_hat=%.6f Var(w)_theory=%.6f  Var(w)/D_hat=%.5f (target 0.5)\n",
              p, r["D_hat"], r["D_theory"], r["Var_w_hat"], r["Var_w_theory"], r["Var_w_hat"]/r["D_hat"]))
}
cat("DONE-VERIFY-VARD\n")

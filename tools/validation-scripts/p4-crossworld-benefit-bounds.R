#!/usr/bin/env Rscript
# ============================================================================
# P4, attempt 2 -- bound the RIGHT estimand.
#
# WHY ATTEMPT 1 BOUNDED THE WRONG THING (p4-kantorovich-partial-id-censored.R).
# That script computed sharp bounds on P(win) over couplings of the TWO ARMS.
# But in win statistics the two arms are INDEPENDENT SUBJECTS BY DESIGN --
# NB(x,x') = P(Y_i(x) > Y_j(x')) with i != j -- so the arm-to-arm coupling is
# known (product), not ambiguous, and NB is POINT-identified from the marginals.
# Bounding it was answering a question nobody asks, which is why the headline
# "width propto (1-D)" conjecture came out inverted and had to be retracted.
#
# THE REAL PARTIAL-IDENTIFICATION PROBLEM (and a point the win-statistics
# literature routinely blurs): the win ratio compares DIFFERENT patients, but is
# habitually read as if it described the SAME patient under two treatments. The
# cross-world "probability of individual benefit"
#       P_ben(eps) := P( Y_i(x+eps)  >  Y_i(x-eps) )      [SAME i]
# is a genuinely different quantity, and it is NOT identified by the marginals:
# it depends on the counterfactual coupling, which no experiment reveals. Its
# sharp identified set is exactly a Frechet-Hoeffding / Kantorovich problem --
# the validated LP machinery from attempt 1, now pointed at the right target.
#
# TWO COUPLING REGIMES (reformulation #3 from the attempt-1 write-up):
#   (a) UNCONSTRAINED: couple the OBSERVED pairs (X_A, delta_A) ~ (X_B, delta_B)
#       freely. The naive Frechet bound.
#   (b) SHARED-CENSORING: couple the LATENT event times (T_A, T_B) and derive
#       BOTH observations from ONE administrative cutoff c -- physically correct
#       for the same patient under the same follow-up. This is a real constraint
#       (delta_A, delta_B are then determined by the T-coupling), so (b) must be
#       nested inside (a). Quantifying how much it tightens is the point.
#
# PREDICTIONS TESTED (not assumed):
#   P1. The win statistic's own NB (product coupling) lies strictly inside both
#       identified sets -- i.e. the win ratio is NOT an estimate of P_ben.
#   P2. Under IDENTICAL marginals (no causal effect, alpha=0) the cross-world
#       set is maximally wide -- the classical "P(benefit) is unidentified"
#       result -- giving a hard validation check with a known answer.
#   P3. Whether decidability D controls the width ON THE RIGHT SCALE, i.e. after
#       normalising by the decided mass (reformulation #1), rather than the raw
#       width that misled attempt 1.
#
# Model matches Claim 6 (§7) exactly so the results compose with it:
# T_x ~ Exponential(lambda(x) = lambda0*exp(alpha*x)), Type-I censoring at c,
# package win rule (A wins iff delta_B=1 & X_A>X_B; loses iff delta_A=1 & X_B>X_A;
# both censored = tie). Public simulation + LP only; no package data.
# ============================================================================
suppressMessages(library(lpSolve))

# --- exact win/lose indicator matrices under the package rule ----------------
# tA, tB: latent event times (vectors of equal-mass quantile points)
# c: shared administrative cutoff
wl_matrices <- function(tA, tB, c) {
  XA <- pmin(tA, c); dA <- as.integer(tA <= c)
  XB <- pmin(tB, c); dB <- as.integer(tB <= c)
  nA <- length(tA); nB <- length(tB)
  W <- outer(XA, XB, ">") & matrix(dB == 1L, nA, nB, byrow = TRUE)   # A beats B
  L <- outer(XA, XB, "<") & matrix(dA == 1L, nA, nB, byrow = FALSE)  # A loses to B
  list(W = W * 1, L = L * 1)
}

# --- LP over couplings with fixed (uniform, equal-mass) marginals ------------
# cost: nA x nB matrix; returns sharp min and max of E_coupling[cost]
lp_range <- function(cost, nA, nB) {
  massA <- rep(1 / nA, nA); massB <- rep(1 / nB, nB)
  smax <- lp.transport(-cost, "min", rep("=", nA), massA, rep("=", nB), massB, integers = NULL)
  smin <- lp.transport(cost, "min", rep("=", nA), massA, rep("=", nB), massB, integers = NULL)
  if (smax$status != 0 || smin$status != 0) stop("LP failed: status ", smax$status, "/", smin$status)
  c(min = sum(smin$solution * cost), max = sum(smax$solution * cost))
}

# --- regime (b): couple LATENT times, shared cutoff c -----------------------
crossworld_shared_c <- function(lamA, lamB, c, nbin = 140) {
  u <- (seq_len(nbin) - 0.5) / nbin
  tA <- qexp(u, lamA); tB <- qexp(u, lamB)          # equal-mass quantile points
  m <- wl_matrices(tA, tB, c)
  NBcost <- m$W - m$L
  rng_nb <- unname(lp_range(NBcost, nbin, nbin))   # unname: downstream c(...) would mangle names
  rng_win <- unname(lp_range(m$W, nbin, nbin))
  # reference couplings
  indep_nb <- mean(NBcost)                          # product coupling = the WIN STATISTIC's NB
  indep_win <- mean(m$W)
  comon_nb <- mean(diag(NBcost))                    # rank-invariance (comonotonic)
  D <- mean(m$W + m$L)                              # decidability under the product coupling
  list(nb_min = rng_nb[1], nb_max = rng_nb[2],
       win_min = rng_win[1], win_max = rng_win[2],
       indep_nb = indep_nb, indep_win = indep_win, comon_nb = comon_nb, D = D)
}

# --- regime (a): couple OBSERVED (X, delta) pairs freely --------------------
# marginals of the observed pair: continuum of event times on [0,c) plus an atom
# at c. Coupling is unconstrained, so delta_A and delta_B need not agree with a
# single shared cutoff -- the naive (looser) Frechet bound.
crossworld_unconstrained <- function(lamA, lamB, c, ncont = 130) {
  pA <- pexp(c, lamA); pB <- pexp(c, lamB)
  uA <- (seq_len(ncont) - 0.5) / ncont * pA
  uB <- (seq_len(ncont) - 0.5) / ncont * pB
  XA <- c(qexp(uA, lamA), c); dA <- c(rep(1L, ncont), 0L); mA <- c(rep(pA / ncont, ncont), 1 - pA)
  XB <- c(qexp(uB, lamB), c); dB <- c(rep(1L, ncont), 0L); mB <- c(rep(pB / ncont, ncont), 1 - pB)
  nA <- length(XA); nB <- length(XB)
  W <- (outer(XA, XB, ">") & matrix(dB == 1L, nA, nB, byrow = TRUE)) * 1
  L <- (outer(XA, XB, "<") & matrix(dA == 1L, nA, nB, byrow = FALSE)) * 1
  NBcost <- W - L
  smax <- lp.transport(-NBcost, "min", rep("=", nA), mA, rep("=", nB), mB, integers = NULL)
  smin <- lp.transport(NBcost, "min", rep("=", nA), mA, rep("=", nB), mB, integers = NULL)
  if (smax$status != 0 || smin$status != 0) stop("LP failed (unconstrained)")
  c(nb_min = sum(smin$solution * NBcost), nb_max = sum(smax$solution * NBcost))
}

alpha <- 0.5; lambda0 <- 1; x0 <- 0; eps <- 0.4
lamA <- lambda0 * exp(alpha * (x0 + eps))   # higher hazard = WORSE outcome
lamB <- lambda0 * exp(alpha * (x0 - eps))

cat("=== P4 attempt 2: sharp bounds on the CROSS-WORLD probability of benefit ===\n\n")

## ---- validation P2: identical marginals must give a maximally wide set -----
cat("[P2] Validation with a KNOWN answer -- identical marginals (alpha=0, no causal\n")
cat("     effect, no censoring). Theory: P(benefit) is completely unidentified, so the\n")
cat("     sharp NB set should approach [-1, 1] while the win statistic reads exactly 0.\n")
v0 <- crossworld_shared_c(1, 1, c = 1e9, nbin = 140)
cat(sprintf("     sharp NB = [%.4f, %.4f]   win-statistic NB (product) = %.4f   comonotonic = %.4f\n",
            v0$nb_min, v0$nb_max, v0$indep_nb, v0$comon_nb))
cat(sprintf("     (discretisation limit is +/-(n-1)/n = +/-%.4f at n=140; comonotonic = 0 = rank invariance)\n\n", 139 / 140))

## ---- main sweep: does censoring/decidability control the width? ------------
cat("[P1/P3] Censoring sweep, PH model matching Claim 6 (alpha=0.5, eps=0.4)\n")
cat("        'win-stat NB' is the PROGRAMME's estimand (independent subjects);\n")
cat("        [nb_min, nb_max] is the sharp CROSS-WORLD (same-patient) set.\n\n")
cat(sprintf("%-7s %-8s %-9s %-11s %-9s %-9s %-10s %-11s\n",
            "c", "D", "win-stat", "comonotonic", "nb_min", "nb_max", "width", "width/D"))
res <- list()
for (cc in c(0.15, 0.4, 0.8, 1.5, 3.0, 8.0)) {
  v <- crossworld_shared_c(lamA, lamB, cc)
  w <- v$nb_max - v$nb_min
  res[[as.character(cc)]] <- c(D = v$D, width = w, indep = v$indep_nb)
  cat(sprintf("%-7.2f %-8.4f %-9.4f %-11.4f %-9.4f %-9.4f %-10.4f %-11.4f\n",
              cc, v$D, v$indep_nb, v$comon_nb, v$nb_min, v$nb_max, w, w / v$D))
}

cat("\n[P3] Width scaling on the RAW vs the DECIDABILITY-NORMALISED scale:\n")
Dv <- sapply(res, function(r) r["D"]); Wv <- sapply(res, function(r) r["width"])
fit_raw <- lm(log(Wv) ~ log(Dv))
cat(sprintf("     log-log slope, raw width vs D           = %+.3f  (R^2=%.3f)\n",
            coef(fit_raw)[2], summary(fit_raw)$r.squared))
cat(sprintf("     width/D across the sweep                = %s\n",
            paste(sprintf("%.3f", Wv / Dv), collapse = "  ")))
cat(sprintf("     CV of raw width = %.3f ;  CV of width/D = %.3f\n",
            sd(Wv) / mean(Wv), sd(Wv / Dv) / mean(Wv / Dv)))

## ---- reformulation #3: does the SHARED-CUTOFF constraint tighten the set? ---
cat("\n[#3] Does enforcing ONE shared administrative cutoff (physically correct for\n")
cat("     the same patient) tighten the set vs coupling observed (X,delta) freely?\n\n")
cat(sprintf("%-7s %-24s %-24s %-10s\n", "c", "shared-c (constrained)", "unconstrained (naive)", "tightening"))
for (cc in c(0.4, 0.8, 1.5, 3.0)) {
  vs <- crossworld_shared_c(lamA, lamB, cc)
  vu <- crossworld_unconstrained(lamA, lamB, cc)
  ws <- vs$nb_max - vs$nb_min; wu <- vu["nb_max"] - vu["nb_min"]
  cat(sprintf("%-7.2f [%+.4f, %+.4f] w=%.3f  [%+.4f, %+.4f] w=%.3f  %-10.3f\n",
              cc, vs$nb_min, vs$nb_max, ws, vu["nb_min"], vu["nb_max"], wu, 1 - ws / wu))
}
cat("\n(tightening = 1 - w_constrained/w_unconstrained; >0 means the shared-cutoff\n")
cat(" constraint carries real identifying information the naive Frechet bound discards.)\n")

## ---- discretisation convergence: are the sharp bounds stable in nbin? -------
cat("\n[conv] Discretisation check -- the LP bounds must be stable in nbin before\n")
cat("       any of the above is trustworthy (lesson: never trust one grid).\n")
cat(sprintf("%-8s %-11s %-11s %-10s\n", "nbin", "nb_min", "nb_max", "width"))
for (nb in c(60, 100, 140, 200)) {
  v <- crossworld_shared_c(lamA, lamB, 1.5, nbin = nb)
  cat(sprintf("%-8d %-11.5f %-11.5f %-10.5f\n", nb, v$nb_min, v$nb_max, v$nb_max - v$nb_min))
}

## ---- closed forms for the sharp endpoints? --------------------------------
# Observed above: comonotonic NB coincides with the LP MINIMUM at every row.
# Mechanism (analytic): under proportional hazards with lambda_A > lambda_B,
# qexp(u,lamA) < qexp(u,lamB) for EVERY u, so under the comonotonic (rank-
# invariance) coupling arm A's latent time is always the smaller one -- A can
# never win, and A loses exactly when its own event is observed. Hence
#       NB_comonotonic = -P(T_A <= c) = -p_A,
# and if that also attains the LP minimum, the sharp LOWER endpoint has a closed
# form. Tested here (not assumed), together with whether the countermonotonic
# coupling attains the UPPER endpoint.
cat("\n[closed] Do the sharp endpoints have closed forms?\n")
cat(sprintf("%-7s %-11s %-11s %-11s %-11s %-11s\n",
            "c", "nb_min(LP)", "-p_A (pred)", "nb_max(LP)", "counter-mono", "p_B"))
for (cc in c(0.4, 0.8, 1.5, 3.0, 8.0)) {
  nbin <- 200
  u <- (seq_len(nbin) - 0.5) / nbin
  tA <- qexp(u, lamA); tB <- qexp(u, lamB)
  v <- crossworld_shared_c(lamA, lamB, cc, nbin = nbin)
  mc <- wl_matrices(tA, rev(tB), cc)                 # countermonotonic pairing
  counter_nb <- mean(diag(mc$W - mc$L))
  pA <- pexp(cc, lamA); pB <- pexp(cc, lamB)
  cat(sprintf("%-7.2f %-11.5f %-11.5f %-11.5f %-11.5f %-11.5f\n",
              cc, v$nb_min, -pA, v$nb_max, counter_nb, pB))
}
cat("\n(col2 vs col3: is the sharp LOWER endpoint exactly -P(T_A<=c)?\n")
cat(" col4 vs col5: does the countermonotonic coupling attain the UPPER endpoint?)\n")

## ---- INDEPENDENT analytic cross-check of the UPPER endpoint ----------------
# In the UNCENSORED limit there are no ties, so NB = 2*P(A>B) - 1, and the sharp
# sup of P(A>B) over couplings with fixed marginals is a classical copula-theory
# result (Makarov / Frank-Nelsen-Schweizer):
#       sup P(A > B) = min( 1, inf_t [ 1 - F_A(t) + F_B(t) ] ).
# Computing this from the analytic exponential CDFs and comparing to the LP is a
# genuinely INDEPENDENT check of the whole LP pipeline -- a different derivation
# route to the same number (lesson E: validate the target more than one way).
cat("\n[xcheck] Uncensored limit: LP upper endpoint vs the classical Makarov bound\n")
tt <- seq(1e-6, 200, length.out = 2000001)
sup_pab <- min(1, min(1 - pexp(tt, lamA) + pexp(tt, lamB)))
nb_max_analytic <- 2 * sup_pab - 1
for (nb in c(200, 400)) {
  v <- crossworld_shared_c(lamA, lamB, c = 1e9, nbin = nb)
  cat(sprintf("   nbin=%-5d  LP nb_max = %.5f   Makarov analytic = %.5f   abs.diff = %.5f\n",
              nb, v$nb_max, nb_max_analytic, abs(v$nb_max - nb_max_analytic)))
}
cat(sprintf("   (sup P(A>B) = %.5f attained at the minimiser of 1-F_A+F_B)\n", sup_pab))
cat("DONE-P4-CROSSWORLD\n")

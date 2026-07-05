#!/usr/bin/env Rscript
# ============================================================================
# P4 first instance: sharp (assumption-free) partial identification of the win
# probability, given ONLY the two arms' marginal laws (no coupling/rank-
# invariance assumption), via an explicit discretized Kantorovich transportation
# LP -- and a direct test of the "D (decidability) as identification budget"
# conjecture from dev/findings-transport-collapsibility.md §6/P4.
#
# BACKGROUND / A CAUGHT ERROR (worth recording). A first attempt hand-derived
# the bounds via a Monge/submodularity shortcut, claiming the comonotonic
# coupling always minimizes and countermonotonic always maximizes P(Y1>Y2)
# given fixed marginals. Direct verification caught this as WRONG in general:
# (i) 1(y1>y2) is NOT globally submodular (a 4-point configuration
# y2<y1<y2'<y1' violates the submodularity inequality that held in a different
# configuration); (ii) for the two Gaussians used to first probe this
# (mu1=0.5,mu2=-0.5, equal sigma), an exact discretized LP shows the TRUE sharp
# maximum (1.0) is achieved by the comonotonic coupling (not counter-), while
# the true sharp MINIMUM (~0.3875, converging with discretization) is achieved
# by neither simple extremal coupling -- a genuinely different optimal
# transport plan. The hand-derived shortcut is retracted; the LP below is the
# validated replacement (cross-checked on a hand-solvable 2-point case first).
#
# THE CENSORED CASE (new, connects to D). For the K=1 Type-I administrative-
# censoring model of Claim 6 (T~Exponential(rate=lambda), X=min(T,c),
# delta=1(T<=c), p=P(T<=c)), the two-arm win rule is: i beats j iff delta_j=1 &
# X_i>X_j; i loses iff delta_i=1 & X_j>X_i; tie if both censored. This gives a
# discretizable (continuum + one atom at X=c) LP. We compute the sharp
# [min,max] identified set for P(win) as the CENSORING LEVEL p varies (hence D
# varies, D=p(2-p) from Claim 6), holding the causal hazard-ratio effect fixed,
# to test whether the sharp-bound WIDTH scales with 1-D (heavier censoring =
# wider assumption-free identified set = smaller "identification budget"), as
# originally conjectured in the super-plan.
# ============================================================================
suppressMessages(library(lpSolve))

# sharp [min,max] identified P(win) for two Type-I censored exponential arms,
# rate lambdaA (arm 1) vs lambdaB (arm 2), common admin censoring at c.
censored_frechet <- function(lambdaA, lambdaB, c, ncont = 120) {
  pA <- 1 - exp(-lambdaA * c); pB <- 1 - exp(-lambdaB * c)
  # continuum bins: ncont equal-probability-mass points in the EVENT region,
  # i.e. quantiles of the truncated Exponential on [0,c) with total mass p.
  uA <- (seq_len(ncont) - 0.5) / ncont * pA   # cumulative-prob points within [0,pA)
  uB <- (seq_len(ncont) - 0.5) / ncont * pB
  xA_cont <- -log(1 - uA) / lambdaA; xB_cont <- -log(1 - uB) / lambdaB
  # bins: [1..ncont] continuum (delta=1, mass pA/ncont each), [ncont+1] atom at c (delta=0, mass 1-pA)
  XA <- c(xA_cont, c); dA <- c(rep(1L, ncont), 0L); massA <- c(rep(pA / ncont, ncont), 1 - pA)
  XB <- c(xB_cont, c); dB <- c(rep(1L, ncont), 0L); massB <- c(rep(pB / ncont, ncont), 1 - pB)
  nA <- length(XA); nB <- length(XB)
  win <- outer(XA, XB, ">") & matrix(dB == 1L, nA, nB, byrow = TRUE)      # A beats B
  lose <- outer(XA, XB, "<") & matrix(dA == 1L, nA, nB, byrow = FALSE)    # A loses to B
  cost_win <- win * 1
  # lp.transport requires supply/demand as plain vectors (not necessarily integer);
  # since masses are unequal across bins here (continuum vs atom), pass them directly.
  row_signs <- rep("=", nA); col_signs <- rep("=", nB)
  # integers=NULL: masses are fractional probabilities, NOT lp.transport's
  # integer-transport-plan default (integers=1:(nc*nr)), which is infeasible
  # here (caught via a status-code check, not a silent wrong answer).
  sol_max <- lp.transport(-cost_win, "min", row_signs, massA, col_signs, massB, integers = NULL)
  sol_min <- lp.transport(cost_win, "min", row_signs, massA, col_signs, massB, integers = NULL)
  if (sol_max$status != 0 || sol_min$status != 0) stop("LP infeasible/failed: status ", sol_max$status, "/", sol_min$status)
  list(pmax = sum(sol_max$solution * cost_win), pmin = sum(sol_min$solution * cost_win),
       pA = pA, pB = pB)
}

# validated point coupling (comonotonic, shared-uniform, the mrwin-simulator style)
comonotonic_point <- function(lambdaA, lambdaB, c, n = 2000000, seed = 1) {
  set.seed(seed)
  U <- runif(n)
  TA <- -log(1 - U) / lambdaA; TB <- -log(1 - U) / lambdaB  # SAME U -> comonotonic
  XA <- pmin(TA, c); dA <- as.integer(TA <= c); XB <- pmin(TB, c); dB <- as.integer(TB <= c)
  mean((dB == 1L) & (XA > XB))
}

cat("== P4: sharp identified set for P(win), K=1 Type-I censored, vs decidability D ==\n")
cat("(Kantorovich LP, discretized ncont=120 continuum bins + 1 atom per arm)\n\n")
lambdaA <- 1.3; lambdaB <- 1.0  # fixed causal hazard-ratio effect (A worse than B)
cat(sprintf("%-8s %-8s %-10s %-12s %-12s %-12s %-14s\n",
            "c", "p", "D=p(2-p)", "LP min", "LP max", "width", "comonotonic-pt"))
for (cc in c(0.15, 0.4, 0.8, 1.5, 3.0, 6.0)) {
  r <- censored_frechet(lambdaA, lambdaB, cc)
  pt <- comonotonic_point(lambdaA, lambdaB, cc)
  D <- r$pA * (2 - r$pA)   # decidability of the BASE (unshifted) population uses p; here approx via pA for context
  cat(sprintf("%-8.2f %-8.4f %-10.4f %-12.5f %-12.5f %-12.5f %-14.5f\n",
              cc, r$pA, D, r$pmin, r$pmax, r$pmax - r$pmin, pt))
}
cat("\n(comonotonic-pt = the ACTUAL mrwin-simulator-style shared-uniform coupling's\n")
cat(" win probability -- checks it lies inside [LP min, LP max] at every c, and\n")
cat(" tracks HOW this specific, causally-meaningful coupling relates to the\n")
cat(" assumption-free sharp bounds as censoring grows.)\n\n")

cat("== Width scaling check: is width propto (1-D), i.e. the tie/censored mass? ==\n")
widths <- c(); Ds <- c(); ones_minus_D <- c()
for (cc in c(0.05, 0.15, 0.4, 0.8, 1.5, 3.0, 6.0, 12.0)) {
  r <- censored_frechet(lambdaA, lambdaB, cc, ncont = 150)
  D <- r$pA * (2 - r$pA)
  widths <- c(widths, r$pmax - r$pmin); Ds <- c(Ds, D); ones_minus_D <- c(ones_minus_D, 1 - D)
}
fit <- lm(log(widths) ~ log(ones_minus_D))
cat(sprintf("width vs (1-D): log-log slope = %.3f (1.0 would mean width propto (1-D) exactly)\n",
            coef(fit)[2]))
cat(sprintf("R^2 = %.4f\n", summary(fit)$r.squared))
for (i in seq_along(Ds)) cat(sprintf("  D=%.4f  1-D=%.4f  width=%.5f  width/(1-D)=%.4f\n",
                                       Ds[i], ones_minus_D[i], widths[i], widths[i] / ones_minus_D[i]))
cat("DONE-P4-KANTOROVICH\n")

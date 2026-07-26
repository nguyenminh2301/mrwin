#!/usr/bin/env Rscript
# ============================================================================
# THE TRANSPORT-VELOCITY STRUCTURE OF CAUSAL WIN STATISTICS (scalar theory).
#
# New result (this session). For a scalar potential-outcome family {F_x} with
# density f_x, the INTERVENTIONAL net-benefit gradient (the do(X+-eps) oracle the
# whole mrwin programme uses) is EXACTLY a density-overlap-weighted average of the
# optimal-transport velocity field:
#
#     B'(0) := d/deps NB(x+eps, x-eps)|_0  =  4 * INT f_x(y)^2 v_x(y) dy,          (*)
#
# where v_x(y) = -d_x F_x(y) / f_x(y) is the transport (Lagrangian) velocity of the
# quantile point at y, i.e. v_x(F_x^{-1}(u)) = d/dx F_x^{-1}(u) (the OT map's
# generator). NB(x,x') = P(Y_i(x)>Y_j(x')) - P(Y_i(x)<Y_j(x')).
#
# CONSEQUENCES tested here:
#  (A) verify (*) numerically on THREE families incl. a genuinely non-location one
#      (v varies with y) -- decisive, pure numeric-analysis, no Monte Carlo.
#  (B) NON-COLLAPSIBILITY: the win gradient is v squashed by the density overlap
#      c(x):=4*INT f_x^2. Across studies with identical causal effect but different
#      outcome-noise LAW, the win gradient changes; and the correct squash factor is
#      the OVERLAP INT f^2, NOT the variance sigma: Gaussian and Laplace at matched
#      sigma give DIFFERENT win gradients (ratio sqrt(2/pi)=0.7979), both recovered
#      by dividing by their OWN 4*INT f^2. => the collapsible target is
#      tau(x) := B'(0) / (4 INT f_x^2) = f^2-weighted mean transport velocity, = mu'(x)
#      under any location model, INVARIANT to the noise law.
#  (C) finite-sample: a plug-in estimator tau_hat = (win gradient)/(4 INT fhat^2)
#      recovers mu' and is invariant to sigma, with Monte-Carlo SE.
#
# Ground truth: closed-form NB for Gaussian/Laplace (analytic), and high-resolution
# numeric integration of the marginals for the non-location family. Public sim only.
# ============================================================================
options(digits = 6)

## ---- numeric machinery -----------------------------------------------------
# NB(x1,x2) = 2 INT f_{x1}(y) F_{x2}(y) dy - 1, on a fine grid.
NB_num <- function(f1, F2, grid, dy) 2 * sum(f1 * F2) * dy - 1

# B'(0) by symmetric finite difference of eps around a base x, families as closures.
Bprime_fd <- function(fam, x, eps = 1e-3, grid = NULL, dy = NULL) {
  if (is.null(grid)) { grid <- seq(fam$lo, fam$hi, length.out = fam$n); dy <- grid[2] - grid[1] }
  fp <- fam$f(grid, x + eps); Fp <- fam$F(grid, x + eps)
  fm <- fam$f(grid, x - eps); Fm <- fam$F(grid, x - eps)
  # NB(x+eps, x-eps): "high" pop = x+eps, "low" pop = x-eps
  nb_pe <- 2 * sum(fp * Fm) * dy - 1
  nb_me <- 2 * sum(fm * Fp) * dy - 1   # NB(x-eps, x+eps) = -nb_pe by antisymmetry (check)
  (nb_pe - nb_me) / (2 * (2 * eps)) * 2  # = d/deps NB(x+eps,x-eps)|_0 ; see note below
}
# Cleaner: directly B(eps)=NB(x+eps,x-eps); B'(0) ~ (B(eps)-B(-eps))/(2eps)=B(eps)/eps (B odd).
Bprime <- function(fam, x, eps = 1e-3) {
  grid <- seq(fam$lo, fam$hi, length.out = fam$n); dy <- grid[2] - grid[1]
  fp <- fam$f(grid, x + eps); Fm <- fam$F(grid, x - eps)
  Beps <- 2 * sum(fp * Fm) * dy - 1
  Beps / eps
}

# RHS of (*): 4 INT f_x^2 v_x, with v_x(y) = -d_x F_x(y)/f_x(y) via finite diff in x.
rhs_star <- function(fam, x, hx = 1e-4) {
  grid <- seq(fam$lo, fam$hi, length.out = fam$n); dy <- grid[2] - grid[1]
  f0 <- fam$f(grid, x)
  Fx_p <- fam$F(grid, x + hx); Fx_m <- fam$F(grid, x - hx)
  dxF <- (Fx_p - Fx_m) / (2 * hx)
  v <- ifelse(f0 > 1e-12, -dxF / f0, 0)
  4 * sum(f0^2 * v) * dy
}
overlap <- function(fam, x) {           # INT f_x^2
  grid <- seq(fam$lo, fam$hi, length.out = fam$n); dy <- grid[2] - grid[1]
  sum(fam$f(grid, x)^2) * dy
}

## ---- families --------------------------------------------------------------
alpha <- 0.8  # causal slope mu'(x) = alpha (location families)
# Gaussian location: Y(x) = alpha*x + sigma Z
gauss_loc <- function(sigma) list(
  f = function(y, x) dnorm(y, alpha * x, sigma),
  F = function(y, x) pnorm(y, alpha * x, sigma),
  lo = -40, hi = 40, n = 200001, sigma = sigma)
# Laplace location: scale b, sd = b*sqrt(2)
laplace_loc <- function(sigma) { b <- sigma / sqrt(2); list(
  f = function(y, x) dexp(abs(y - alpha * x), rate = 1 / b) / (2 * b) * (2 * b) / (2 * b), # = (1/(2b))e^{-|.|/b}
  F = function(y, x) ifelse(y < alpha * x, 0.5 * exp((y - alpha * x) / b), 1 - 0.5 * exp(-(y - alpha * x) / b)),
  lo = -60, hi = 60, n = 240001, sigma = sigma, b = b) }
# fix laplace density properly
laplace_loc <- function(sigma) { b <- sigma / sqrt(2); list(
  f = function(y, x) (1 / (2 * b)) * exp(-abs(y - alpha * x) / b),
  F = function(y, x) ifelse(y < alpha * x, 0.5 * exp((y - alpha * x) / b), 1 - 0.5 * exp(-(y - alpha * x) / b)),
  lo = -60, hi = 60, n = 240001, sigma = sigma, b = b) }
# Non-location: variance grows with x. Y(x)=alpha*x + sigma0 e^{gamma x} Z. v_x(y)=alpha+gamma(y-alpha x).
gauss_scale <- function(sigma0, gamma) list(
  f = function(y, x) dnorm(y, alpha * x, sigma0 * exp(gamma * x)),
  F = function(y, x) pnorm(y, alpha * x, sigma0 * exp(gamma * x)),
  lo = -60, hi = 60, n = 240001, sigma0 = sigma0, gamma = gamma)

## ---- (A) verify the local structure formula (*) ----------------------------
cat("== (A) local structure: B'(0) vs 4*INT f^2 v  (should MATCH) ==\n")
cat(sprintf("%-26s %-12s %-14s %-14s %-10s\n","family","x","B'(0) [num]","4*INT f^2 v","rel.err"))
checkA <- function(name, fam, x, analytic = NA) {
  lhs <- Bprime(fam, x, eps = 1e-3); rhs <- rhs_star(fam, x)
  cat(sprintf("%-26s %-12.3f %-14.6f %-14.6f %-10.2e", name, x, lhs, rhs, abs(lhs - rhs) / abs(rhs)))
  if (!is.na(analytic)) cat(sprintf("  analytic=%.6f", analytic))
  cat("\n"); invisible(c(lhs = lhs, rhs = rhs))
}
# Gaussian location analytic B'(0) = 2 alpha /(sigma sqrt(pi))
checkA("Gaussian-loc sigma=1",  gauss_loc(1.0),  0.0, 2 * alpha / (1.0 * sqrt(pi)))
checkA("Gaussian-loc sigma=2",  gauss_loc(2.0),  0.0, 2 * alpha / (2.0 * sqrt(pi)))
checkA("Laplace-loc  sigma=1",  laplace_loc(1.0), 0.0)
checkA("Gauss-scale g=0.3 x=0",  gauss_scale(1.0, 0.3), 0.0)   # non-location: v varies with y
checkA("Gauss-scale g=0.3 x=0.5",gauss_scale(1.0, 0.3), 0.5)

## ---- (B) non-collapsibility + the density-overlap correction ---------------
cat("\n== (B) win gradient is NON-collapsible; correct factor is INT f^2 (not sigma) ==\n")
cat(sprintf("%-14s %-8s %-14s %-14s %-16s %-14s\n",
            "family","sigma","win grad B'(0)","4*INT f^2","tau=B'(0)/(4INTf^2)","(target mu'=alpha)"))
for (sig in c(0.5, 1.0, 2.0, 4.0)) {
  for (fm in list(list("Gaussian", gauss_loc(sig)), list("Laplace", laplace_loc(sig)))) {
    fam <- fm[[2]]; b <- Bprime(fam, 0.0); ov4 <- 4 * overlap(fam, 0.0)
    cat(sprintf("%-14s %-8.2f %-14.6f %-14.6f %-16.6f %-14.3f\n",
                fm[[1]], sig, b, ov4, b / ov4, alpha))
  }
}
cat(sprintf("\nAt matched sigma=1: win-grad(Gauss)/win-grad(Laplace) = %.4f  (predicted sqrt(2/pi)=%.4f)\n",
            Bprime(gauss_loc(1), 0) / Bprime(laplace_loc(1), 0), sqrt(2 / pi)))

## ---- (C) finite-sample plug-in estimator of the collapsible effect ---------
cat("\n== (C) finite-sample: raw win-gradient vs overlap-corrected, Monte Carlo ==\n")
# Estimate the interventional win gradient by the paired do(X+-eps) oracle on samples;
# estimate INT f^2 by a plug-in (Gaussian-kernel U-statistic: INT fhat^2 = mean over
# pairs of K_h(Yi-Yj)/... ). We use the exact analytic B'(0) target and check the
# overlap estimate + correction on samples.
intf2_hat <- function(y) {                         # unbiased-ish kernel estimate of INT f^2
  n <- length(y); h <- 1.06 * sd(y) * n^(-1/5)     # Silverman
  # INT fhat^2 = (1/n^2) sum_ij phi_{sqrt2 h}(yi-yj) for Gaussian kernel; use all pairs
  # subsample pairs for speed
  m <- min(n, 1500); idx <- sample(n, m)
  yy <- y[idx]; d <- outer(yy, yy, "-")
  mean(dnorm(d, 0, sqrt(2) * h))
}
win_grad_sample <- function(gen, x, eps = 0.25, n = 20000, seed = 1) {
  set.seed(seed)
  z <- rnorm(n)                                     # shared latent for +/- (paired do)
  yp <- gen(x + eps, z); ym <- gen(x - eps, z)
  # NB between the +eps population and -eps population over all i,j pairs:
  # P(yp_i > ym_j) - P(yp_i < ym_j), via sort-based rank (O(n log n))
  allv <- c(yp, ym); r <- rank(allv, ties.method = "average")
  Rp <- sum(r[1:n]); # sum of ranks of yp among all 2n
  # Mann-Whitney: U = Rp - n(n+1)/2 = #(yp_i > ym_j) (ties .5). NB = (2U/n^2) - 1
  U <- Rp - n * (n + 1) / 2
  nb <- 2 * U / n^2 - 1
  # B(eps)=NB(x+eps,x-eps); B'(0)=B(eps)/eps (B odd). Use B'(0) convention (matches Part A/B).
  list(nb = nb, grad = nb / eps, y = c(yp, ym))
}
gen_gauss <- function(sigma) function(x, z) alpha * x + sigma * z
cat(sprintf("%-10s %-8s %-16s %-16s %-16s\n","family","sigma","raw win grad","INT f^2 hat","corrected tau"))
for (sig in c(0.5, 1.0, 2.0)) {
  gr <- c(); ov <- c(); tc <- c()
  for (s in 1:20) {
    r <- win_grad_sample(gen_gauss(sig), 0, eps = 0.15, n = 12000, seed = 100 + s)
    o <- intf2_hat(r$y)
    gr <- c(gr, r$grad); ov <- c(ov, o); tc <- c(tc, r$grad / (4 * o))
  }
  cat(sprintf("%-10s %-8.2f %-16s %-16s %-16s\n","Gaussian", sig,
              sprintf("%.4f(%.4f)", mean(gr), sd(gr)/sqrt(20)),
              sprintf("%.4f", mean(ov)),
              sprintf("%.4f(%.4f)", mean(tc), sd(tc)/sqrt(20))))
}
cat(sprintf("target mu' = alpha = %.3f  (corrected tau should hit this at every sigma)\n", alpha))
cat("DONE-SCALAR\n")

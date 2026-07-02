# One-sample Anderson-Rubin inference for win-ratio MR (Paper 03, eq. 4-9).
#
# Packages the one-sample pairwise-U-statistic AR test that Paper 03's Table 1
# (AR vs Wald coverage) and Section 8.3 (over-ID pleiotropy) validate, and the
# degeneracy diagnostic of Section 6. Prior to this file these existed only as
# closures inside tools/validation-scripts/p3-weak-iv-robust-probe.R and
# p3-degeneracy-supci-probe.R; the numerical logic here is ported unchanged from
# those validated closures (lesson G: reuse validated machinery, do not re-derive).
#
# Closed-form building blocks (dev/p3-weak-iv-robust-winmr.md eq. 3-6), reusing
# only the shipped per-subject win-loss kernel `mrwin_subject_win_loss_cpp`:
#   W_i   = subject i's unweighted net win count vs everyone (weights_o = 1)
#   rho_i = subject i's Z-weighted net win count vs everyone (weights_o = Z)
#   U_n^h = sum_i Z_i W_i / choose(n,2),   U_n^x = 2 sum_i Z_i(X_i-Xbar) / (n-1)
#   ghat_i(beta) = Z_i w(O_i) - rho_i/(n-1) - beta*(X_i-Xbar)(Z_i-Zbar)
#   zeta1(beta) = Var_i(ghat_i(beta));  AR(beta0) = n*U_n(beta0)^2/(4*zeta1(beta0))
# Both sums are antisymmetric-kernel identities (the same algebra as the trio
# identity of the companion within-family paper): U_n^h, U_n^x reduce to a single
# O(n) combination once the O(n log n / n^2) per-subject aggregates are computed.

.mrwin_ar1_blocks <- function(time, status, X, Z) {
  subj <- getFromNamespace("mrwin_subject_win_loss_cpp", "mrwin")
  n <- length(X)
  cw <- subj(time, status, time, status, rep(1, n)); s_i <- cw[, 1] - cw[, 2]
  cz <- subj(time, status, time, status, Z); r_i <- cz[, 1] - cz[, 2]
  a <- 2 * sum(Z * s_i) / (n * (n - 1))
  b <- 2 * (n * sum(X * Z) - sum(X) * sum(Z)) / (n * (n - 1))
  gh <- (Z * s_i - r_i) / (n - 1)
  SX <- sum(X); SZ <- sum(Z); SXZ <- sum(X * Z)
  gx <- ((n - 1) * X * Z - X * (SZ - Z) - Z * (SX - X) + (SXZ - X * Z)) / (n - 1)
  list(a = a, b = b, c0 = stats::var(gh), c1 = stats::cov(gh, gx), c2 = stats::var(gx),
       n = n, s_i = s_i)
}
.mrwin_ar1_stat <- function(m, beta) {
  m$n * (m$a - m$b * beta)^2 / (4 * (m$c0 - 2 * beta * m$c1 + beta^2 * m$c2))
}

#' One-sample identification-robust (Anderson-Rubin) inference for win-ratio MR
#'
#' Weak-instrument-robust inference for the causal net-benefit gradient from
#' individual-level one-sample data, inverting the pairwise Anderson-Rubin moment
#' `AR(b) = n*(a - b*beta_hat_denom)^2 / (4*zeta1_hat(b))`, a degree-2 U-statistic
#' analogue of the linear-IV AR test, valid under any instrument strength (it
#' never divides by `Cov(Z,X)`). Because both the numerator and the plug-in
#' variance `zeta1_hat(b)` are quadratic in `b`, the confidence set has a closed
#' form (the same quadratic-inversion geometry as the Fieller interval), returned
#' exactly rather than by grid search.
#'
#' A degeneracy diagnostic is also reported: the pairwise moment can become a
#' degenerate U-statistic when its first (Hajek) projection variance
#' `zeta1(beta_hat)` is small, in which case the naive chi-squared critical value
#' over-rejects (see [mrwin_ar_onesample_supci_test()]). This function estimates
#' a percentile-bootstrap confidence interval for the concentration parameter
#' `c^2 = n*zeta1(beta_hat)` (a point estimate alone is unreliable near the
#' boundary) and flags `degenerate = TRUE` when that interval's lower bound is
#' not clearly bounded away from zero.
#'
#' @param time,status N x K matrices of event times and 0/1 status, in priority
#'   order (column 1 = highest priority), as elsewhere in the package.
#' @param X Numeric exposure (length N).
#' @param Z Numeric scalar instrument, e.g. a polygenic score (length N).
#' @param alpha Test level (default 0.05, i.e. a 95% confidence set).
#' @param degeneracy_check Logical; if `TRUE` (default) also compute the
#'   bootstrap degeneracy diagnostic described above.
#' @param boot_reps Number of bootstrap resamples for the degeneracy diagnostic
#'   (default 200); each resample re-runs the O(N log N / N^2) per-subject
#'   kernel, so this is the dominant cost when `degeneracy_check = TRUE`.
#' @param alpha1 One-sided level for the bootstrap CI on `c^2` (default 0.025;
#'   the CI's lower endpoint, not the point estimate, drives the flag).
#' @param seed Optional RNG seed for the bootstrap (reproducibility).
#' @return A list: `gamma` (point estimate = `a/b`), `ci` (length-2 bounds, may be
#'   `c(-Inf, Inf)`), `ci_level`, `bounded`, `empty` (TRUE if no `b` satisfies the
#'   test -- a numerical edge case, not expected in practice), `degenerate`
#'   (logical, `NA` if `degeneracy_check = FALSE`), `c2_hat`, `c2_boot_ci`
#'   (length-2 bootstrap CI for `c^2 = n*zeta1(gamma)`), `n`.
#' @seealso [mrwin_ar_onesample_overid()], [mrwin_ar_onesample_supci_test()],
#'   [mrwin_winmr_ar()] (the two-sample / summary-data analogue)
#' @export
mrwin_ar_onesample <- function(time, status, X, Z, alpha = 0.05,
                                degeneracy_check = TRUE, boot_reps = 200L,
                                alpha1 = 0.025, seed = NULL) {
  time <- as.matrix(time); status <- as.matrix(status)
  X <- as.numeric(X); Z <- as.numeric(Z); n <- length(X)
  if (nrow(time) != n || nrow(status) != n || length(Z) != n) {
    stop("`time`/`status`/`Z` rows must match length(X).", call. = FALSE)
  }
  m <- .mrwin_ar1_blocks(time, status, X, Z)
  q <- stats::qchisq(1 - alpha, 1)
  A <- m$n * m$b^2 - 4 * q * m$c2
  B <- -2 * m$n * m$a * m$b + 8 * q * m$c1
  C <- m$n * m$a^2 - 4 * q * m$c0
  disc <- B^2 - 4 * A * C
  gamma <- m$a / m$b
  empty <- FALSE; bounded <- FALSE; ci <- c(-Inf, Inf)
  if (abs(A) < .Machine$double.eps^0.5) {
    if (abs(B) < .Machine$double.eps^0.5) {
      empty <- C > 0                       # constant sign, no root
    } else {
      root <- -C / B                       # linear: one-sided half-line
      ci <- if (B > 0) c(-Inf, root) else c(root, Inf)
    }
  } else if (disc < 0) {
    if (A > 0) empty <- TRUE else ci <- c(-Inf, Inf)
  } else {
    r1 <- (-B - sqrt(disc)) / (2 * A); r2 <- (-B + sqrt(disc)) / (2 * A)
    lo <- min(r1, r2); hi <- max(r1, r2)
    if (A > 0) { ci <- c(lo, hi); bounded <- TRUE } else { ci <- c(-Inf, Inf) }
    # A<0 with real roots: the *rejected* region is (lo,hi); the accepted set is
    # the two-ray complement, i.e. unbounded -- report as unbounded, honestly.
  }
  c2_hat <- max(m$n * (m$c0 - 2 * gamma * m$c1 + gamma^2 * m$c2), 0)
  degenerate <- NA; c2_boot_ci <- c(NA_real_, NA_real_)
  if (isTRUE(degeneracy_check)) {
    if (!is.null(seed)) {
      old <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
      set.seed(seed)
      on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
    }
    boots <- vapply(seq_len(boot_reps), function(b) {
      idx <- sample.int(n, n, replace = TRUE)
      mb <- tryCatch(.mrwin_ar1_blocks(time[idx, , drop = FALSE], status[idx, , drop = FALSE], X[idx], Z[idx]),
                     error = function(e) NULL)
      if (is.null(mb)) return(NA_real_)
      max(mb$n * (mb$c0 - 2 * gamma * mb$c1 + gamma^2 * mb$c2), 0)
    }, numeric(1))
    boots <- boots[is.finite(boots)]
    if (length(boots) >= max(10L, boot_reps / 4)) {
      c2_boot_ci <- stats::quantile(boots, c(alpha1, 1 - alpha1), na.rm = TRUE, names = FALSE)
      degenerate <- c2_boot_ci[1] < 1  # lower CI bound not clearly bounded away from the c^2=0 boundary
    }
  }
  list(gamma = gamma, ci = ci, ci_level = 1 - alpha, bounded = bounded, empty = empty,
       degenerate = degenerate, c2_hat = c2_hat, c2_boot_ci = c2_boot_ci, n = n)
}

#' One-sample over-identification (pleiotropy) test for multi-variant win-ratio MR
#'
#' The vector-instrument (e.g. multi-SNP) analogue of [mrwin_ar_onesample()]'s
#' point estimate, giving a weak-instrument-robust over-identification test: the
#' minimized Anderson-Rubin statistic `min_b AR(b)`, computed with an efficient
#' (two-step GMM) weight matrix, which is `chi^2_{L-1}` under a valid instrument
#' (a degree-2-U-statistic analogue of the Sargan/Hansen-J test). Rejection
#' signals a violated exclusion restriction (pleiotropy); the test is, by
#' construction, blind to pleiotropy proportional to instrument strength
#' (InSIDE-violating, absorbed into the point estimate) and powered against
#' pleiotropy independent of instrument strength (InSIDE-satisfying).
#'
#' @param time,status N x K matrices of event times and 0/1 status.
#' @param X Numeric exposure (length N).
#' @param G N x L instrument matrix (e.g. L SNP dosages).
#' @return A list: `Q` (the minimized AR statistic), `Q_df` (`L-1`), `Q_p`,
#'   `gamma0` (the preliminary point estimate used to form the weight matrix),
#'   `n_snp`.
#' @section Many-instrument warning:
#' This construction estimates an `L x L` weight matrix from `n` observations
#' and is validated only for `L` small relative to `n`. A stress simulation
#' (`tools/validation-scripts/p3-manyweak-instrument-probe.R`) shows type-I
#' error under a valid instrument is calibrated at `L/n = 0.016` (0.02) and
#' `L/n = 0.04` (0.09, mildly inflated), but breaks down sharply beyond that:
#' 0.71 at `L/n = 0.2` and 1.00 (always rejects) at `L/n = 0.5`. This mirrors
#' the many-weak-instrument problem known for the linear-IV analogue of this
#' test (Kleibergen 2002, Econometrica), for which a correction exists in the
#' linear case; no such correction is implemented here. This function warns
#' when `L/n > 0.02` and errors when `L/n > 0.15`, thresholds set from the
#' simulation above, not a formal bound.
#' @seealso [mrwin_ar_onesample()], [mrwin_winmr_ar()] (the two-sample analogue,
#'   `Q = min_b AR(b)` there is the same construction on summary statistics)
#' @export
mrwin_ar_onesample_overid <- function(time, status, X, G) {
  time <- as.matrix(time); status <- as.matrix(status); X <- as.numeric(X)
  G <- as.matrix(G); n <- length(X); L <- ncol(G)
  if (nrow(time) != n || nrow(status) != n || nrow(G) != n) {
    stop("`time`/`status`/`G` rows must match length(X).", call. = FALSE)
  }
  if (L < 2L) stop("`G` must have at least 2 columns for an over-ID test.", call. = FALSE)
  if (L / n > 0.15) {
    stop(sprintf(
      "L/n = %.3f exceeds 0.15: this test's weight-matrix estimation is known to badly over-reject in this regime (simulated type-I ~1.0 at L/n=0.5); refusing to return a result. See ?mrwin_ar_onesample_overid.",
      L / n), call. = FALSE)
  }
  if (L / n > 0.02) {
    warning(sprintf(
      "L/n = %.3f exceeds 0.02: this test is validated as calibrated only for small L/n (type-I 0.02 at L/n=0.016; already mildly inflated to 0.09 at L/n=0.04); treat this result as approximate. See ?mrwin_ar_onesample_overid.",
      L / n), call. = FALSE)
  }
  subj <- getFromNamespace("mrwin_subject_win_loss_cpp", "mrwin")
  cw <- subj(time, status, time, status, rep(1, n)); s_i <- cw[, 1] - cw[, 2]
  Uh <- 2 * as.numeric(crossprod(G, s_i)) / (n * (n - 1))
  Ux <- 2 * (n * as.numeric(crossprod(G, X)) - colSums(G) * sum(X)) / (n * (n - 1))
  GH <- matrix(0, n, L)
  for (l in seq_len(L)) {
    cz <- subj(time, status, time, status, G[, l]); r <- cz[, 1] - cz[, 2]
    GH[, l] <- (G[, l] * s_i - r) / (n - 1)
  }
  SX <- sum(X); GX <- matrix(0, n, L)
  for (l in seq_len(L)) {
    Zl <- G[, l]; SZ <- sum(Zl); SXZ <- sum(X * Zl)
    GX[, l] <- ((n - 1) * X * Zl - X * (SZ - Zl) - Zl * (SX - X) + (SXZ - X * Zl)) / (n - 1)
  }
  gamma0 <- sum(Ux * Uh) / sum(Ux * Ux)
  Gp <- GH - gamma0 * GX
  Sig <- 4 * stats::cov(Gp)
  Si <- tryCatch(solve(Sig + diag(1e-8 * mean(diag(Sig)), L)), error = function(e) NULL)
  if (is.null(Si)) stop("Weight matrix is numerically singular; too few informative variants.", call. = FALSE)
  num <- as.numeric(t(Ux) %*% Si %*% Uh); den <- as.numeric(t(Ux) %*% Si %*% Ux)
  Q <- n * (as.numeric(t(Uh) %*% Si %*% Uh) - num^2 / den)
  list(Q = Q, Q_df = L - 1L, Q_p = stats::pchisq(Q, L - 1L, lower.tail = FALSE),
       gamma0 = gamma0, n_snp = L)
}

#' Degeneracy-robust Anderson-Rubin point-null test (least-favorable / Andrews-Cheng)
#'
#' When [mrwin_ar_onesample()] flags `degenerate = TRUE`, the basic chi-squared
#' Anderson-Rubin test can over-reject: the pairwise U-statistic moment can
#' degenerate (its first-order/Hajek variance vanishes) and acquire a
#' non-Gaussian limit, `n*U_n -> N(0,4c^2) + sum_k lambda_k(Z_k^2-1)`, so a single
#' chi-squared critical value is no longer valid. This function implements the
#' least-favorable (Andrews-Cheng) remedy validated in
#' `tools/validation-scripts/p3-degeneracy-supci-probe.R`: because the limiting
#' quantile is monotone increasing in `c^2`, the critical value is taken at the
#' upper endpoint of a bootstrap confidence interval for `c^2`, combined with the
#' estimated spectrum of the degenerate part via Monte Carlo simulation of the
#' limit law. This is **uniformly valid but conservative near the degenerate
#' boundary** (a genuine validity-power frontier, not a tuning artifact -- see
#' Paper 03 Section 6): `c^2` is weakly identified there, so the valid upper
#' limit is necessarily wide.
#'
#' Unlike [mrwin_ar_onesample()], this function tests a single specified null
#' value `beta0` (as validated) rather than returning a confidence interval: the
#' O(N^2) pairwise kernel matrix and its eigendecomposition needed for the
#' spectrum are too expensive to recompute at every point of a confidence-set
#' grid scan. It is intended for the flagged, extreme-tail case, not as the
#' default inference (the ordinary [mrwin_ar_onesample()] CI is exact and cheap
#' in the non-degenerate regime that covers ordinary practice).
#'
#' @param time,status N x K matrices of event times and 0/1 status.
#' @param X Numeric exposure (length N).
#' @param Z Numeric scalar instrument (length N).
#' @param beta0 The null value of the causal gradient to test.
#' @param alpha Overall test level (default 0.05); split `alpha = alpha1 + alpha2`
#'   (Bonferroni) between the bootstrap CI for `c^2` and the simulated critical
#'   value.
#' @param alpha1 One-sided level for the bootstrap CI on `c^2` (default 0.025).
#' @param boot_reps Bootstrap resamples for the `c^2` CI (default 300).
#' @param mc_reps Monte Carlo draws for the limiting-law quantile (default 4000).
#' @param r_eigen Number of leading eigenvalues (by magnitude) of the degenerate
#'   part retained (default 6; the tail is dropped, a documented approximation).
#' @param seed Optional RNG seed (reproducibility).
#' @return A list: `stat` (`(n*U_n(beta0))^2`), `crit` (the least-favorable
#'   critical value), `reject` (logical), `c2_hi` (the bootstrap upper bound
#'   driving the critical value), `n`. Cost is O(N^2) in memory (the pairwise
#'   kernel matrix at `beta0`) and O(N^3) for the eigendecomposition; a
#'   subquadratic spectrum estimator is an open question (Paper 03 Section 9),
#'   so this function is recommended for moderate N (a few thousand), not
#'   biobank-scale data.
#' @seealso [mrwin_ar_onesample()]
#' @export
mrwin_ar_onesample_supci_test <- function(time, status, X, Z, beta0, alpha = 0.05,
                                          alpha1 = 0.025, boot_reps = 300L,
                                          mc_reps = 4000L, r_eigen = 6L, seed = NULL) {
  time <- as.matrix(time); status <- as.matrix(status)
  X <- as.numeric(X); Z <- as.numeric(Z); n <- length(X)
  if (nrow(time) != n || nrow(status) != n || length(Z) != n) {
    stop("`time`/`status`/`Z` rows must match length(X).", call. = FALSE)
  }
  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
    set.seed(seed)
    on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
  }
  hpair <- function(t1, s1, t2, s2) {
    K <- ncol(t1); res <- integer(nrow(t1))
    for (k in seq_len(K)) {
      dec <- res != 0
      iw <- !dec & s2[, k] == 1 & t1[, k] > t2[, k]
      il <- !dec & !iw & s1[, k] == 1 & t2[, k] > t1[, k]
      res[iw] <- 1L; res[il] <- -1L
    }
    res
  }
  Phi <- matrix(0, n, n)
  for (i in seq_len(n - 1L)) {
    j <- (i + 1L):n
    hij <- hpair(time[rep(i, length(j)), , drop = FALSE], status[rep(i, length(j)), , drop = FALSE],
                 time[j, , drop = FALSE], status[j, , drop = FALSE])
    val <- (hij - beta0 * (X[i] - X[j])) * (Z[i] - Z[j])
    Phi[i, j] <- val; Phi[j, i] <- val
  }
  Un <- sum(Phi[upper.tri(Phi)]) / choose(n, 2)
  gi <- rowSums(Phi) / (n - 1)
  Pt <- Phi - outer(gi, rep(1, n)) - outer(rep(1, n), gi) + Un; diag(Pt) <- 0
  ev <- eigen(Pt, symmetric = TRUE, only.values = TRUE)$values / n
  lam <- ev[order(abs(ev), decreasing = TRUE)][seq_len(min(r_eigen, length(ev)))]

  c2point <- function(idx) {
    ni <- length(idx)
    ti <- time[idx, , drop = FALSE]; si <- status[idx, , drop = FALSE]
    Xi <- X[idx]; Zi <- Z[idx]
    mb <- .mrwin_ar1_blocks(ti, si, Xi, Zi)
    max(ni * (mb$c0 - 2 * beta0 * mb$c1 + beta0^2 * mb$c2), 0)
  }
  boots <- vapply(seq_len(boot_reps), function(b) {
    tryCatch(c2point(sample.int(n, n, replace = TRUE)), error = function(e) NA_real_)
  }, numeric(1))
  boots <- boots[is.finite(boots)]
  c2_hi <- as.numeric(stats::quantile(boots, 1 - alpha1, na.rm = TRUE))
  alpha2 <- alpha - alpha1
  Z0 <- stats::rnorm(mc_reps); Zk <- matrix(stats::rnorm(mc_reps * length(lam)), mc_reps, length(lam))
  Rc <- 2 * sqrt(c2_hi) * Z0 + as.numeric((Zk^2 - 1) %*% lam)
  crit <- as.numeric(stats::quantile(Rc^2, 1 - alpha2, na.rm = TRUE))
  stat <- (n * Un)^2
  list(stat = stat, crit = crit, reject = stat > crit, c2_hi = c2_hi, n = n)
}

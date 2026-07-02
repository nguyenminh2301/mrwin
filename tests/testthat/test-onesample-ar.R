# One-sample Anderson-Rubin inference (Paper 03, eq. 4-9): mrwin_ar_onesample(),
# mrwin_ar_onesample_overid(), mrwin_ar_onesample_supci_test(). The core algebra
# is ported unchanged from the validated closures in
# tools/validation-scripts/p3-weak-iv-robust-probe.R; the tests below include a
# direct differential check against that reference (lesson H: a new numerical
# path is not trusted until its parity/calibration probe passes).

# reference closures, copied verbatim from the validated probe script
.ar1_ref_blocks <- function(time, status, X, Z) {
  subj <- getFromNamespace("mrwin_subject_win_loss_cpp", "mrwin")
  n <- length(X)
  cu <- subj(time, status, time, status, rep(1, n)); s_i <- cu[, 1] - cu[, 2]
  cz <- subj(time, status, time, status, Z); r_i <- cz[, 1] - cz[, 2]
  a <- 2 * sum(Z * s_i) / (n * (n - 1))
  b <- 2 * (n * sum(X * Z) - sum(X) * sum(Z)) / (n * (n - 1))
  gh <- (Z * s_i - r_i) / (n - 1); SX <- sum(X); SZ <- sum(Z); SXZ <- sum(X * Z)
  gx <- ((n - 1) * X * Z - X * (SZ - Z) - Z * (SX - X) + (SXZ - X * Z)) / (n - 1)
  list(a = a, b = b, c0 = stats::var(gh), c1 = stats::cov(gh, gx), c2 = stats::var(gx), n = n)
}

test_that("mrwin_ar_onesample matches the validated reference algebra exactly (weak instrument -> unbounded)", {
  cfg <- mrwin_config(n_outcome = 3000, m_snps = 1, seed = 42)
  sim <- mrwin_simulate(cfg, seed = 42)
  X <- sim$X; Z <- as.numeric(sim$G[, 1])
  ref <- .ar1_ref_blocks(sim$time, sim$status, X, Z)
  pkg <- mrwin_ar_onesample(sim$time, sim$status, X, Z, degeneracy_check = FALSE)
  expect_equal(pkg$gamma, ref$a / ref$b, tolerance = 1e-8)
  expect_false(pkg$bounded)                    # single noisy SNP: weak instrument
  expect_equal(pkg$ci, c(-Inf, Inf))
})

test_that("mrwin_ar_onesample: closed-form CI matches quadratic inversion of the reference AR stat (strong instrument -> bounded)", {
  cfg <- mrwin_config(n_outcome = 6000, m_snps = 20, seed = 11)
  sim <- mrwin_simulate(cfg, seed = 11)
  Zc <- as.numeric(sim$G %*% sim$true_betas)   # the actual combined score driving X: a strong instrument
  ref <- .ar1_ref_blocks(sim$time, sim$status, sim$X, Zc)
  q <- stats::qchisq(0.95, 1)
  ARstat_ref <- function(beta) ref$n * (ref$a - ref$b * beta)^2 / (4 * (ref$c0 - 2 * beta * ref$c1 + beta^2 * ref$c2))
  pkg <- mrwin_ar_onesample(sim$time, sim$status, sim$X, Zc, degeneracy_check = FALSE)
  expect_true(pkg$bounded)
  # the CI endpoints must each be (approximately) roots of AR(beta) = qchisq(.95,1)
  expect_equal(ARstat_ref(pkg$ci[1]), q, tolerance = 1e-4)
  expect_equal(ARstat_ref(pkg$ci[2]), q, tolerance = 1e-4)
  expect_true(pkg$gamma > pkg$ci[1] && pkg$gamma < pkg$ci[2])   # point estimate inside its own CI
})

test_that("mrwin_ar_onesample degeneracy bootstrap returns a sane c^2 CI and a logical flag", {
  cfg <- mrwin_config(n_outcome = 1500, m_snps = 1, seed = 3)
  sim <- mrwin_simulate(cfg, seed = 3)
  r <- mrwin_ar_onesample(sim$time, sim$status, sim$X, as.numeric(sim$G[, 1]),
                           boot_reps = 40L, seed = 9)
  expect_true(is.logical(r$degenerate))
  expect_false(is.na(r$degenerate))
  expect_length(r$c2_boot_ci, 2)
  expect_true(r$c2_boot_ci[1] <= r$c2_boot_ci[2])
  expect_true(r$c2_hat >= 0)
})

test_that("mrwin_ar_onesample_overid matches the validated reference over-ID statistic exactly", {
  cfg <- mrwin_config(n_outcome = 2500, m_snps = 15, seed = 21)
  sim <- mrwin_simulate(cfg, seed = 21)
  subj <- getFromNamespace("mrwin_subject_win_loss_cpp", "mrwin")
  n <- nrow(sim$G); L <- ncol(sim$G)
  s_i <- { cu <- subj(sim$time, sim$status, sim$time, sim$status, rep(1, n)); cu[, 1] - cu[, 2] }
  Uh <- 2 * as.numeric(crossprod(sim$G, s_i)) / (n * (n - 1))
  Ux <- 2 * (n * as.numeric(crossprod(sim$G, sim$X)) - colSums(sim$G) * sum(sim$X)) / (n * (n - 1))
  GH <- matrix(0, n, L)
  for (l in 1:L) { cz <- subj(sim$time, sim$status, sim$time, sim$status, sim$G[, l]); r <- cz[, 1] - cz[, 2]
    GH[, l] <- (sim$G[, l] * s_i - r) / (n - 1) }
  SX <- sum(sim$X); GX <- matrix(0, n, L)
  for (l in 1:L) { Zl <- sim$G[, l]; SZ <- sum(Zl); SXZ <- sum(sim$X * Zl)
    GX[, l] <- ((n - 1) * sim$X * Zl - sim$X * (SZ - Zl) - Zl * (SX - sim$X) + (SXZ - sim$X * Zl)) / (n - 1) }
  b0 <- sum(Ux * Uh) / sum(Ux * Ux); Gp <- GH - b0 * GX
  Si <- solve(4 * stats::cov(Gp) + diag(1e-8 * mean(diag(4 * stats::cov(Gp))), L))
  num <- as.numeric(t(Ux) %*% Si %*% Uh); den <- as.numeric(t(Ux) %*% Si %*% Ux)
  Qref <- n * (as.numeric(t(Uh) %*% Si %*% Uh) - num^2 / den)

  pkg <- mrwin_ar_onesample_overid(sim$time, sim$status, sim$X, sim$G)
  expect_equal(pkg$Q, Qref, tolerance = 1e-6)
  expect_equal(pkg$Q_df, L - 1L)
  expect_true(pkg$Q_p >= 0 && pkg$Q_p <= 1)
})

test_that("mrwin_ar_onesample_overid errors informatively on a single-column instrument", {
  cfg <- mrwin_config(n_outcome = 300, m_snps = 1, seed = 1)
  sim <- mrwin_simulate(cfg, seed = 1)
  expect_error(mrwin_ar_onesample_overid(sim$time, sim$status, sim$X, sim$G), "at least 2 columns")
})

test_that("mrwin_ar_onesample_supci_test runs end-to-end and returns a well-formed test result", {
  cfg <- mrwin_config(n_outcome = 400, m_snps = 1, seed = 7)
  sim <- mrwin_simulate(cfg, seed = 7)
  r <- mrwin_ar_onesample_supci_test(sim$time, sim$status, sim$X, as.numeric(sim$G[, 1]),
                                      beta0 = 0, boot_reps = 60L, mc_reps = 500L, seed = 5)
  expect_true(is.finite(r$stat) && r$stat >= 0)
  expect_true(is.finite(r$crit) && r$crit > 0)
  expect_true(is.logical(r$reject))
  expect_true(is.finite(r$c2_hi) && r$c2_hi >= 0)
})

test_that("mrwin_ar_onesample_supci_test does not reject its own point estimate (sanity, not a calibration claim)", {
  # A full calibration study (type-I error across the degeneracy continuum) is a
  # Monte Carlo claim, not a unit-test claim; it lives in
  # tools/validation-scripts/p3-supci-onesample-calibration.R and is reported in
  # the manuscript. This test only checks the packaged function is well-behaved:
  # testing the null AT its own basic-AR point estimate should essentially never
  # reject (the least-favorable critical value only gets larger than the naive
  # chi-sq one, never smaller, and the point estimate trivially satisfies U_n=0
  # up to the small-sample gap between the closed-form gamma and beta0).
  cfg <- mrwin_config(n_outcome = 600, m_snps = 1, seed = 17)
  sim <- mrwin_simulate(cfg, seed = 17)
  Z <- as.numeric(sim$G[, 1])
  basic <- mrwin_ar_onesample(sim$time, sim$status, sim$X, Z, degeneracy_check = FALSE)
  r <- mrwin_ar_onesample_supci_test(sim$time, sim$status, sim$X, Z, beta0 = basic$gamma,
                                      boot_reps = 60L, mc_reps = 500L, seed = 3)
  expect_false(isTRUE(r$reject))
})

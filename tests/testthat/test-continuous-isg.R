# WP15 / M1: continuous kernel-smoothed ISG. Step-1 internal check — the boxcar
# kernel with b = 1/D and centres k/D must reduce EXACTLY to the decile estimator
# (the same contrasts, covariance, and pooled Fieller inference).

test_that("continuous ISG reduces exactly to the decile estimator (boxcar, b=1/D)", {
  for (cfgi in list(c(N = 1000, D = 5, seed = 11),
                    c(N = 2000, D = 10, seed = 22))) {
    cfg <- mrwin_config(n_outcome = cfgi["N"], m_snps = 60, seed = cfgi["seed"])
    dat <- mrwin_simulate(cfg, seed = cfgi["seed"])
    D <- cfgi["D"]
    est <- mrwin_estimate(dat$time, dat$status, dat$G, dat$X, dat$true_betas, n_strata = D)
    ana <- mrwin_analytic_inference(est, dat$X)
    cont <- mrwin_continuous_isg(est$kernel, est$score, dat$X,
                                 bandwidth = 1 / D, kernel_fn = "boxcar",
                                 centers = (1:(D - 1)) / D)
    expect_equal(unname(cont$log_theta), unname(est$log_theta), tolerance = 1e-10)
    expect_equal(unname(cont$delta_x), unname(est$delta_x), tolerance = 1e-10)
    expect_equal(unname(cont$cov_u), unname(ana$cov_u), tolerance = 1e-10)
    expect_equal(cont$delta_gls, ana$delta_gls, tolerance = 1e-8)
    expect_equal(unname(cont$ci95_delta_fieller), unname(ana$ci95_delta_fieller), tolerance = 1e-8)
  }
})

test_that("continuous ISG runs with a smooth kernel and a centre grid", {
  cfg <- mrwin_config(n_outcome = 1500, m_snps = 60, seed = 7)
  dat <- mrwin_simulate(cfg, seed = 7)
  est <- mrwin_estimate(dat$time, dat$status, dat$G, dat$X, dat$true_betas, n_strata = 5)
  cont <- mrwin_continuous_isg(est$kernel, est$score, dat$X,
                               bandwidth = 0.12, n_centers = 20, kernel_fn = "epanechnikov")
  expect_true(is.finite(cont$delta_gls))
  expect_true(cont$n_centers >= 1L)
  expect_equal(dim(cont$cov_u), c(2L * cont$n_centers, 2L * cont$n_centers))
  expect_true(isSymmetric(unname(cont$cov_u), tol = 1e-8))
})

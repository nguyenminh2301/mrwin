# Direct coverage for the top-level analytic engine mrwin_analytic_bootstrap(),
# which mrwin(inference = "analytic") dispatches to.

test_that("mrwin_analytic_bootstrap is output-compatible and shares the bootstrap point", {
  cfg <- mrwin_config(n_outcome = 700, m_snps = 12, seed = 414)
  dat <- mrwin_simulate(cfg, seed = 414)

  ab <- mrwin_analytic_bootstrap(
    dat$time, dat$status, dat$G, dat$X, dat$true_betas,
    sigma_beta = 0, n_strata = 5, B_gwas = 4, seed = 414
  )
  mb <- mrwin_multiplier_bootstrap(
    dat$time, dat$status, dat$G, dat$X, dat$true_betas,
    sigma_beta = 0, n_strata = 5, B = 4, seed = 414, adjustment = "none"
  )

  expect_s3_class(ab, "mrwin_bootstrap")
  expect_true(all(c(
    "delta_gls", "se_delta_gls", "dscwr", "ci95_delta_fieller",
    "ci95_dscwr_fieller", "fieller_p_value", "covariance_method", "contrast_plan"
  ) %in% names(ab)))

  # The point estimate is deterministic in the data and does not depend on the
  # inference engine, so analytic and bootstrap must agree on it exactly.
  expect_equal(ab$point_log_theta, mb$point_log_theta)
  expect_equal(ab$point_delta_x, mb$point_delta_x)
  expect_equal(ab$delta_isg, mb$delta_isg)
  expect_identical(as.integer(ab$strata), as.integer(mb$strata))

  # sigma_beta = 0 means no GWAS-resample iterations.
  expect_identical(ab$B, 0L)
  expect_match(ab$covariance_method, "analytic")
})

test_that("mrwin_analytic_bootstrap adds the GWAS-resample term when sigma_beta > 0", {
  cfg <- mrwin_config(n_outcome = 600, m_snps = 12, seed = 415)
  dat <- mrwin_simulate(cfg, seed = 415)
  ab <- mrwin_analytic_bootstrap(
    dat$time, dat$status, dat$G, dat$X, dat$true_betas,
    sigma_beta = 0.15, n_strata = 5, B_gwas = 200, seed = 415
  )
  expect_identical(ab$B, 200L)
  expect_true(is.finite(ab$se_delta_gls))
  expect_match(ab$covariance_method, "gwas")
})

test_that("mrwin_analytic_bootstrap honours the stratification argument", {
  cfg <- mrwin_config(n_outcome = 600, m_snps = 10, seed = 416)
  dat <- mrwin_simulate(cfg, seed = 416)
  ab <- mrwin_analytic_bootstrap(
    dat$time, dat$status, dat$G, dat$X, dat$true_betas,
    sigma_beta = 0, n_strata = 5, B_gwas = 2, seed = 416,
    stratification = "doubly_ranked"
  )
  score <- as.numeric(dat$G %*% dat$true_betas)
  expect_identical(
    as.integer(ab$strata),
    mrwin_doubly_ranked_strata(score, dat$X, n_strata = 5)$strata
  )
})

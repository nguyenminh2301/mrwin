wp4_oracle <- function() {
  path <- system.file("extdata/oracle/wp4_estimator_small.R", package = "mrwin")
  if (identical(path, "")) {
    path <- file.path("..", "..", "inst", "extdata", "oracle", "wp4_estimator_small.R")
  }
  env <- new.env(parent = baseenv())
  sys.source(path, env)
  env$wp4_estimator_small
}

test_that("simulation and point estimator return adjacent contrasts", {
  cfg <- mrwin_config(n_outcome = 80, m_snps = 6, seed = 1)
  dat <- mrwin_simulate(cfg)

  fit <- mrwin_estimate(
    time = dat$time,
    status = dat$status,
    G = dat$G,
    X = dat$X,
    beta_hat = dat$true_betas,
    n_strata = 4
  )

  expect_length(fit$log_theta, 3)
  expect_length(fit$cwr, 3)
  expect_length(fit$delta_x, 3)
  expect_length(fit$delta_isg, 3)
  expect_named(fit$delta_isg, c("s2_vs_s1", "s3_vs_s2", "s4_vs_s3"))
  expect_true(all(is.finite(fit$log_theta)))
})

test_that("WP4 point estimator matches Python oracle fixture", {
  oracle <- wp4_oracle()
  fit <- mrwin_estimate(
    time = oracle$time,
    status = oracle$status,
    G = oracle$G,
    X = oracle$X,
    beta_hat = oracle$beta_hat,
    n_strata = oracle$n_strata
  )

  expect_equal(fit$score, oracle$expected$score, tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(fit$strata, oracle$expected$strata)
  expect_equal(fit$wins, oracle$expected$wins, tolerance = 1e-12)
  expect_equal(fit$losses, oracle$expected$losses, tolerance = 1e-12)
  expect_equal(fit$total, oracle$expected$total, tolerance = 1e-12)
  expect_equal(fit$log_theta, oracle$expected$log_theta, tolerance = 1e-12)
  expect_equal(fit$cwr, oracle$expected$cwr, tolerance = 1e-12)
  expect_equal(fit$delta_x, oracle$expected$delta_x, tolerance = 1e-12)
  expect_equal(fit$delta_isg, oracle$expected$delta_isg, tolerance = 1e-12)
})

test_that("WP4 GLS pooling and Q statistic match Python oracle fixture", {
  oracle <- wp4_oracle()
  pooled <- mrwin_gls_pool(oracle$expected$delta_isg, oracle$sigma_isg, shrink = FALSE)

  expect_equal(pooled$delta_gls, oracle$gls_no_shrink$delta_gls, tolerance = 1e-12)
  expect_equal(pooled$se_delta_gls, oracle$gls_no_shrink$se_delta_gls, tolerance = 1e-12)
  expect_equal(pooled$dscwr, oracle$gls_no_shrink$dscwr, tolerance = 1e-12)
  expect_equal(pooled$q, oracle$gls_no_shrink$q, tolerance = 1e-12)
  expect_equal(pooled$q_df, oracle$gls_no_shrink$q_df)
  expect_equal(sum(pooled$gls_weights), 1, tolerance = 1e-12)
})

test_that("GLS pooling matches a hand-calculated diagonal example", {
  delta <- c(a = 2, b = 4)
  sigma <- diag(c(1, 4))
  pooled <- mrwin_gls_pool(delta, sigma, shrink = FALSE)

  expect_equal(pooled$delta_gls, 2.4, tolerance = 1e-12)
  expect_equal(pooled$se_delta_gls, sqrt(0.8), tolerance = 1e-12)
  expect_equal(pooled$q, 0.8, tolerance = 1e-12)
  expect_equal(pooled$q_df, 1)
  expect_equal(pooled$q_p_value, stats::pchisq(0.8, df = 1, lower.tail = FALSE))
  expect_equal(pooled$gls_weights, c(a = 0.8, b = 0.2), tolerance = 1e-12)
})

test_that("single-contrast GLS reports no heterogeneity test", {
  pooled <- mrwin_gls_pool(c(a = 2), matrix(4, nrow = 1), shrink = FALSE)

  expect_equal(pooled$delta_gls, 2)
  expect_equal(pooled$q, 0)
  expect_equal(pooled$q_df, 0)
  expect_true(is.na(pooled$q_p_value))
})

test_that("WP4 estimator and GLS validate invalid inputs", {
  oracle <- wp4_oracle()

  expect_error(
    mrwin_estimate(
      oracle$time,
      oracle$status,
      oracle$G,
      oracle$X,
      oracle$beta_hat,
      n_strata = nrow(oracle$G) + 1L
    ),
    "cannot exceed"
  )
  expect_error(
    mrwin_estimate(
      oracle$time,
      oracle$status,
      oracle$G,
      oracle$X,
      oracle$beta_hat,
      n_strata = 2.5
    ),
    "single integer"
  )
  expect_error(
    mrwin_estimate(
      oracle$time,
      oracle$status,
      oracle$G,
      oracle$X,
      oracle$beta_hat,
      n_strata = oracle$n_strata,
      weights = c(rep(1, nrow(oracle$G) - 1L), -1)
    ),
    "finite and non-negative"
  )
  expect_error(
    mrwin_gls_pool(c(1, NA), diag(2), shrink = FALSE),
    "finite values"
  )
  expect_error(
    mrwin_gls_pool(c(1, 2), matrix(c(1, 0.5, 0.1, 1), nrow = 2), shrink = FALSE),
    "symmetric"
  )
  expect_error(
    mrwin_gls_pool(c(1, 2), matrix(c(1, 2, 2, 1), nrow = 2), shrink = FALSE),
    "positive semi-definite"
  )
})

test_that("multiplier bootstrap returns GLS DS-CWR", {
  cfg <- mrwin_config(n_outcome = 80, m_snps = 6, seed = 2)
  dat <- mrwin_simulate(cfg)

  boot <- mrwin_multiplier_bootstrap(
    time = dat$time,
    status = dat$status,
    G = dat$G,
    X = dat$X,
    beta_hat = dat$true_betas,
    sigma_beta = rep(0.01, length(dat$true_betas)),
    n_strata = 4,
    B = 20,
    seed = 3
  )

  expect_true(is.finite(boot$delta_gls))
  expect_true(is.finite(boot$se_delta_gls))
  expect_true(is.finite(boot$dscwr))
  expect_equal(nrow(boot$sigma_isg), 3)
})

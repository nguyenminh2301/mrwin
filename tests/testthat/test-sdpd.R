test_that("MR-Egger returns intercept diagnostics", {
  beta_x <- c(0.1, 0.2, 0.3, 0.4)
  beta_y <- 0.02 + 0.5 * beta_x
  se_y <- rep(0.1, 4)

  fit <- mrwin_mr_egger(beta_x, beta_y, se_y)

  expect_equal(unname(fit$intercept), 0.02, tolerance = 1e-8)
  expect_true(is.finite(fit$intercept_p_value))
})

test_that("pleiotropy bounded CI widens sampling CI", {
  ci <- mrwin_pleiotropy_bounded_ci(delta_hat = -0.1, se_delta = 0.02, bias_radius = 0.05)

  expect_lt(ci$ci_delta_pleiotropy_bounded[1], ci$ci_delta_sampling[1])
  expect_gt(ci$ci_delta_pleiotropy_bounded[2], ci$ci_delta_sampling[2])
})

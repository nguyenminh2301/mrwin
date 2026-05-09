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
  expect_length(fit$delta_x, 3)
  expect_length(fit$delta_isg, 3)
  expect_true(all(is.finite(fit$log_theta)))
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

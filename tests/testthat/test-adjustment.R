wp6_toy <- function() {
  n <- 16L
  score <- seq(-2, 2, length.out = n)
  B <- 6L
  list(
    time = matrix(seq_len(n), ncol = 1L),
    status = matrix(rep(1, n), ncol = 1L),
    G = matrix(score, ncol = 1L),
    X = seq(0, 3, length.out = n),
    beta = 1,
    Z = matrix(c(
      -2.0, -1.8, -1.6, -1.4,
      -1.0, -0.8, 0.8, 1.0,
      0.2, 0.4, 0.6, 0.8,
      1.2, 1.4, 1.6, 1.8
    ), ncol = 1L, dimnames = list(NULL, "z")),
    beta_draws = matrix(rep(1, B), nrow = B, byrow = TRUE),
    multiplier_weights = matrix(c(
      1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
      1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
      1.2, 0.8, 1.1, 0.9, 1.3, 0.7, 1.4, 0.6,
      0.9, 1.1, 1.5, 0.7, 1.2, 0.8, 1.3, 0.6,
      0.7, 1.4, 0.8, 1.3, 0.9, 1.2, 1.0, 1.1,
      1.5, 1.0, 0.6, 1.4, 0.8, 1.3, 0.7, 1.2,
      0.9, 1.1, 1.5, 0.7, 1.2, 0.8, 1.3, 0.6,
      1.3, 0.6, 1.2, 0.8, 1.1, 0.9, 1.5, 0.7,
      1.5, 1.0, 0.6, 1.4, 0.8, 1.3, 0.7, 1.2,
      0.9, 1.1, 1.5, 0.7, 1.2, 0.8, 1.3, 0.6,
      1.3, 0.6, 1.2, 0.8, 1.1, 0.9, 1.5, 0.7,
      1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
    ), nrow = B, byrow = TRUE)
  )
}

test_that("WP6 ordinal IPTW reports ESS, balance, and dropped strata", {
  toy <- wp6_toy()
  strata <- rep(seq_len(4L), each = 4L)

  adj <- mrwin_propensity_weights(
    strata = strata,
    covariates = toy$Z,
    ess_fraction = 0.8
  )

  expect_s3_class(adj, "mrwin_propensity_weights")
  expect_equal(adj$dropped_strata, 2L)
  expect_equal(adj$active_strata, c(1L, 3L, 4L))
  expect_true(all(is.finite(adj$weights)))
  expect_true(adj$ess["s2"] < adj$ess_threshold["s2"])
  expect_true(all(adj$ess[c("s1", "s3", "s4")] >= adj$ess_threshold[c("s1", "s3", "s4")]))
  expect_equal(adj$balance$covariate, "z")
  expect_true(is.finite(adj$balance$before_max_abs_smd))
  expect_true(is.finite(adj$balance$after_max_abs_smd))
})

test_that("WP6 point estimator bridges across positivity-filtered strata", {
  toy <- wp6_toy()

  fit <- mrwin_estimate(
    time = toy$time,
    status = toy$status,
    G = toy$G,
    X = toy$X,
    beta_hat = toy$beta,
    n_strata = 4L,
    active_strata = c(1L, 3L, 4L)
  )

  expect_equal(fit$active_strata, c(1L, 3L, 4L))
  expect_equal(fit$contrast_plan$label, c("s3_vs_s1", "s4_vs_s3"))
  expect_equal(fit$contrast_plan$bridged, c(TRUE, FALSE))
  expect_named(fit$delta_isg, c("s3_vs_s1", "s4_vs_s3"))
})

test_that("WP6 bootstrap refits IPTW and preserves bridged contrast schema", {
  toy <- wp6_toy()

  boot <- mrwin_multiplier_bootstrap(
    time = toy$time,
    status = toy$status,
    G = toy$G,
    X = toy$X,
    beta_hat = toy$beta,
    sigma_beta = 0,
    n_strata = 4L,
    B = nrow(toy$beta_draws),
    beta_draws = toy$beta_draws,
    multiplier_weights = toy$multiplier_weights,
    covariates = toy$Z,
    adjustment = "ordinal_iptw",
    ess_fraction = 0.8
  )

  expect_s3_class(boot, "mrwin_bootstrap")
  expect_equal(boot$active_strata, c(1L, 3L, 4L))
  expect_equal(boot$dropped_strata, 2L)
  expect_equal(boot$contrast_plan$label, c("s3_vs_s1", "s4_vs_s3"))
  expect_true(boot$adjustment$has_positivity_failure)
  expect_true(boot$adjustment$has_bridging)
  expect_equal(ncol(boot$bootstrap_log_theta), 2L)
  expect_equal(ncol(boot$bootstrap_ess), 4L)
  expect_true(all(boot$bootstrap_dropped_strata[, "s2"]))
  expect_equal(boot$n_valid, nrow(toy$beta_draws))
  expect_equal(boot$n_invalid, 0L)
})

test_that("WP6 adjustment validates covariates, truncation, and future GPS mode", {
  toy <- wp6_toy()
  strata <- rep(seq_len(4L), each = 4L)

  expect_error(
    mrwin_propensity_weights(strata, matrix(1, nrow = length(strata), ncol = 1L)),
    "positive variance"
  )
  expect_error(
    mrwin_propensity_weights(strata, toy$Z, truncate = c(0.9, 0.1)),
    "`iptw_truncation`"
  )
  expect_error(
    mrwin_multiplier_bootstrap(
      toy$time,
      toy$status,
      toy$G,
      toy$X,
      toy$beta,
      n_strata = 4L,
      B = 6L,
      adjustment = "ordinal_iptw"
    ),
    "`covariates` are required"
  )
  expect_error(
    mrwin(
      endpoint = mrwin_endpoint(toy$time, toy$status),
      genotype = toy$G,
      exposure = toy$X,
      gwas = mrwin_gwas(toy$beta, 0.01),
      covariates = toy$Z,
      controls = mrwin_controls(n_strata = 4L, bootstrap = 6L, adjustment = "gps")
    ),
    "planned for a later work package"
  )
})

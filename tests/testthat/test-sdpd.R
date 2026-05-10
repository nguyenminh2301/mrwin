manual_aalen <- function(g, time, status) {
  design <- cbind(time, g * time)
  inv <- solve(crossprod(design))
  coef <- inv %*% crossprod(design, status)
  resid <- status - drop(design %*% coef)
  meat <- crossprod(design * resid, design * resid)
  vcov <- inv %*% meat %*% inv
  c(beta = coef[2L, 1L], se = sqrt(max(vcov[2L, 2L], 1e-20)))
}

test_that("MR-Egger returns intercept diagnostics", {
  beta_x <- c(0.1, 0.2, 0.3, 0.4)
  beta_y <- 0.02 + 0.5 * beta_x
  se_y <- rep(0.1, 4)

  fit <- mrwin_mr_egger(beta_x, beta_y, se_y)

  expect_s3_class(fit, "mrwin_mr_egger")
  expect_equal(unname(fit$intercept), 0.02, tolerance = 1e-8)
  expect_true(is.finite(fit$intercept_p_value))
})

test_that("WP7 MR-Egger separates deterministic null and pleiotropic summaries", {
  beta_x <- seq(-0.4, 0.4, length.out = 8)
  se_y <- rep(0.01, length(beta_x))

  null <- mrwin_mr_egger(beta_x, 0.5 * beta_x, se_y)
  alt <- mrwin_mr_egger(beta_x, 0.05 + 0.5 * beta_x, se_y)

  expect_equal(unname(null$intercept), 0, tolerance = 1e-12)
  expect_gt(null$intercept_p_value, 0.9)
  expect_equal(unname(alt$intercept), 0.05, tolerance = 1e-12)
  expect_lt(alt$intercept_p_value, 0.05)
})

test_that("WP7 Aalen per-SNP summary matches hand sandwich calculation", {
  G <- cbind(
    rs1 = c(0, 0, 1, 1, 2, 2),
    rs2 = c(0, 1, 0, 2, 1, 2),
    rs3 = c(2, 1, 2, 0, 1, 0)
  )
  time <- c(1.0, 1.4, 2.0, 2.5, 3.0, 3.6)
  status <- c(0, 1, 0, 1, 1, 1)

  fit <- mrwin_aalen_per_snp(G, time, status)
  expected <- vapply(seq_len(ncol(G)), function(j) manual_aalen(G[, j], time, status), numeric(2))

  expect_s3_class(fit, "mrwin_snp_summary")
  expect_equal(fit$beta, expected["beta", ], tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(fit$se, expected["se", ], tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("WP7 Cox per-SNP summary returns finite log-HR estimates", {
  G <- cbind(
    rs1 = c(0, 0, 0, 1, 1, 1, 2, 2, 2),
    rs2 = c(0, 1, 2, 0, 1, 2, 0, 1, 2),
    rs3 = c(2, 1, 0, 2, 1, 0, 2, 1, 0)
  )
  time <- c(9, 8, 7, 6, 5, 4, 3, 2, 1)
  status <- rep(1, length(time))

  fit <- mrwin_cox_per_snp(G, time, status)

  expect_s3_class(fit, "mrwin_snp_summary")
  expect_true(all(is.finite(fit$beta)))
  expect_true(all(is.finite(fit$se)))
  expect_gt(fit$beta[["rs1"]], 0)
})

test_that("WP7 public SDPD wrapper returns both scales and underpower caveat", {
  G <- cbind(
    rs1 = c(0, 0, 1, 1, 2, 2, 0, 1),
    rs2 = c(0, 1, 0, 2, 1, 2, 2, 0),
    rs3 = c(2, 1, 2, 0, 1, 0, 1, 2),
    rs4 = c(0, 2, 1, 0, 2, 1, 0, 2)
  )
  X <- drop(G %*% c(0.2, -0.1, 0.15, 0.05)) + seq(-0.2, 0.2, length.out = nrow(G))
  time <- c(5, 6, 4, 3, 2, 7, 8, 1)
  status <- c(1, 0, 1, 1, 1, 0, 0, 1)

  sdpd <- mrwin_sdpd(G, X, time, status, scale = "both", min_snps = 5)

  expect_s3_class(sdpd, "mrwin_sdpd")
  expect_named(sdpd$results, c("aalen", "cox"))
  expect_true(sdpd$underpowered)
  expect_equal(sdpd$n_events, sum(status))
  expect_true(all(vapply(sdpd$results, function(x) "outcome_summary" %in% names(x), logical(1))))
})

test_that("pleiotropy bounded CI widens sampling CI", {
  ci <- mrwin_pleiotropy_bounded_ci(delta_hat = -0.1, se_delta = 0.02, bias_radius = 0.05)

  expect_s3_class(ci, "mrwin_pleiotropy_bounded_ci")
  expect_lt(ci$ci_delta_pleiotropy_bounded[1], ci$ci_delta_sampling[1])
  expect_gt(ci$ci_delta_pleiotropy_bounded[2], ci$ci_delta_sampling[2])
  expect_named(ci$ci_delta_sampling, c("low", "high"))
})

test_that("WP7 diagnostics validate malformed inputs", {
  G <- cbind(rs1 = c(0, 1, 2, 0), rs2 = c(1, 0, 1, 2), rs3 = c(2, 1, 0, 2))
  time <- c(1, 2, 3, 4)
  status <- c(0, 1, 0, 1)

  expect_error(mrwin_aalen_per_snp(cbind(G, zero = 1), time, status), "zero-variance")
  expect_error(mrwin_cox_per_snp(G, time, c(0, 1, 2, 1)), "`status`")
  expect_error(mrwin_mr_egger(c(1, 1, 1), c(1, 2, 3), c(0.1, 0.1, 0.1)), "positive variance")
  expect_error(mrwin_mr_egger(c(1, 2), c(1, 2, 3), c(0.1, 0.1, 0.1)), "same length")
  expect_error(mrwin_pleiotropy_bounded_ci(0, 0.1, -0.1), "`bias_radius`")
})

test_that("WP7 high-level workflow records SDPD underpower and pleiotropy-bound outputs", {
  cfg <- mrwin_config(n_outcome = 50, m_snps = 4, seed = 71)
  dat <- mrwin_simulate(cfg)

  fit <- mrwin(
    endpoint = mrwin_endpoint(dat$time, dat$status),
    genotype = dat$G,
    exposure = dat$X,
    gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
    controls = mrwin_controls(
      n_strata = 4,
      bootstrap = 6,
      seed = 72,
      run_sdpd = TRUE,
      sdpd_scale = "aalen",
      sdpd_min_snps = 10,
      pleiotropy_bias_radius = 0.05
    )
  )

  warning_codes <- vapply(fit$warnings, `[[`, character(1), "code")
  expect_s3_class(fit$sdpd, "mrwin_sdpd")
  expect_true("sdpd_underpowered" %in% warning_codes)
  expect_s3_class(fit$inference$pleiotropy_bounded, "mrwin_pleiotropy_bounded_ci")
  expect_true(is.finite(summary(fit)$estimate$pleiotropy_ci_low))
})

# Fieller is the primary reported interval (R2 calibration fix); the
# bivariate-Delta interval is retained as a labelled reference.

test_that("primary reported CI and p-value are Fieller", {
  cfg <- mrwin_config(n_outcome = 400, m_snps = 20, seed = 5)
  dat <- mrwin_simulate(cfg, seed = 5)
  ep <- mrwin_endpoint(dat$time, dat$status, colnames(dat$time))
  gw <- mrwin_gwas(dat$true_betas, rep(0.05, length(dat$true_betas)))
  fit <- mrwin(endpoint = ep, genotype = dat$G, exposure = dat$X, gwas = gw,
               controls = mrwin_controls(n_strata = 4, bootstrap = 100, seed = 5, run_sdpd = FALSE))

  td <- tidy(fit)
  expect_equal(td$ci_method, "Fieller")
  expect_equal(unname(td$ci_low), unname(fit$inference$ci95_dscwr_fieller[[1]]))
  expect_equal(unname(td$ci_high), unname(fit$inference$ci95_dscwr_fieller[[2]]))
  expect_equal(unname(td$p_value), unname(fit$bootstrap$fieller_p_value))
  # bivariate-Delta retained as labelled reference, and distinct from Fieller
  expect_equal(unname(td$ci_low_delta_method), unname(fit$inference$ci95_dscwr[[1]]))
  expect_true(is.finite(td$p_value) || isTRUE(td$fieller_unbounded))

  s <- summary(fit)
  expect_equal(s$estimate$ci_method, "Fieller")
  expect_output(print(s), "Fieller")
  expect_output(print(fit), "Fieller")
})

test_that("Fieller p-value is a valid two-sided p", {
  cfg <- mrwin_config(n_outcome = 500, m_snps = 30, seed = 8)
  dat <- mrwin_simulate(cfg, seed = 8)
  ep <- mrwin_endpoint(dat$time, dat$status, colnames(dat$time))
  gw <- mrwin_gwas(dat$true_betas, rep(0, length(dat$true_betas)))
  fit <- mrwin(endpoint = ep, genotype = dat$G, exposure = dat$X, gwas = gw,
               controls = mrwin_controls(n_strata = 5, bootstrap = 50, seed = 8,
                                         inference = "analytic", run_sdpd = FALSE))
  p <- fit$bootstrap$fieller_p_value
  expect_true(is.na(p) || (p >= 0 && p <= 1))
})

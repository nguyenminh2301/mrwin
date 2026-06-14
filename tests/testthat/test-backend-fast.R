# WP13 / S1 integration: backend = "fast" routes through mrwin() and matches the
# sparse/dense backends (fast path for K=1, dense fallback for K>1).

.fit_with_backend <- function(endpoint, dat, backend, seed) {
  mrwin(
    endpoint = endpoint,
    genotype = dat$G,
    exposure = dat$X,
    gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
    controls = mrwin_controls(
      n_strata = 4, bootstrap = 20, seed = seed,
      backend = backend, run_sdpd = FALSE
    )
  )
}

test_that("backend='fast' equals 'sparse' on a single-endpoint (K=1) cohort", {
  cfg <- mrwin_config(n_outcome = 300, m_snps = 12, seed = 41)
  dat <- mrwin_simulate(cfg, seed = 41)
  ep <- mrwin_endpoint(
    dat$time[, 1, drop = FALSE], dat$status[, 1, drop = FALSE], "death"
  )
  fit_sparse <- .fit_with_backend(ep, dat, "sparse", 7)
  fit_fast <- .fit_with_backend(ep, dat, "fast", 7)

  expect_equal(fit_fast$point$delta_gls, fit_sparse$point$delta_gls, tolerance = 1e-8)
  expect_equal(fit_fast$point$dscwr, fit_sparse$point$dscwr, tolerance = 1e-8)
  expect_equal(fit_fast$point$log_theta, fit_sparse$point$log_theta, tolerance = 1e-8)
  expect_equal(fit_fast$inference$se_delta_gls, fit_sparse$inference$se_delta_gls,
               tolerance = 1e-8)
  expect_equal(fit_fast$heterogeneity$q, fit_sparse$heterogeneity$q, tolerance = 1e-8)
})

test_that("backend='fast' falls back and equals 'sparse' on K=3 endpoints", {
  cfg <- mrwin_config(n_outcome = 300, m_snps = 12, seed = 42)
  dat <- mrwin_simulate(cfg, seed = 42)
  ep <- mrwin_endpoint(dat$time, dat$status, colnames(dat$time))
  fit_sparse <- .fit_with_backend(ep, dat, "sparse", 9)
  fit_fast <- .fit_with_backend(ep, dat, "fast", 9)

  expect_equal(fit_fast$point$delta_gls, fit_sparse$point$delta_gls, tolerance = 1e-10)
  expect_equal(fit_fast$inference$se_delta_gls, fit_sparse$inference$se_delta_gls,
               tolerance = 1e-10)
})

test_that("mrwin_controls accepts backend='fast'", {
  ctrl <- mrwin_controls(backend = "fast")
  expect_equal(ctrl$backend, "fast")
})

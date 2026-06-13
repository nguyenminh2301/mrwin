# Rust backend parity tests.
#
# These are skipped unless the optional compiled extendr library is present
# (i.e. the package was installed with a Rust toolchain). When available, they
# assert that backend = "rust" reproduces backend = "sparse" to 1e-10, which in
# turn matches the dense backend (see test-backend-sparse.R).

test_that("rust backend matches sparse backend when available", {
  skip_if_not(
    isTRUE(mrwin:::.mrwin_rust_available()),
    "compiled Rust backend not installed"
  )

  cfg <- mrwin_config(n_outcome = 300, m_snps = 12, seed = 101)
  dat <- mrwin_simulate(cfg, seed = 101)

  mk <- function(backend) {
    mrwin(
      endpoint = mrwin_endpoint(dat$time, dat$status, colnames(dat$time)),
      genotype = dat$G,
      exposure = dat$X,
      gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
      controls = mrwin_controls(
        n_strata = 5, bootstrap = 50, seed = 7,
        backend = backend, run_sdpd = FALSE
      )
    )
  }

  fit_sparse <- mk("sparse")
  fit_rust <- mk("rust")

  expect_equal(fit_rust$point$delta_gls, fit_sparse$point$delta_gls, tolerance = 1e-10)
  expect_equal(fit_rust$point$dscwr, fit_sparse$point$dscwr, tolerance = 1e-10)
  expect_equal(fit_rust$inference$se_delta_gls, fit_sparse$inference$se_delta_gls, tolerance = 1e-10)
  expect_equal(
    unname(fit_rust$point$log_theta),
    unname(fit_sparse$point$log_theta),
    tolerance = 1e-10
  )
})

test_that("rust backend errors clearly when the library is absent", {
  skip_if(
    isTRUE(mrwin:::.mrwin_rust_available()),
    "compiled Rust backend is installed; absence path not exercised"
  )
  expect_error(
    mrwin_rust_bootstrap(
      time = matrix(1, 4, 2), status = matrix(0L, 4, 2),
      G = matrix(0, 4, 2), X = rep(0, 4), beta_hat = c(0, 0)
    ),
    "requires the compiled extendr library"
  )
})

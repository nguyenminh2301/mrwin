test_that("WP1 constructors validate core API objects", {
  endpoint <- mrwin_endpoint(
    time = c("t_death", "t_hf"),
    status = c("d_death", "d_hf"),
    priority = c("death", "hf")
  )
  expect_s3_class(endpoint, "mrwin_endpoint_spec")
  expect_equal(endpoint$priority, c("death", "hf"))

  gwas <- mrwin_gwas(beta = c(0.1, 0.2), se = c(0.01, 0.02), snp = c("rs1", "rs2"))
  expect_s3_class(gwas, "mrwin_gwas_spec")
  expect_equal(gwas$snp, c("rs1", "rs2"))

  controls <- mrwin_controls(n_strata = 4, bootstrap = 10, seed = 1)
  expect_s3_class(controls, "mrwin_controls")
  expect_equal(controls$n_strata, 4)
})

test_that("mrwin high-level workflow returns an mrwin_fit", {
  cfg <- mrwin_config(n_outcome = 80, m_snps = 6, seed = 11)
  dat <- mrwin_simulate(cfg)
  df <- data.frame(
    X = dat$X,
    dat$time,
    dat$status,
    dat$G
  )
  names(df) <- c(
    "X",
    "t_death", "t_hf", "t_renal",
    "d_death", "d_hf", "d_renal",
    paste0("snp", seq_len(ncol(dat$G)))
  )

  fit <- mrwin(
    data = df,
    endpoint = mrwin_endpoint(
      time = c("t_death", "t_hf", "t_renal"),
      status = c("d_death", "d_hf", "d_renal"),
      priority = c("death", "hf", "renal")
    ),
    genotype = paste0("snp", seq_len(ncol(dat$G))),
    exposure = "X",
    gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
    controls = mrwin_controls(n_strata = 4, bootstrap = 10, seed = 12, run_sdpd = FALSE)
  )

  expect_s3_class(fit, "mrwin_fit")
  expect_true(is.finite(fit$point$dscwr))
  expect_true(is.data.frame(summary(fit)$estimate))
})

test_that("mrwin blocks unsupported WP1 options explicitly", {
  cfg <- mrwin_config(n_outcome = 40, m_snps = 4, seed = 13)
  dat <- mrwin_simulate(cfg)

  expect_error(
    mrwin(
      endpoint = mrwin_endpoint(dat$time, dat$status),
      genotype = dat$G,
      exposure = dat$X,
      gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
      controls = mrwin_controls(n_strata = 4, bootstrap = 5, backend = "rcpp", run_sdpd = FALSE)
    ),
    "rcpp"
  )
})

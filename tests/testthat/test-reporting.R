test_that("print.mrwin_fit shows backend and Q statistic", {
  cfg <- mrwin_config(n_outcome = 80, m_snps = 6, seed = 201)
  dat <- mrwin_simulate(cfg)
  fit <- mrwin(
    endpoint = mrwin_endpoint(dat$time, dat$status, colnames(dat$time)),
    genotype = dat$G, exposure = dat$X,
    gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
    controls = mrwin_controls(n_strata = 4, bootstrap = 20, seed = 202, run_sdpd = FALSE)
  )
  out <- capture.output(print(fit))
  expect_true(any(grepl("Backend:", out)))
  expect_true(any(grepl("Q:", out)))
  expect_true(any(grepl("Adjustment:", out)))
})

test_that("summary.mrwin_fit includes p_value, fieller, diagnostics", {
  cfg <- mrwin_config(n_outcome = 80, m_snps = 6, seed = 203)
  dat <- mrwin_simulate(cfg)
  fit <- mrwin(
    endpoint = mrwin_endpoint(dat$time, dat$status, colnames(dat$time)),
    genotype = dat$G, exposure = dat$X,
    gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
    controls = mrwin_controls(n_strata = 4, bootstrap = 20, seed = 204, run_sdpd = FALSE)
  )
  s <- summary(fit)

  expect_true("p_value" %in% names(s$estimate))
  expect_true(is.finite(s$estimate$p_value))
  expect_s3_class(s$fieller, "data.frame")
  expect_true("unbounded" %in% names(s$fieller))
  expect_s3_class(s$diagnostics, "data.frame")
  expect_true("n_valid_bootstrap" %in% names(s$diagnostics))
  expect_true("cwr" %in% names(s$adjacent))
})

test_that("print.summary.mrwin_fit produces formatted output", {
  cfg <- mrwin_config(n_outcome = 80, m_snps = 6, seed = 205)
  dat <- mrwin_simulate(cfg)
  fit <- mrwin(
    endpoint = mrwin_endpoint(dat$time, dat$status, colnames(dat$time)),
    genotype = dat$G, exposure = dat$X,
    gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
    controls = mrwin_controls(n_strata = 4, bootstrap = 20, seed = 206, run_sdpd = FALSE)
  )
  out <- capture.output(print(summary(fit)))
  expect_true(any(grepl("Main Estimate", out)))
  expect_true(any(grepl("DS-CWR:", out)))
  expect_true(any(grepl("Heterogeneity", out)))
  expect_true(any(grepl("Fieller", out)))
})

test_that("plot.mrwin_fit supports isg, forest, and bootstrap types", {
  cfg <- mrwin_config(n_outcome = 80, m_snps = 6, seed = 207)
  dat <- mrwin_simulate(cfg)
  fit <- mrwin(
    endpoint = mrwin_endpoint(dat$time, dat$status, colnames(dat$time)),
    genotype = dat$G, exposure = dat$X,
    gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
    controls = mrwin_controls(n_strata = 4, bootstrap = 20, seed = 208, run_sdpd = FALSE)
  )

  expect_no_error(plot(fit, type = "isg"))
  expect_no_error(plot(fit, type = "forest"))
  expect_no_error(plot(fit, type = "bootstrap"))
  expect_error(plot(fit, type = "invalid"), "should be one of")
})

test_that("tidy.mrwin_fit returns a one-row data frame", {
  cfg <- mrwin_config(n_outcome = 80, m_snps = 6, seed = 209)
  dat <- mrwin_simulate(cfg)
  fit <- mrwin(
    endpoint = mrwin_endpoint(dat$time, dat$status, colnames(dat$time)),
    genotype = dat$G, exposure = dat$X,
    gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
    controls = mrwin_controls(n_strata = 4, bootstrap = 20, seed = 210, run_sdpd = FALSE)
  )
  td <- tidy(fit)

  expect_s3_class(td, "data.frame")
  expect_equal(nrow(td), 1)
  expect_true("estimate" %in% names(td))
  expect_true("p_value" %in% names(td))
  expect_true("ci_low" %in% names(td))
  expect_true("ci_high" %in% names(td))
  expect_true("n" %in% names(td))
  expect_true("m_snps" %in% names(td))
  expect_true("q_stat" %in% names(td))
  expect_equal(td$term, "DS-CWR")
  expect_true(is.finite(td$estimate))
  expect_true(is.finite(td$p_value))
})

test_that("mrwin_report produces text output", {
  cfg <- mrwin_config(n_outcome = 80, m_snps = 6, seed = 211)
  dat <- mrwin_simulate(cfg)
  fit <- mrwin(
    endpoint = mrwin_endpoint(dat$time, dat$status, colnames(dat$time)),
    genotype = dat$G, exposure = dat$X,
    gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
    controls = mrwin_controls(n_strata = 4, bootstrap = 20, seed = 212, run_sdpd = FALSE)
  )

  out <- capture.output(mrwin_report(fit, format = "text"))
  expect_true(any(grepl("mrwin Analysis Report", out)))
  expect_true(any(grepl("Data Summary", out)))
  expect_true(any(grepl("Main Estimate", out)))
  expect_true(any(grepl("DS-CWR", out)))
})

test_that("mrwin_report produces markdown output", {
  cfg <- mrwin_config(n_outcome = 80, m_snps = 6, seed = 213)
  dat <- mrwin_simulate(cfg)
  fit <- mrwin(
    endpoint = mrwin_endpoint(dat$time, dat$status, colnames(dat$time)),
    genotype = dat$G, exposure = dat$X,
    gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
    controls = mrwin_controls(n_strata = 4, bootstrap = 20, seed = 214, run_sdpd = FALSE)
  )

  out <- capture.output(mrwin_report(fit, format = "markdown"))
  expect_true(any(grepl("# mrwin Analysis Report", out)))
  expect_true(any(grepl("## Data Summary", out)))
  expect_true(any(grepl("## Main Estimate", out)))
  expect_true(any(grepl("| DS-CWR |", out)))
})

test_that("mrwin_report writes to file", {
  cfg <- mrwin_config(n_outcome = 80, m_snps = 6, seed = 215)
  dat <- mrwin_simulate(cfg)
  fit <- mrwin(
    endpoint = mrwin_endpoint(dat$time, dat$status, colnames(dat$time)),
    genotype = dat$G, exposure = dat$X,
    gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
    controls = mrwin_controls(n_strata = 4, bootstrap = 20, seed = 216, run_sdpd = FALSE)
  )

  tmp <- tempfile(fileext = ".md")
  on.exit(unlink(tmp))
  mrwin_report(fit, file = tmp, format = "markdown")
  expect_true(file.exists(tmp))
  content <- readLines(tmp)
  expect_true(any(grepl("mrwin Analysis Report", content)))
})

test_that("mrwin_report validates input", {
  expect_error(mrwin_report(list()), "mrwin_fit")
})

test_that("warnings show codes in print output", {
  cfg <- mrwin_config(n_outcome = 80, m_snps = 6, seed = 217)
  dat <- mrwin_simulate(cfg)
  fit <- mrwin(
    endpoint = mrwin_endpoint(dat$time, dat$status, colnames(dat$time)),
    genotype = dat$G, exposure = dat$X,
    gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
    controls = mrwin_controls(n_strata = 4, bootstrap = 20, seed = 218, run_sdpd = TRUE, sdpd_scale = "both")
  )
  out <- capture.output(print(fit))
  expect_true(any(grepl("Warnings", out)) || length(fit$warnings) == 0)
})

test_that("tidy works with sparse backend", {
  cfg <- mrwin_config(n_outcome = 80, m_snps = 6, seed = 219)
  dat <- mrwin_simulate(cfg)
  fit <- mrwin(
    endpoint = mrwin_endpoint(dat$time, dat$status, colnames(dat$time)),
    genotype = dat$G, exposure = dat$X,
    gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
    controls = mrwin_controls(n_strata = 4, bootstrap = 20, seed = 220, backend = "sparse", run_sdpd = FALSE)
  )
  td <- tidy(fit)
  expect_equal(nrow(td), 1)
  expect_true(is.finite(td$estimate))
})

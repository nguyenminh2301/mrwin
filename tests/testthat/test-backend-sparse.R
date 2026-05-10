test_that("sparse estimate matches dense estimate on the same data", {
  cfg <- mrwin_config(n_outcome = 100, m_snps = 8, seed = 42)
  dat <- mrwin_simulate(cfg)
  n_strata <- 4L

  strata_obj <- mrwin_prs_strata(dat$G, dat$true_betas, n_strata = n_strata)
  strata <- strata_obj$strata
  contrast_plan <- getFromNamespace(".mrwin_make_contrast_plan", "mrwin")(
    seq_len(n_strata), n_strata = n_strata
  )

  pair_kernels <- mrwin_precompute_pair_kernels(
    dat$time, dat$status, strata, contrast_plan
  )

  dense_fit <- mrwin_estimate(
    time = dat$time, status = dat$status,
    G = dat$G, X = dat$X, beta_hat = dat$true_betas,
    n_strata = n_strata
  )
  sparse_fit <- mrwin_sparse_estimate(
    time = dat$time, status = dat$status,
    G = dat$G, X = dat$X, beta_hat = dat$true_betas,
    n_strata = n_strata, pair_kernels = pair_kernels
  )

  expect_equal(sparse_fit$log_theta, dense_fit$log_theta, tolerance = 1e-12)
  expect_equal(sparse_fit$cwr, dense_fit$cwr, tolerance = 1e-12)
  expect_equal(sparse_fit$delta_x, dense_fit$delta_x, tolerance = 1e-12)
  expect_equal(sparse_fit$delta_isg, dense_fit$delta_isg, tolerance = 1e-12)
  expect_equal(sparse_fit$wins, dense_fit$wins, tolerance = 1e-12)
  expect_equal(sparse_fit$losses, dense_fit$losses, tolerance = 1e-12)
  expect_equal(sparse_fit$total, dense_fit$total, tolerance = 1e-12)
})

test_that("sparse estimate with weights matches dense estimate", {
  cfg <- mrwin_config(n_outcome = 80, m_snps = 6, seed = 7)
  dat <- mrwin_simulate(cfg)
  n_strata <- 3L
  weights <- runif(nrow(dat$G), 0.5, 2.0)

  strata_obj <- mrwin_prs_strata(dat$G, dat$true_betas, n_strata = n_strata)
  strata <- strata_obj$strata
  contrast_plan <- getFromNamespace(".mrwin_make_contrast_plan", "mrwin")(
    seq_len(n_strata), n_strata = n_strata
  )

  pair_kernels <- mrwin_precompute_pair_kernels(
    dat$time, dat$status, strata, contrast_plan
  )

  dense_fit <- mrwin_estimate(
    time = dat$time, status = dat$status,
    G = dat$G, X = dat$X, beta_hat = dat$true_betas,
    n_strata = n_strata, weights = weights
  )
  sparse_fit <- mrwin_sparse_estimate(
    time = dat$time, status = dat$status,
    G = dat$G, X = dat$X, beta_hat = dat$true_betas,
    n_strata = n_strata, weights = weights,
    pair_kernels = pair_kernels
  )

  expect_equal(sparse_fit$log_theta, dense_fit$log_theta, tolerance = 1e-12)
  expect_equal(sparse_fit$delta_x, dense_fit$delta_x, tolerance = 1e-12)
  expect_equal(sparse_fit$delta_isg, dense_fit$delta_isg, tolerance = 1e-12)
})

test_that("sparse bootstrap produces finite results", {
  cfg <- mrwin_config(n_outcome = 80, m_snps = 6, seed = 11)
  dat <- mrwin_simulate(cfg)

  boot <- mrwin_sparse_bootstrap(
    time = dat$time, status = dat$status,
    G = dat$G, X = dat$X,
    beta_hat = dat$true_betas,
    sigma_beta = rep(0.01, length(dat$true_betas)),
    n_strata = 4L, B = 20L, seed = 12
  )

  expect_s3_class(boot, "mrwin_bootstrap")
  expect_true(is.finite(boot$delta_gls))
  expect_true(is.finite(boot$se_delta_gls))
  expect_true(is.finite(boot$dscwr))
  expect_equal(nrow(boot$sigma_isg), 3)
  expect_equal(ncol(boot$sigma_isg), 3)
})

test_that("sparse bootstrap is reproducible with a fixed seed", {
  cfg <- mrwin_config(n_outcome = 80, m_snps = 6, seed = 33)
  dat <- mrwin_simulate(cfg)

  boot1 <- mrwin_sparse_bootstrap(
    time = dat$time, status = dat$status,
    G = dat$G, X = dat$X,
    beta_hat = dat$true_betas,
    sigma_beta = rep(0.01, length(dat$true_betas)),
    n_strata = 4L, B = 15L, seed = 34
  )
  boot2 <- mrwin_sparse_bootstrap(
    time = dat$time, status = dat$status,
    G = dat$G, X = dat$X,
    beta_hat = dat$true_betas,
    sigma_beta = rep(0.01, length(dat$true_betas)),
    n_strata = 4L, B = 15L, seed = 34
  )

  expect_equal(boot1$bootstrap_log_theta, boot2$bootstrap_log_theta)
  expect_equal(boot1$bootstrap_delta_x, boot2$bootstrap_delta_x)
  expect_equal(boot1$delta_gls, boot2$delta_gls)
  expect_equal(boot1$ci95_delta_fieller, boot2$ci95_delta_fieller)
})

test_that("sparse and dense backends produce matching mrwin() results", {
  cfg <- mrwin_config(n_outcome = 100, m_snps = 8, seed = 55)
  dat <- mrwin_simulate(cfg)

  ctrl_dense <- mrwin_controls(
    n_strata = 4L, bootstrap = 30L, seed = 56,
    backend = "dense", run_sdpd = FALSE
  )
  ctrl_sparse <- mrwin_controls(
    n_strata = 4L, bootstrap = 30L, seed = 56,
    backend = "sparse", run_sdpd = FALSE
  )

  fit_dense <- mrwin(
    endpoint = mrwin_endpoint(dat$time, dat$status, colnames(dat$time)),
    genotype = dat$G, exposure = dat$X,
    gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
    controls = ctrl_dense
  )
  fit_sparse <- mrwin(
    endpoint = mrwin_endpoint(dat$time, dat$status, colnames(dat$time)),
    genotype = dat$G, exposure = dat$X,
    gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
    controls = ctrl_sparse
  )

  expect_equal(fit_dense$point$delta_gls, fit_sparse$point$delta_gls, tolerance = 1e-10)
  expect_equal(fit_dense$point$dscwr, fit_sparse$point$dscwr, tolerance = 1e-10)
  expect_equal(fit_dense$point$log_theta, fit_sparse$point$log_theta, tolerance = 1e-10)
  expect_equal(fit_dense$point$delta_x, fit_sparse$point$delta_x, tolerance = 1e-10)
  expect_equal(fit_dense$point$delta_isg, fit_sparse$point$delta_isg, tolerance = 1e-10)
  expect_equal(fit_dense$inference$se_delta_gls, fit_sparse$inference$se_delta_gls, tolerance = 1e-10)
  expect_equal(fit_dense$heterogeneity$q, fit_sparse$heterogeneity$q, tolerance = 1e-10)
})

test_that("mrwin() blocks rcpp backend gracefully", {
  cfg <- mrwin_config(n_outcome = 50, m_snps = 4, seed = 99)
  dat <- mrwin_simulate(cfg)

  expect_error(
    mrwin(
      endpoint = mrwin_endpoint(dat$time, dat$status, colnames(dat$time)),
      genotype = dat$G, exposure = dat$X,
      gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
      controls = mrwin_controls(backend = "rcpp")
    ),
    "rcpp"
  )
})

test_that("precomputed pair kernels can be reused across calls", {
  cfg <- mrwin_config(n_outcome = 60, m_snps = 5, seed = 77)
  dat <- mrwin_simulate(cfg)
  n_strata <- 3L

  strata_obj <- mrwin_prs_strata(dat$G, dat$true_betas, n_strata = n_strata)
  contrast_plan <- getFromNamespace(".mrwin_make_contrast_plan", "mrwin")(
    seq_len(n_strata), n_strata = n_strata
  )

  pk <- mrwin_precompute_pair_kernels(
    dat$time, dat$status, strata_obj$strata, contrast_plan
  )

  expect_true(is.list(pk))
  expect_length(pk, nrow(contrast_plan))

  for (i in seq_along(pk)) {
    if (!is.null(pk[[i]])) {
      expect_true(is.matrix(pk[[i]]$kernel))
      expect_true(all(pk[[i]]$kernel %in% c(-1L, 0L, 1L)))
    }
  }
})

test_that("mrwin_verify_sparse_dense_parity reports all close", {
  result <- mrwin_verify_sparse_dense_parity(
    n = 80, m_snps = 6, n_strata = 4, B = 20, seed = 88
  )

  expect_s3_class(result, "mrwin_parity_check")
  expect_true(result$all_close)
  expect_true(all(vapply(result$checks, function(x) x < 1e-10, logical(1))))
})

test_that("benchmark runs without error on small configs", {
  bench <- mrwin_benchmark(
    n_vec = c(50, 100),
    m_snps = 5,
    n_strata = 3,
    B = 10,
    backends = c("dense", "sparse"),
    seed = 100,
    n_iter = 1
  )

  expect_s3_class(bench, "mrwin_benchmark")
  expect_equal(nrow(bench$results), 4)
  expect_true(all(bench$results$success))
  expect_true(all(bench$results$time_seconds >= 0))

  s <- summary(bench)
  expect_true(is.data.frame(s))
  expect_true(nrow(s) > 0)
})

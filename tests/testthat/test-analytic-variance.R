# WP17 / M3: analytic influence-function covariance (sampling part).

test_that("analytic covariance is symmetric PSD with correct dimensions", {
  cfg <- mrwin_config(n_outcome = 400, m_snps = 10, seed = 101)
  dat <- mrwin_simulate(cfg, seed = 101)
  est <- mrwin_estimate(dat$time, dat$status, dat$G, dat$X, dat$true_betas, n_strata = 5)
  cov_u <- mrwin_analytic_covariance(est$kernel, est$strata, dat$X, est$contrast_plan)

  dm1 <- nrow(est$contrast_plan)
  expect_equal(dim(cov_u), c(2L * dm1, 2L * dm1))
  expect_true(isSymmetric(unname(cov_u), tol = 1e-8))
  ev <- eigen(cov_u, symmetric = TRUE, only.values = TRUE)$values
  expect_gte(min(ev), -1e-8)
})

test_that("analytic inference matches the fixed-strata multiplier bootstrap", {
  cfg <- mrwin_config(n_outcome = 800, m_snps = 12, seed = 202)
  dat <- mrwin_simulate(cfg, seed = 202)
  est <- mrwin_estimate(dat$time, dat$status, dat$G, dat$X, dat$true_betas, n_strata = 5)
  ana <- mrwin_analytic_inference(est, dat$X)

  boot <- mrwin_multiplier_bootstrap(
    time = dat$time, status = dat$status, G = dat$G, X = dat$X,
    beta_hat = dat$true_betas, sigma_beta = 0, n_strata = 5, B = 3000,
    seed = 202, adjustment = "none"
  )
  # SE (the quantity M3 replaces) agrees within Monte-Carlo + first-order error.
  expect_equal(ana$se_delta_gls, boot$se_delta_gls, tolerance = 0.12)
  expect_true(is.finite(ana$delta_gls))
  # CI width is 2*z*se, so this re-confirms the variance match.
  expect_equal(
    unname(ana$ci95_delta[["high"]] - ana$ci95_delta[["low"]]),
    unname(boot$ci95_delta[["high"]] - boot$ci95_delta[["low"]]),
    tolerance = 0.12
  )
})

test_that("analytic sampling + GWAS resample matches the full bootstrap (sigma_beta>0)", {
  cfg <- mrwin_config(n_outcome = 700, m_snps = 14, seed = 303)
  dat <- mrwin_simulate(cfg, seed = 303)
  est <- mrwin_estimate(dat$time, dat$status, dat$G, dat$X, dat$true_betas, n_strata = 5)
  ana <- mrwin_analytic_inference(
    est, dat$X, G = dat$G, beta_hat = dat$true_betas,
    sigma_beta = 0.2, n_strata = 5, B_gwas = 3000, seed = 303
  )
  expect_true(ana$gwas_included)
  # GWAS term should be a real, non-trivial share of the variance here
  expect_gt(sum(diag(ana$cov_gwas)), 0)

  boot <- mrwin_multiplier_bootstrap(
    time = dat$time, status = dat$status, G = dat$G, X = dat$X,
    beta_hat = dat$true_betas, sigma_beta = 0.2, n_strata = 5, B = 3000,
    seed = 303, adjustment = "none"
  )
  expect_equal(ana$se_delta_gls, boot$se_delta_gls, tolerance = 0.12)
})

test_that("GWAS resample covariance is symmetric with correct dimensions", {
  cfg <- mrwin_config(n_outcome = 400, m_snps = 10, seed = 44)
  dat <- mrwin_simulate(cfg, seed = 44)
  est <- mrwin_estimate(dat$time, dat$status, dat$G, dat$X, dat$true_betas, n_strata = 4)
  cg <- mrwin_gwas_resample_covariance(
    est$kernel, dat$G, dat$X, dat$true_betas, sigma_beta = 0.15,
    n_strata = 4, contrast_plan = est$contrast_plan, B_gwas = 500, seed = 44
  )
  dm1 <- nrow(est$contrast_plan)
  expect_equal(dim(cg), c(2L * dm1, 2L * dm1))
  expect_true(isSymmetric(unname(cg), tol = 1e-8))
})

test_that("analytic inference requires the dense kernel", {
  est <- structure(list(kernel = NULL), class = "mrwin_estimate")
  expect_error(mrwin_analytic_inference(est, 1:3), "dense kernel")
})

test_that("mrwin(inference='analytic') runs end-to-end and methods work", {
  cfg <- mrwin_config(n_outcome = 500, m_snps = 12, seed = 77)
  dat <- mrwin_simulate(cfg, seed = 77)
  ep <- mrwin_endpoint(dat$time, dat$status, colnames(dat$time))
  gw <- mrwin_gwas(dat$true_betas, rep(0.05, length(dat$true_betas)))
  fit <- mrwin(
    endpoint = ep, genotype = dat$G, exposure = dat$X, gwas = gw,
    controls = mrwin_controls(n_strata = 5, bootstrap = 500, seed = 77,
                              inference = "analytic", run_sdpd = FALSE)
  )
  expect_s3_class(fit, "mrwin_fit")
  expect_true(is.finite(fit$point$delta_gls))
  expect_true(is.finite(fit$inference$se_delta_gls))
  expect_equal(fit$bootstrap$covariance_method, "analytic_sampling_plus_gwas_resample")
  expect_output(print(fit))
  expect_no_error(summary(fit))
  expect_no_error(tidy(fit))
})

test_that("mrwin(inference='analytic') agrees with the bootstrap on se", {
  cfg <- mrwin_config(n_outcome = 700, m_snps = 12, seed = 88)
  dat <- mrwin_simulate(cfg, seed = 88)
  ep <- mrwin_endpoint(dat$time, dat$status, colnames(dat$time))
  gw <- mrwin_gwas(dat$true_betas, rep(0.05, length(dat$true_betas)))
  ctrl <- function(inf) mrwin_controls(n_strata = 5, bootstrap = 2500, seed = 88,
                                       inference = inf, run_sdpd = FALSE)
  fit_a <- mrwin(endpoint = ep, genotype = dat$G, exposure = dat$X, gwas = gw, controls = ctrl("analytic"))
  fit_b <- mrwin(endpoint = ep, genotype = dat$G, exposure = dat$X, gwas = gw, controls = ctrl("bootstrap"))
  expect_equal(fit_a$inference$se_delta_gls, fit_b$inference$se_delta_gls, tolerance = 0.15)
})

test_that("inference='analytic' supports IPTW end-to-end", {
  cfg <- mrwin_config(n_outcome = 300, m_snps = 8, seed = 9)
  dat <- mrwin_simulate(cfg, seed = 9)
  ep <- mrwin_endpoint(dat$time, dat$status, colnames(dat$time))
  gw <- mrwin_gwas(dat$true_betas, rep(0.05, length(dat$true_betas)))
  fit <- mrwin(endpoint = ep, genotype = dat$G, exposure = dat$X, gwas = gw,
               covariates = cbind(z = as.numeric(scale(dat$U))),
               controls = mrwin_controls(n_strata = 4, bootstrap = 100, seed = 9,
                                         inference = "analytic", adjustment = "ordinal_iptw",
                                         run_sdpd = FALSE))
  expect_s3_class(fit, "mrwin_fit")
  expect_true(is.finite(fit$inference$se_delta_gls))
  expect_no_error(summary(fit))
})

test_that("IPTW analytic agrees with the IPTW bootstrap up to the estimated-weights gap", {
  cfg <- mrwin_config(n_outcome = 400, m_snps = 10, seed = 66)
  dat <- mrwin_simulate(cfg, seed = 66)
  ep <- mrwin_endpoint(dat$time, dat$status, colnames(dat$time))
  gw <- mrwin_gwas(dat$true_betas, rep(0, length(dat$true_betas)))  # sigma_beta=0 isolates the weighted IF
  Z <- cbind(u = as.numeric(scale(dat$U)))
  ctrl <- function(inf) mrwin_controls(n_strata = 5, bootstrap = 800, seed = 66, inference = inf,
                                       adjustment = "ordinal_iptw", run_sdpd = FALSE)
  fit_a <- mrwin(endpoint = ep, genotype = dat$G, exposure = dat$X, gwas = gw, covariates = Z, controls = ctrl("analytic"))
  fit_b <- mrwin(endpoint = ep, genotype = dat$G, exposure = dat$X, gwas = gw, covariates = Z, controls = ctrl("bootstrap"))
  # weights-fixed IF omits the (small) estimated-weights correction -> ~few % gap
  expect_equal(fit_a$inference$se_delta_gls, fit_b$inference$se_delta_gls, tolerance = 0.12)
})

test_that("analytic covariance flags zero win/loss contrasts", {
  # a degenerate kernel with no wins in a contrast triggers a clear error
  cfg <- mrwin_config(n_outcome = 200, m_snps = 8, seed = 5)
  dat <- mrwin_simulate(cfg, seed = 5)
  est <- mrwin_estimate(dat$time, dat$status, dat$G, dat$X, dat$true_betas, n_strata = 4)
  k0 <- est$kernel
  k0[k0 == 1L] <- 0L  # remove all wins
  expect_error(
    mrwin_analytic_covariance(k0, est$strata, dat$X, est$contrast_plan),
    "positive win and loss"
  )
})

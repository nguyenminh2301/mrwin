# WP16 / M2: doubly-ranked stratification (Tian, Burgess et al., 2023).

test_that("doubly-ranked assignment matches a hand-computed example", {
  prs <- c(10, 20, 30, 40, 50, 60)
  exposure <- c(5, 9, 1, 8, 2, 6)
  res <- mrwin_doubly_ranked_strata(prs, exposure, n_strata = 3)
  # pre-blocks {1,2,3} and {4,5,6} by PRS; within each, exposure rank = stratum
  expect_equal(res$strata, c(2L, 3L, 1L, 3L, 1L, 2L))
  expect_equal(res$score, prs)
})

test_that("mean exposure increases with stratum and instrument stays balanced", {
  set.seed(11)
  n <- 1000L
  prs <- rnorm(n)
  # exposure depends on prs (valid instrument) plus noise
  exposure <- 0.7 * prs + rnorm(n)
  res <- mrwin_doubly_ranked_strata(prs, exposure, n_strata = 10L)

  mean_exp <- tapply(exposure, res$strata, mean)
  expect_true(all(diff(mean_exp) > 0))          # strictly increasing exposure
  mean_prs <- tapply(prs, res$strata, mean)
  # instrument means are much more balanced across strata than exposure means
  expect_lt(diff(range(mean_prs)), 0.5 * diff(range(mean_exp)))
})

test_that("full blocks produce balanced stratum counts", {
  set.seed(5)
  n <- 500L  # 500 / 10 = 50 full blocks, no remainder
  prs <- rnorm(n)
  exposure <- rnorm(n)
  res <- mrwin_doubly_ranked_strata(prs, exposure, n_strata = 10L)
  counts <- as.integer(table(factor(res$strata, levels = 1:10)))
  expect_equal(counts, rep(50L, 10L))
})

test_that("assignment is deterministic and handles a remainder block", {
  set.seed(3)
  n <- 47L  # not a multiple of 10
  prs <- rnorm(n)
  exposure <- rnorm(n)
  a <- mrwin_doubly_ranked_strata(prs, exposure, n_strata = 10L)
  b <- mrwin_doubly_ranked_strata(prs, exposure, n_strata = 10L)
  expect_identical(a$strata, b$strata)
  expect_true(all(a$strata >= 1L & a$strata <= 10L))
})

test_that("doubly-ranked validates inputs", {
  expect_error(mrwin_doubly_ranked_strata(1:5, 1:4, 2), "same length")
  expect_error(mrwin_doubly_ranked_strata(1:5, 1:5, 1), "at least 2")
  expect_error(mrwin_doubly_ranked_strata(1:3, 1:3, 10), "cannot exceed")
  expect_error(
    mrwin_doubly_ranked_strata(c(1, 2, NA, 4), c(1, 2, 3, 4), 2),
    "must be finite"
  )
})

test_that(".mrwin_assign_strata dispatches and requires X for doubly-ranked", {
  set.seed(7)
  N <- 300L; M <- 5L
  G <- matrix(rbinom(N * M, 2, 0.3), N, M)
  beta <- rnorm(M, 0, 0.3)
  score <- as.numeric(G %*% beta)
  X <- 0.7 * scale(score)[, 1] + rnorm(N)

  assign_strata <- getFromNamespace(".mrwin_assign_strata", "mrwin")
  prs <- assign_strata(G, beta, X, 5L, "prs_rank")
  dr <- assign_strata(G, beta, X, 5L, "doubly_ranked")
  expect_identical(prs$strata, mrwin_prs_strata(G, beta, n_strata = 5L)$strata)
  expect_identical(dr$strata, mrwin_doubly_ranked_strata(score, X, n_strata = 5L)$strata)
  expect_false(identical(prs$strata, dr$strata))
  expect_error(
    assign_strata(G, beta, NULL, 5L, "doubly_ranked"),
    "`X` is required"
  )
})

test_that("stratification wires end-to-end through mrwin() on every backend", {
  set.seed(42)
  N <- 800L; M <- 6L
  G <- matrix(rbinom(N * M, 2, 0.3), N, M)
  beta <- rnorm(M, 0, 0.3)
  score <- as.numeric(G %*% beta)
  X <- 0.8 * scale(score)[, 1] + rnorm(N)
  lp <- 0.4 * X
  time <- matrix(rexp(N, rate = exp(lp - mean(lp))), N, 1)
  status <- matrix(1L, N, 1)

  s_dr <- mrwin_doubly_ranked_strata(score, X, n_strata = 5L)$strata

  run <- function(strat, backend = "dense", inference = "bootstrap") {
    ctl <- mrwin_controls(
      n_strata = 5L, bootstrap = 30L, seed = 1L, backend = backend,
      inference = inference, stratification = strat, run_sdpd = FALSE
    )
    suppressWarnings(mrwin(
      endpoint = mrwin_endpoint(time, status), genotype = G, exposure = X,
      beta_gwas = beta, se_gwas = rep(0, M), controls = ctl
    ))
  }

  f_prs <- run("prs_rank")
  f_dr <- run("doubly_ranked")
  expect_identical(
    as.integer(f_prs$bootstrap$strata),
    mrwin_prs_strata(G, beta, n_strata = 5L)$strata
  )
  expect_identical(as.integer(f_dr$bootstrap$strata), s_dr)
  expect_false(identical(f_prs$bootstrap$strata, f_dr$bootstrap$strata))

  for (b in c("dense", "sparse", "fast")) {
    ff <- run("doubly_ranked", backend = b)
    expect_identical(as.integer(ff$bootstrap$strata), s_dr)
  }
  fa <- run("doubly_ranked", inference = "analytic")
  expect_identical(as.integer(fa$bootstrap$strata), s_dr)
})

test_that("mrwin() flags doubly-ranked as invalid for the between-stratum estimand", {
  # Doubly-ranked balances the instrument across strata, so the adjacent-stratum
  # DS-CWR contrast is confounder-driven (external calibration: type-I ~1.0).
  # mrwin() must keep it runnable (reproducibility / future LACE work) but emit a
  # structured `doubly_ranked_invalid` warning. prs_rank must not warn.
  set.seed(1)
  N <- 400L; M <- 5L
  G <- matrix(rbinom(N * M, 2, 0.3), N, M)
  beta <- rnorm(M, 0, 0.3)
  X <- 0.8 * scale(as.numeric(G %*% beta))[, 1] + rnorm(N)
  time <- matrix(rexp(N, rate = exp(0.3 * X - mean(0.3 * X))), N, 1)
  status <- matrix(1L, N, 1)
  fit_of <- function(strat) {
    ctl <- mrwin_controls(n_strata = 4L, bootstrap = 20L, seed = 1L,
                          stratification = strat, run_sdpd = FALSE)
    mrwin(endpoint = mrwin_endpoint(time, status), genotype = G, exposure = X,
          beta_gwas = beta, se_gwas = rep(0, M), controls = ctl)
  }
  expect_warning(fit <- fit_of("doubly_ranked"), "doubly.ranked|instrument")
  codes <- vapply(fit$warnings, function(w) w$code, character(1))
  expect_true("doubly_ranked_invalid" %in% codes)
  f_prs <- fit_of("prs_rank")
  prs_codes <- vapply(f_prs$warnings, function(w) w$code, character(1))
  expect_false("doubly_ranked_invalid" %in% prs_codes)
})

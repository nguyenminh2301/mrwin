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

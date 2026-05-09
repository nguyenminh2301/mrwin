.wp2_valid_inputs <- function() {
  time <- matrix(c(
    5, 3,
    4, 2,
    6, 4,
    7, 5,
    8, 6,
    9, 7
  ), ncol = 2, byrow = TRUE)
  status <- matrix(c(
    1, 1,
    0, 1,
    1, 0,
    0, 0,
    1, 1,
    0, 1
  ), ncol = 2, byrow = TRUE)
  G <- matrix(c(
    0, 1,
    1, 0,
    2, 1,
    0, 2,
    1, 1,
    2, 0
  ), ncol = 2, byrow = TRUE)
  X <- c(0.1, 0.4, 0.9, 1.1, 1.4, 1.8)
  list(time = time, status = status, G = G, X = X)
}

test_that("mrwin_validate_data returns normalized validated data", {
  x <- .wp2_valid_inputs()
  validated <- mrwin_validate_data(
    endpoint = mrwin_endpoint(x$time, x$status),
    genotype = x$G,
    exposure = x$X,
    gwas = mrwin_gwas(c(0.2, 0.1), c(0.01, 0.02)),
    controls = mrwin_controls(n_strata = 3, bootstrap = 5)
  )

  expect_s3_class(validated, "mrwin_validated_data")
  expect_equal(nrow(validated$G), 6)
  expect_equal(length(validated$issues$errors), 0)
})

test_that("mrwin_validate_data catches endpoint value errors", {
  x <- .wp2_valid_inputs()
  x$status[1, 1] <- 2
  expect_error(
    mrwin_validate_data(
      endpoint = mrwin_endpoint(x$time, x$status),
      genotype = x$G,
      exposure = x$X,
      gwas = mrwin_gwas(c(0.2, 0.1), c(0.01, 0.02)),
      controls = mrwin_controls(n_strata = 3, bootstrap = 5)
    ),
    "status_not_binary"
  )

  x <- .wp2_valid_inputs()
  x$time[1, 1] <- -1
  expect_error(
    mrwin_validate_data(
      endpoint = mrwin_endpoint(x$time, x$status),
      genotype = x$G,
      exposure = x$X,
      gwas = mrwin_gwas(c(0.2, 0.1), c(0.01, 0.02)),
      controls = mrwin_controls(n_strata = 3, bootstrap = 5)
    ),
    "time_negative"
  )
})

test_that("mrwin_validate_data catches genotype and GWAS errors", {
  x <- .wp2_valid_inputs()
  x$G[, 2] <- 1
  expect_error(
    mrwin_validate_data(
      endpoint = mrwin_endpoint(x$time, x$status),
      genotype = x$G,
      exposure = x$X,
      gwas = mrwin_gwas(c(0.2, 0.1), c(0.01, 0.02)),
      controls = mrwin_controls(n_strata = 3, bootstrap = 5)
    ),
    "zero_variance_genotype"
  )

  x <- .wp2_valid_inputs()
  expect_error(
    mrwin_validate_data(
      endpoint = mrwin_endpoint(x$time, x$status),
      genotype = x$G,
      exposure = x$X,
      gwas = mrwin_gwas(c(0.2, 0.1, 0.3), c(0.01, 0.02, 0.03)),
      controls = mrwin_controls(n_strata = 3, bootstrap = 5)
    ),
    "gwas_beta_mismatch"
  )
})

test_that("mrwin_validate_data catches terminal order violations", {
  x <- .wp2_valid_inputs()
  x$time[1, ] <- c(3, 5)
  x$status[1, ] <- c(1, 1)

  expect_error(
    mrwin_validate_data(
      endpoint = mrwin_endpoint(x$time, x$status, c("death", "hf")),
      genotype = x$G,
      exposure = x$X,
      gwas = mrwin_gwas(c(0.2, 0.1), c(0.01, 0.02)),
      controls = mrwin_controls(n_strata = 3, bootstrap = 5)
    ),
    "terminal_event_order_violation"
  )
})

test_that("mrwin_validate_data can drop incomplete rows explicitly", {
  x <- .wp2_valid_inputs()
  x$X[2] <- NA_real_
  validated <- mrwin_validate_data(
    endpoint = mrwin_endpoint(x$time, x$status),
    genotype = x$G,
    exposure = x$X,
    gwas = mrwin_gwas(c(0.2, 0.1), c(0.01, 0.02)),
    controls = mrwin_controls(n_strata = 3, bootstrap = 5),
    complete_cases = "drop"
  )

  expect_equal(nrow(validated$G), 5)
  expect_equal(validated$issues$warnings[[1]]$code, "rows_dropped")
})

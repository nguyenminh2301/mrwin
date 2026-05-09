test_that("hierarchical kernel is antisymmetric with zero diagonal", {
  time <- matrix(c(
    5, 7,
    3, 8,
    6, 4
  ), ncol = 2, byrow = TRUE)
  status <- matrix(c(
    0, 1,
    1, 0,
    0, 1
  ), ncol = 2, byrow = TRUE)

  h <- mrwin_kernel(time, status)

  expect_equal(diag(h), rep(0L, 3))
  expect_equal(h, -t(h))
  expect_true(all(h %in% c(-1L, 0L, 1L)))
})

test_that("stratum win loss counts agree with kernel block", {
  h <- matrix(c(
    0, 1, -1,
    -1, 0, 1,
    1, -1, 0
  ), nrow = 3, byrow = TRUE)

  sums <- mrwin_stratum_win_loss(h, idx_high = c(1, 2), idx_low = 3)

  expect_equal(unname(sums["wins"]), 1)
  expect_equal(unname(sums["losses"]), 1)
  expect_equal(unname(sums["total"]), 2)
})

test_that("pair kernel matches the corresponding dense kernel block", {
  time <- matrix(c(
    5, 7, 9,
    3, 8, 6,
    6, 4, 5,
    8, 6, 7,
    2, 9, 8
  ), ncol = 3, byrow = TRUE)
  status <- matrix(c(
    0, 1, 1,
    1, 0, 1,
    0, 1, 0,
    1, 1, 1,
    1, 0, 0
  ), ncol = 3, byrow = TRUE)

  high <- c(2, 4)
  low <- c(1, 3, 5)
  dense <- mrwin_kernel(time, status)
  sparse <- mrwin_pair_kernel(
    time[high, , drop = FALSE],
    status[high, , drop = FALSE],
    time[low, , drop = FALSE],
    status[low, , drop = FALSE]
  )

  expect_equal(sparse, dense[high, low, drop = FALSE])
})

test_that("pair win loss matches dense weighted aggregation", {
  time <- matrix(c(
    5, 7,
    3, 8,
    6, 4,
    8, 6
  ), ncol = 2, byrow = TRUE)
  status <- matrix(c(
    0, 1,
    1, 0,
    0, 1,
    1, 1
  ), ncol = 2, byrow = TRUE)
  high <- c(2, 4)
  low <- c(1, 3)
  weights <- c(1.5, 0.8, 1.2, 2.0)
  dense <- mrwin_kernel(time, status)

  dense_sums <- mrwin_stratum_win_loss(dense, high, low, weights = weights)
  sparse_sums <- mrwin_pair_win_loss(
    time[high, , drop = FALSE],
    status[high, , drop = FALSE],
    time[low, , drop = FALSE],
    status[low, , drop = FALSE],
    weights_high = weights[high],
    weights_low = weights[low]
  )

  expect_equal(sparse_sums, dense_sums)
  expect_equal(
    mrwin_pair_log_cwr(
      time[high, , drop = FALSE],
      status[high, , drop = FALSE],
      time[low, , drop = FALSE],
      status[low, , drop = FALSE],
      weights_high = weights[high],
      weights_low = weights[low]
    ),
    mrwin_stratum_log_cwr(dense, high, low, weights = weights)
  )
})

test_that("sparse adjacent summaries match dense summaries by strata", {
  time <- matrix(c(
    5, 7,
    3, 8,
    6, 4,
    8, 6,
    2, 9,
    9, 3
  ), ncol = 2, byrow = TRUE)
  status <- matrix(c(
    0, 1,
    1, 0,
    0, 1,
    1, 1,
    1, 0,
    0, 1
  ), ncol = 2, byrow = TRUE)
  strata <- c(1, 2, 1, 3, 2, 3)
  weights <- c(1, 2, 1.5, 1, 0.5, 2.5)
  dense <- mrwin_kernel(time, status)
  sparse <- mrwin_sparse_adjacent_win_loss(time, status, strata, weights)

  for (row in seq_len(nrow(sparse))) {
    high <- sparse$high[row]
    low <- sparse$low[row]
    dense_sums <- mrwin_stratum_win_loss(
      dense,
      which(strata == high),
      which(strata == low),
      weights = weights
    )
    expect_equal(unname(as.numeric(sparse[row, c("wins", "losses", "total")])), unname(dense_sums))
  }
})

test_that("pair kernel validates incompatible inputs", {
  expect_error(
    mrwin_pair_kernel(
      matrix(1, nrow = 2, ncol = 2),
      matrix(1, nrow = 2, ncol = 2),
      matrix(1, nrow = 2, ncol = 3),
      matrix(1, nrow = 2, ncol = 3)
    ),
    "same number of priority columns"
  )
})

test_that("sparse adjacent summaries validate inputs before pairing", {
  time <- matrix(c(1, 2, 3, 4), ncol = 2)
  status <- matrix(c(0, 1, 0, 1), ncol = 2)

  expect_error(
    mrwin_sparse_adjacent_win_loss(time, status, c(1, NA)),
    "non-missing integer strata"
  )
  expect_error(
    mrwin_sparse_adjacent_win_loss(time, status, c(1, 2.5)),
    "non-missing integer strata"
  )
  expect_error(
    mrwin_sparse_adjacent_win_loss(time, matrix(c(0, 2, 0, 1), ncol = 2), c(1, 3)),
    "binary 0/1"
  )
})

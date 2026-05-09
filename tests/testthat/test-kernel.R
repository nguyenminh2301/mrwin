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

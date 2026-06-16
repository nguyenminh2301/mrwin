# WP13 / S1 parity: fast single-endpoint (K=1) win/loss == dense reference.

test_that("fast pair win/loss matches dense on a hand example", {
  time <- matrix(c(5, 3, 6, 8, 2), ncol = 1)
  status <- matrix(c(0, 1, 0, 1, 1), ncol = 1)
  high <- c(2, 4)
  low <- c(1, 3, 5)
  dense <- mrwin_pair_win_loss(
    time[high, , drop = FALSE], status[high, , drop = FALSE],
    time[low, , drop = FALSE], status[low, , drop = FALSE]
  )
  fast <- mrwin_fast_pair_win_loss(
    time[high, , drop = FALSE], status[high, , drop = FALSE],
    time[low, , drop = FALSE], status[low, , drop = FALSE]
  )
  expect_equal(fast, dense)
})

test_that("fast pair win/loss matches dense over random K=1 cohorts (heavy ties)", {
  set.seed(20260614)
  for (rep in seq_len(60)) {
    n_high <- sample(1:30, 1)
    n_low <- sample(1:30, 1)
    tie_support <- sample(c(2L, 4L, 8L, 50L), 1)
    weighted <- runif(1) < 0.5

    t_high <- matrix(sample(seq_len(tie_support), n_high, replace = TRUE), ncol = 1)
    t_low <- matrix(sample(seq_len(tie_support), n_low, replace = TRUE), ncol = 1)
    s_high <- matrix(rbinom(n_high, 1, 0.5), ncol = 1)
    s_low <- matrix(rbinom(n_low, 1, 0.5), ncol = 1)
    w_high <- if (weighted) rexp(n_high) else NULL
    w_low <- if (weighted) rexp(n_low) else NULL

    dense <- mrwin_pair_win_loss(t_high, s_high, t_low, s_low,
                                 weights_high = w_high, weights_low = w_low)
    fast <- mrwin_fast_pair_win_loss(t_high, s_high, t_low, s_low,
                                     weights_high = w_high, weights_low = w_low)
    expect_equal(fast, dense, tolerance = 1e-10)
  }
})

test_that("fast pair win/loss handles all-censored and all-event columns", {
  set.seed(7)
  t_high <- matrix(sample(1:5, 20, replace = TRUE), ncol = 1)
  t_low <- matrix(sample(1:5, 25, replace = TRUE), ncol = 1)
  for (pair in list(c(0, 0), c(1, 1), c(0, 1), c(1, 0))) {
    s_high <- matrix(pair[1], nrow = 20, ncol = 1)
    s_low <- matrix(pair[2], nrow = 25, ncol = 1)
    dense <- mrwin_pair_win_loss(t_high, s_high, t_low, s_low)
    fast <- mrwin_fast_pair_win_loss(t_high, s_high, t_low, s_low)
    expect_equal(fast, dense)
  }
})

test_that("fast win/loss is antisymmetric across the pair", {
  set.seed(99)
  t_high <- matrix(sample(1:6, 30, replace = TRUE), ncol = 1)
  t_low <- matrix(sample(1:6, 28, replace = TRUE), ncol = 1)
  s_high <- matrix(rbinom(30, 1, 0.5), ncol = 1)
  s_low <- matrix(rbinom(28, 1, 0.5), ncol = 1)
  hl <- mrwin_fast_pair_win_loss(t_high, s_high, t_low, s_low)
  lh <- mrwin_fast_pair_win_loss(t_low, s_low, t_high, s_high)
  expect_equal(unname(hl[["wins"]]), unname(lh[["losses"]]))
  expect_equal(unname(hl[["losses"]]), unname(lh[["wins"]]))
})

test_that("fast adjacent summaries match the sparse dense summaries", {
  set.seed(2024)
  for (rep in seq_len(40)) {
    n <- sample(6:150, 1)
    n_strata <- sample(c(2L, 5L, 10L), 1)
    tie_support <- sample(c(3L, 8L, n), 1)
    weighted <- runif(1) < 0.5
    time <- matrix(sample(seq_len(tie_support), n, replace = TRUE), ncol = 1)
    status <- matrix(rbinom(n, 1, 0.5), ncol = 1)
    strata <- sample(seq_len(n_strata), n, replace = TRUE)
    weights <- if (weighted) rexp(n) else NULL

    dense <- mrwin_sparse_adjacent_win_loss(time, status, strata, weights)
    fast <- mrwin_fast_adjacent_win_loss(time, status, strata, weights)
    # align by (high, low) and compare
    expect_equal(nrow(fast), nrow(dense))
    if (nrow(dense) > 0L) {
      key_d <- paste(dense$high, dense$low, sep = "-")
      key_f <- paste(fast$high, fast$low, sep = "-")
      ord <- match(key_d, key_f)
      expect_false(any(is.na(ord)))
      expect_equal(fast[ord, c("wins", "losses", "total")],
                   dense[, c("wins", "losses", "total")],
                   tolerance = 1e-10, ignore_attr = TRUE)
    }
  }
})

test_that("fast pair win/loss matches dense over random K=2 cohorts (heavy ties)", {
  set.seed(20260614)
  for (rep in seq_len(80)) {
    n_high <- sample(1:25, 1)
    n_low <- sample(1:25, 1)
    sup <- sample(c(2L, 3L, 5L, 30L), 1)
    weighted <- runif(1) < 0.5

    t_high <- matrix(sample(seq_len(sup), 2 * n_high, replace = TRUE), ncol = 2)
    t_low <- matrix(sample(seq_len(sup), 2 * n_low, replace = TRUE), ncol = 2)
    s_high <- matrix(rbinom(2 * n_high, 1, 0.5), ncol = 2)
    s_low <- matrix(rbinom(2 * n_low, 1, 0.5), ncol = 2)
    w_high <- if (weighted) rexp(n_high) else NULL
    w_low <- if (weighted) rexp(n_low) else NULL

    dense <- mrwin_pair_win_loss(t_high, s_high, t_low, s_low,
                                 weights_high = w_high, weights_low = w_low)
    fast <- mrwin_fast_pair_win_loss(t_high, s_high, t_low, s_low,
                                     weights_high = w_high, weights_low = w_low)
    expect_equal(fast, dense, tolerance = 1e-9)
  }
})

test_that("fast K=2 adjacent summaries match dense across strata", {
  set.seed(606)
  for (rep in seq_len(25)) {
    n <- sample(6:120, 1)
    n_strata <- sample(c(2L, 4L), 1)
    sup <- sample(c(3L, 8L), 1)
    weighted <- runif(1) < 0.5
    time <- matrix(sample(seq_len(sup), 2 * n, replace = TRUE), ncol = 2)
    status <- matrix(rbinom(2 * n, 1, 0.5), ncol = 2)
    strata <- sample(seq_len(n_strata), n, replace = TRUE)
    weights <- if (weighted) rexp(n) else NULL

    dense <- mrwin_sparse_adjacent_win_loss(time, status, strata, weights)
    fast <- mrwin_fast_adjacent_win_loss(time, status, strata, weights)
    expect_equal(nrow(fast), nrow(dense))
    if (nrow(dense) > 0L) {
      key_d <- paste(dense$high, dense$low, sep = "-")
      key_f <- paste(fast$high, fast$low, sep = "-")
      ord <- match(key_d, key_f)
      expect_false(any(is.na(ord)))
      expect_equal(fast[ord, c("wins", "losses", "total")],
                   dense[, c("wins", "losses", "total")],
                   tolerance = 1e-9, ignore_attr = TRUE)
    }
  }
})

test_that("fast pair win/loss matches dense over random K=3 cohorts (heavy ties)", {
  set.seed(20260614)
  for (rep in seq_len(60)) {
    n_high <- sample(1:20, 1)
    n_low <- sample(1:20, 1)
    sup <- sample(c(2L, 3L, 5L, 40L), 1)
    weighted <- runif(1) < 0.5

    t_high <- matrix(sample(seq_len(sup), 3 * n_high, replace = TRUE), ncol = 3)
    t_low <- matrix(sample(seq_len(sup), 3 * n_low, replace = TRUE), ncol = 3)
    s_high <- matrix(rbinom(3 * n_high, 1, 0.5), ncol = 3)
    s_low <- matrix(rbinom(3 * n_low, 1, 0.5), ncol = 3)
    w_high <- if (weighted) rexp(n_high) else NULL
    w_low <- if (weighted) rexp(n_low) else NULL

    dense <- mrwin_pair_win_loss(t_high, s_high, t_low, s_low,
                                 weights_high = w_high, weights_low = w_low)
    fast <- mrwin_fast_pair_win_loss(t_high, s_high, t_low, s_low,
                                     weights_high = w_high, weights_low = w_low)
    expect_equal(fast, dense, tolerance = 1e-9)
  }
})

test_that("fast K=3 adjacent summaries match dense across strata", {
  set.seed(909)
  for (rep in seq_len(15)) {
    n <- sample(8:90, 1)
    n_strata <- sample(c(2L, 4L), 1)
    sup <- sample(c(3L, 8L), 1)
    weighted <- runif(1) < 0.5
    time <- matrix(sample(seq_len(sup), 3 * n, replace = TRUE), ncol = 3)
    status <- matrix(rbinom(3 * n, 1, 0.5), ncol = 3)
    strata <- sample(seq_len(n_strata), n, replace = TRUE)
    weights <- if (weighted) rexp(n) else NULL

    dense <- mrwin_sparse_adjacent_win_loss(time, status, strata, weights)
    fast <- mrwin_fast_adjacent_win_loss(time, status, strata, weights)
    expect_equal(nrow(fast), nrow(dense))
    if (nrow(dense) > 0L) {
      key_d <- paste(dense$high, dense$low, sep = "-")
      ord <- match(key_d, paste(fast$high, fast$low, sep = "-"))
      expect_false(any(is.na(ord)))
      expect_equal(fast[ord, c("wins", "losses", "total")],
                   dense[, c("wins", "losses", "total")],
                   tolerance = 1e-9, ignore_attr = TRUE)
    }
  }
})

test_that("fast pair path rejects K>3 input", {
  expect_error(
    mrwin_fast_pair_win_loss(
      matrix(1, nrow = 2, ncol = 4), matrix(0, nrow = 2, ncol = 4),
      matrix(1, nrow = 2, ncol = 4), matrix(0, nrow = 2, ncol = 4)
    ),
    "K in \\{1, 2, 3\\}"
  )
})

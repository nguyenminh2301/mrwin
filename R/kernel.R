mrwin_kernel <- function(time, status, block_size = 4000L) {
  time <- as.matrix(time)
  status <- as.matrix(status)
  if (!all(dim(time) == dim(status))) {
    stop("`time` and `status` must have the same dimensions.", call. = FALSE)
  }
  if (any(!is.finite(time))) {
    stop("`time` must contain finite observed times.", call. = FALSE)
  }

  n <- nrow(time)
  k <- ncol(time)
  h <- matrix(0L, nrow = n, ncol = n)
  block_size <- max(1L, as.integer(block_size))

  starts <- seq.int(1L, n, by = block_size)
  for (a in starts) {
    b <- min(a + block_size - 1L, n)
    rows <- a:b
    nr <- length(rows)
    decided <- matrix(FALSE, nrow = nr, ncol = n)
    decided[cbind(seq_len(nr), rows)] <- TRUE

    for (priority in seq_len(k)) {
      tk_i <- matrix(time[rows, priority], nrow = nr, ncol = n)
      tk_j <- matrix(time[, priority], nrow = nr, ncol = n, byrow = TRUE)
      dk_i <- matrix(status[rows, priority], nrow = nr, ncol = n)
      dk_j <- matrix(status[, priority], nrow = nr, ncol = n, byrow = TRUE)

      i_wins <- dk_j == 1L & tk_i > tk_j
      i_loses <- dk_i == 1L & tk_j > tk_i
      new_win <- i_wins & !decided
      new_loss <- i_loses & !decided

      h_block <- h[rows, , drop = FALSE]
      h_block[new_win] <- 1L
      h_block[new_loss] <- -1L
      h[rows, ] <- h_block
      decided <- decided | new_win | new_loss
    }
  }

  diag(h) <- 0L
  h
}

mrwin_stratum_win_loss <- function(kernel, idx_high, idx_low, weights = NULL) {
  idx_high <- as.integer(idx_high)
  idx_low <- as.integer(idx_low)
  if (length(idx_high) == 0L || length(idx_low) == 0L) {
    return(c(wins = NA_real_, losses = NA_real_, total = NA_real_))
  }

  h_block <- kernel[idx_high, idx_low, drop = FALSE]
  if (is.null(weights)) {
    wins <- sum(h_block == 1L)
    losses <- sum(h_block == -1L)
    total <- length(idx_high) * length(idx_low)
  } else {
    weights <- as.numeric(weights)
    pair_weights <- outer(weights[idx_high], weights[idx_low], "*")
    wins <- sum(pair_weights * (h_block == 1L))
    losses <- sum(pair_weights * (h_block == -1L))
    total <- sum(pair_weights)
  }

  c(wins = as.numeric(wins), losses = as.numeric(losses), total = as.numeric(total))
}

mrwin_stratum_log_cwr <- function(kernel, idx_high, idx_low, weights = NULL, floor = 1e-12) {
  sums <- mrwin_stratum_win_loss(kernel, idx_high, idx_low, weights = weights)
  log(max(sums[["wins"]], floor) / max(sums[["losses"]], floor))
}

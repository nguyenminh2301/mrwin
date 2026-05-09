mrwin_kernel <- function(time, status, block_size = 4000L) {
  time <- as.matrix(time)
  status <- as.matrix(status)
  if (!all(dim(time) == dim(status))) {
    stop("`time` and `status` must have the same dimensions.", call. = FALSE)
  }
  if (any(!is.finite(time))) {
    stop("`time` must contain finite observed times.", call. = FALSE)
  }
  if (!all(status %in% c(0, 1))) {
    stop("`status` values must be binary 0/1.", call. = FALSE)
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

mrwin_pair_kernel <- function(time_high, status_high, time_low, status_low) {
  time_high <- as.matrix(time_high)
  status_high <- as.matrix(status_high)
  time_low <- as.matrix(time_low)
  status_low <- as.matrix(status_low)
  .mrwin_validate_pair_inputs(time_high, status_high, time_low, status_low)

  n_high <- nrow(time_high)
  n_low <- nrow(time_low)
  k <- ncol(time_high)
  h <- matrix(0L, nrow = n_high, ncol = n_low)
  decided <- matrix(FALSE, nrow = n_high, ncol = n_low)

  for (priority in seq_len(k)) {
    tk_high <- matrix(time_high[, priority], nrow = n_high, ncol = n_low)
    tk_low <- matrix(time_low[, priority], nrow = n_high, ncol = n_low, byrow = TRUE)
    dk_high <- matrix(status_high[, priority], nrow = n_high, ncol = n_low)
    dk_low <- matrix(status_low[, priority], nrow = n_high, ncol = n_low, byrow = TRUE)

    high_wins <- dk_low == 1L & tk_high > tk_low
    high_loses <- dk_high == 1L & tk_low > tk_high
    new_win <- high_wins & !decided
    new_loss <- high_loses & !decided

    h[new_win] <- 1L
    h[new_loss] <- -1L
    decided <- decided | new_win | new_loss
  }

  h
}

mrwin_pair_win_loss <- function(
    time_high,
    status_high,
    time_low,
    status_low,
    weights_high = NULL,
    weights_low = NULL
) {
  h <- mrwin_pair_kernel(time_high, status_high, time_low, status_low)
  if (is.null(weights_high) && is.null(weights_low)) {
    wins <- sum(h == 1L)
    losses <- sum(h == -1L)
    total <- length(h)
  } else {
    if (is.null(weights_high) || is.null(weights_low)) {
      stop("`weights_high` and `weights_low` must be supplied together.", call. = FALSE)
    }
    weights_high <- as.numeric(weights_high)
    weights_low <- as.numeric(weights_low)
    if (length(weights_high) != nrow(h) || length(weights_low) != ncol(h)) {
      stop("Pair weights must match the corresponding stratum sizes.", call. = FALSE)
    }
    if (any(!is.finite(weights_high)) || any(!is.finite(weights_low)) ||
        any(weights_high < 0) || any(weights_low < 0)) {
      stop("Pair weights must be finite and non-negative.", call. = FALSE)
    }
    pair_weights <- outer(weights_high, weights_low, "*")
    wins <- sum(pair_weights * (h == 1L))
    losses <- sum(pair_weights * (h == -1L))
    total <- sum(pair_weights)
  }

  c(wins = as.numeric(wins), losses = as.numeric(losses), total = as.numeric(total))
}

mrwin_pair_log_cwr <- function(
    time_high,
    status_high,
    time_low,
    status_low,
    weights_high = NULL,
    weights_low = NULL,
    floor = 1e-12
) {
  sums <- mrwin_pair_win_loss(
    time_high = time_high,
    status_high = status_high,
    time_low = time_low,
    status_low = status_low,
    weights_high = weights_high,
    weights_low = weights_low
  )
  log(max(sums[["wins"]], floor) / max(sums[["losses"]], floor))
}

mrwin_stratum_win_loss <- function(kernel, idx_high, idx_low, weights = NULL) {
  kernel <- as.matrix(kernel)
  idx_high <- as.integer(idx_high)
  idx_low <- as.integer(idx_low)
  if (length(idx_high) == 0L || length(idx_low) == 0L) {
    return(c(wins = NA_real_, losses = NA_real_, total = NA_real_))
  }
  if (any(idx_high < 1L | idx_high > nrow(kernel)) || any(idx_low < 1L | idx_low > ncol(kernel))) {
    stop("Stratum indices are outside kernel dimensions.", call. = FALSE)
  }

  h_block <- kernel[idx_high, idx_low, drop = FALSE]
  if (is.null(weights)) {
    wins <- sum(h_block == 1L)
    losses <- sum(h_block == -1L)
    total <- length(idx_high) * length(idx_low)
  } else {
    weights <- as.numeric(weights)
    if (length(weights) != nrow(kernel)) {
      stop("`weights` length must match the number of kernel rows.", call. = FALSE)
    }
    if (any(!is.finite(weights)) || any(weights < 0)) {
      stop("`weights` must be finite and non-negative.", call. = FALSE)
    }
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

mrwin_sparse_adjacent_win_loss <- function(time, status, strata, weights = NULL) {
  time <- as.matrix(time)
  status <- as.matrix(status)
  strata <- suppressWarnings(as.numeric(strata))
  if (!all(dim(time) == dim(status))) {
    stop("`time` and `status` must have the same dimensions.", call. = FALSE)
  }
  if (any(!is.finite(time))) {
    stop("`time` must contain finite observed times.", call. = FALSE)
  }
  if (!all(status %in% c(0, 1))) {
    stop("`status` values must be binary 0/1.", call. = FALSE)
  }
  if (length(strata) != nrow(time)) {
    stop("`strata` length must match the number of endpoint rows.", call. = FALSE)
  }
  if (any(is.na(strata)) || any(!is.finite(strata)) || any(strata != floor(strata))) {
    stop("`strata` must contain non-missing integer strata.", call. = FALSE)
  }
  strata <- as.integer(strata)
  if (is.null(weights)) {
    weights <- rep(1, length(strata))
  } else {
    weights <- as.numeric(weights)
    if (length(weights) != length(strata)) {
      stop("`weights` length must match `strata`.", call. = FALSE)
    }
  }
  if (any(!is.finite(weights)) || any(weights < 0)) {
    stop("`weights` must be finite and non-negative.", call. = FALSE)
  }

  levels <- sort(unique(strata))
  if (length(levels) < 2L) {
    return(data.frame(
      high = integer(),
      low = integer(),
      wins = numeric(),
      losses = numeric(),
      total = numeric()
    ))
  }

  rows <- vector("list", max(length(levels) - 1L, 0L))
  pos <- 1L
  for (d in levels[-1L]) {
    low <- d - 1L
    if (!low %in% levels) {
      next
    }
    idx_high <- which(strata == d)
    idx_low <- which(strata == low)
    sums <- mrwin_pair_win_loss(
      time[idx_high, , drop = FALSE],
      status[idx_high, , drop = FALSE],
      time[idx_low, , drop = FALSE],
      status[idx_low, , drop = FALSE],
      weights_high = weights[idx_high],
      weights_low = weights[idx_low]
    )
    rows[[pos]] <- data.frame(
      high = d,
      low = low,
      wins = sums[["wins"]],
      losses = sums[["losses"]],
      total = sums[["total"]]
    )
    pos <- pos + 1L
  }

  if (pos == 1L) {
    return(data.frame(
      high = integer(),
      low = integer(),
      wins = numeric(),
      losses = numeric(),
      total = numeric()
    ))
  }
  do.call(rbind, rows[seq_len(pos - 1L)])
}

.mrwin_validate_pair_inputs <- function(time_high, status_high, time_low, status_low) {
  if (!all(dim(time_high) == dim(status_high))) {
    stop("`time_high` and `status_high` must have the same dimensions.", call. = FALSE)
  }
  if (!all(dim(time_low) == dim(status_low))) {
    stop("`time_low` and `status_low` must have the same dimensions.", call. = FALSE)
  }
  if (ncol(time_high) != ncol(time_low)) {
    stop("High and low strata must have the same number of priority columns.", call. = FALSE)
  }
  if (any(!is.finite(time_high)) || any(!is.finite(time_low))) {
    stop("Pair-kernel times must be finite.", call. = FALSE)
  }
  if (!all(status_high %in% c(0, 1)) || !all(status_low %in% c(0, 1))) {
    stop("Pair-kernel status values must be binary 0/1.", call. = FALSE)
  }
  invisible(TRUE)
}

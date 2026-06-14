# WP13 / Phase II S1: subquadratic single-endpoint (K=1) win/loss.
#
# For a single endpoint the weighted win/loss counts between two strata
# factorise into weighted dominance counts that a descending-time sweep
# evaluates in O(N log N) instead of the dense O(n_high * n_low) cross-product,
# with bit-for-bit identical results. See inst/spec/wp13-fast-kernel.md and
# inst/spec/implementation-plan.md (S1).
#
# Reduction (high set H = stratum d, low set L = stratum d-1; strict times):
#   W   = sum_{j in L, status_j=1} w_j * ( sum_{i in H : t_i > t_j} w_i )
#   Lo  = sum_{i in H, status_i=1} w_i * ( sum_{j in L : t_j > t_i} w_j )
#   Tot = ( sum_{i in H} w_i ) * ( sum_{j in L} w_j )

mrwin_fast_pair_win_loss <- function(
    time_high,
    status_high,
    time_low,
    status_low,
    weights_high = NULL,
    weights_low = NULL
) {
  time_high <- as.matrix(time_high)
  status_high <- as.matrix(status_high)
  time_low <- as.matrix(time_low)
  status_low <- as.matrix(status_low)
  if (ncol(time_high) != 1L || ncol(time_low) != 1L) {
    stop(
      "`mrwin_fast_pair_win_loss` is the single-endpoint (K=1) fast path; use ",
      "backend = \"sparse\" for K>1 until the WP13 S2 hierarchical fast path lands.",
      call. = FALSE
    )
  }
  .mrwin_validate_pair_inputs(time_high, status_high, time_low, status_low)

  n_high <- nrow(time_high)
  n_low <- nrow(time_low)
  if (is.null(weights_high) && is.null(weights_low)) {
    wh <- rep(1, n_high)
    wl <- rep(1, n_low)
  } else {
    if (is.null(weights_high) || is.null(weights_low)) {
      stop("`weights_high` and `weights_low` must be supplied together.", call. = FALSE)
    }
    wh <- as.numeric(weights_high)
    wl <- as.numeric(weights_low)
    if (length(wh) != n_high || length(wl) != n_low) {
      stop("Pair weights must match the corresponding stratum sizes.", call. = FALSE)
    }
    if (any(!is.finite(wh)) || any(!is.finite(wl)) || any(wh < 0) || any(wl < 0)) {
      stop("Pair weights must be finite and non-negative.", call. = FALSE)
    }
  }

  th <- as.numeric(time_high[, 1L])
  tl <- as.numeric(time_low[, 1L])
  sh <- as.integer(status_high[, 1L])
  sl <- as.integer(status_low[, 1L])

  times <- c(th, tl)
  side_high <- c(rep(TRUE, n_high), rep(FALSE, n_low))
  status <- c(sh, sl)
  wts <- c(wh, wl)

  ord <- order(times, decreasing = TRUE)
  t_s <- times[ord]
  hi_s <- side_high[ord]
  st_s <- status[ord]
  w_s <- wts[ord]

  n <- length(t_s)
  if (n == 0L) {
    return(c(wins = 0, losses = 0, total = sum(wh) * sum(wl)))
  }
  # group id by equal time over the descending order: strictly-greater-time
  # accumulation = exclusive prefix across groups.
  grp <- cumsum(c(TRUE, t_s[-1L] != t_s[-n]))
  high_w <- ifelse(hi_s, w_s, 0)
  low_w <- ifelse(hi_s, 0, w_s)
  grp_high <- as.numeric(tapply(high_w, grp, sum))
  grp_low <- as.numeric(tapply(low_w, grp, sum))
  excl_high <- cumsum(grp_high) - grp_high
  excl_low <- cumsum(grp_low) - grp_low
  acc_high_entry <- excl_high[grp]
  acc_low_entry <- excl_low[grp]

  win_mask <- (!hi_s) & (st_s == 1L)   # low entry with event -> high wins
  loss_mask <- (hi_s) & (st_s == 1L)   # high entry with event -> high loses
  wins <- sum(w_s[win_mask] * acc_high_entry[win_mask])
  losses <- sum(w_s[loss_mask] * acc_low_entry[loss_mask])

  c(
    wins = as.numeric(wins),
    losses = as.numeric(losses),
    total = as.numeric(sum(wh) * sum(wl))
  )
}

mrwin_fast_adjacent_win_loss <- function(time, status, strata, weights = NULL) {
  time <- as.matrix(time)
  status <- as.matrix(status)
  if (ncol(time) != 1L) {
    stop(
      "`mrwin_fast_adjacent_win_loss` is the single-endpoint (K=1) fast path; ",
      "use `mrwin_sparse_adjacent_win_loss` for K>1.",
      call. = FALSE
    )
  }
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
  empty <- data.frame(
    high = integer(), low = integer(),
    wins = numeric(), losses = numeric(), total = numeric()
  )
  if (length(levels) < 2L) {
    return(empty)
  }

  rows <- vector("list", length(levels) - 1L)
  pos <- 1L
  for (d in levels[-1L]) {
    low <- d - 1L
    if (!low %in% levels) {
      next
    }
    idx_high <- which(strata == d)
    idx_low <- which(strata == low)
    sums <- mrwin_fast_pair_win_loss(
      time[idx_high, , drop = FALSE],
      status[idx_high, , drop = FALSE],
      time[idx_low, , drop = FALSE],
      status[idx_low, , drop = FALSE],
      weights_high = weights[idx_high],
      weights_low = weights[idx_low]
    )
    rows[[pos]] <- data.frame(
      high = d, low = low,
      wins = sums[["wins"]], losses = sums[["losses"]], total = sums[["total"]]
    )
    pos <- pos + 1L
  }

  if (pos == 1L) {
    return(empty)
  }
  do.call(rbind, rows[seq_len(pos - 1L)])
}

# Backend dispatch used by the sparse bootstrap: fast path when fast=TRUE and the
# endpoint is single-priority, otherwise the dense pair kernel. Both return an
# identical c(wins, losses, total), so swapping is parity-preserving by contract.
.mrwin_pair_win_loss_backend <- function(
    time_high, status_high, time_low, status_low,
    weights_high = NULL, weights_low = NULL, fast = FALSE
) {
  use_fast <- isTRUE(fast) && ncol(as.matrix(time_high)) == 1L
  if (use_fast) {
    mrwin_fast_pair_win_loss(
      time_high, status_high, time_low, status_low,
      weights_high = weights_high, weights_low = weights_low
    )
  } else {
    mrwin_pair_win_loss(
      time_high, status_high, time_low, status_low,
      weights_high = weights_high, weights_low = weights_low
    )
  }
}

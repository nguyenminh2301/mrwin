# WP13 / Phase II S1+S2: subquadratic hierarchical win/loss.
#
# K=1 (S1): the dense O(n_high * n_low) pair sweep is replaced by an O(N log N)
# descending-time sweep over weighted dominance counts.
# K=2 (S2): decompose into level-1 separations (a K=1 sweep) plus, for pairs tied
# at level 1, level-2 separations split over the four (status1) regimes; the two
# mixed regimes use a Fenwick-based 2D dominance count, O(N log N) overall.
# Both return win/loss/total identical to the dense kernel. K>=3 is not yet a
# fast path and falls back to the dense backend (see dev/wp13-fast-kernel.md).

# Single-endpoint weighted sweep on numeric column vectors. Returns c(wins, losses).
.mrwin_fast_1d <- function(th, sh, wh, tl, sl, wl) {
  n_high <- length(th)
  n_low <- length(tl)
  if (n_high == 0L || n_low == 0L) {
    return(c(0, 0))
  }
  times <- c(th, tl)
  side_high <- c(rep(TRUE, n_high), rep(FALSE, n_low))
  status <- c(as.integer(sh), as.integer(sl))
  wts <- c(wh, wl)

  ord <- order(times, decreasing = TRUE)
  t_s <- times[ord]
  hi_s <- side_high[ord]
  st_s <- status[ord]
  w_s <- wts[ord]
  n <- length(t_s)

  grp <- cumsum(c(TRUE, t_s[-1L] != t_s[-n]))
  high_w <- ifelse(hi_s, w_s, 0)
  low_w <- ifelse(hi_s, 0, w_s)
  grp_high <- as.numeric(tapply(high_w, grp, sum))
  grp_low <- as.numeric(tapply(low_w, grp, sum))
  excl_high <- cumsum(grp_high) - grp_high
  excl_low <- cumsum(grp_low) - grp_low
  acc_high_entry <- excl_high[grp]
  acc_low_entry <- excl_low[grp]

  win_mask <- (!hi_s) & (st_s == 1L)
  loss_mask <- (hi_s) & (st_s == 1L)
  wins <- sum(w_s[win_mask] * acc_high_entry[win_mask])
  losses <- sum(w_s[loss_mask] * acc_low_entry[loss_mask])
  c(wins, losses)
}

# Weighted 2D dominance count via a Fenwick tree.
#   sum_{i,j} aw_i bw_j over A x B with
#     x:  ax_i <= bx_j (x_le) or ax_i >= bx_j (!x_le)
#     y:  ay_i >  by_j (y_a_gt_b) or ay_i < by_j (!y_a_gt_b)   [strict]
.mrwin_dom2d_count <- function(ax, ay, aw, bx, by, bw, x_le, y_a_gt_b) {
  na <- length(ax)
  nb <- length(bx)
  if (na == 0L || nb == 0L) {
    return(0)
  }
  ys <- sort(unique(c(ay, by)))
  m <- length(ys)
  ay_r <- match(ay, ys)
  by_r <- match(by, ys)

  tree <- numeric(m + 1L)
  bit_add <- function(i, v) {
    while (i <= m) {
      tree[i] <<- tree[i] + v
      i <- i + bitwAnd(i, -i)
    }
  }
  bit_prefix <- function(i) {
    s <- 0
    while (i > 0L) {
      s <- s + tree[i]
      i <- i - bitwAnd(i, -i)
    }
    s
  }

  if (x_le) {
    ao <- order(ax)
    bo <- order(bx)
  } else {
    ao <- order(ax, decreasing = TRUE)
    bo <- order(bx, decreasing = TRUE)
  }
  axo <- ax[ao]; ayro <- as.integer(ay_r[ao]); awo <- aw[ao]
  bxo <- bx[bo]; byro <- as.integer(by_r[bo]); bwo <- bw[bo]

  total <- 0
  res <- 0
  ai <- 1L
  for (jj in seq_len(nb)) {
    bxj <- bxo[jj]
    if (x_le) {
      while (ai <= na && axo[ai] <= bxj) {
        bit_add(ayro[ai], awo[ai]); total <- total + awo[ai]; ai <- ai + 1L
      }
    } else {
      while (ai <= na && axo[ai] >= bxj) {
        bit_add(ayro[ai], awo[ai]); total <- total + awo[ai]; ai <- ai + 1L
      }
    }
    r <- byro[jj]
    s <- if (y_a_gt_b) total - bit_prefix(r) else bit_prefix(r - 1L)
    res <- res + bwo[jj] * s
  }
  res
}

# K=2 assembly. th/tl are N x 2 matrices; sh/sl are N x 2 status matrices.
.mrwin_fast_pair_2d <- function(th, sh, tl, sl, wh, wl) {
  th1 <- th[, 1L]; th2 <- th[, 2L]
  sh1 <- as.integer(sh[, 1L]); sh2 <- as.integer(sh[, 2L])
  tl1 <- tl[, 1L]; tl2 <- tl[, 2L]
  sl1 <- as.integer(sl[, 1L]); sl2 <- as.integer(sl[, 2L])

  base <- .mrwin_fast_1d(th1, sh1, wh, tl1, sl1, wl)
  W <- base[1L]
  Lo <- base[2L]

  # regime (high d1=0, low d1=0): always tie at 1 -> level-2 1D on subsets
  Hi <- which(sh1 == 0L)
  Lj <- which(sl1 == 0L)
  if (length(Hi) && length(Lj)) {
    s <- .mrwin_fast_1d(th2[Hi], sh2[Hi], wh[Hi], tl2[Lj], sl2[Lj], wl[Lj])
    W <- W + s[1L]; Lo <- Lo + s[2L]
  }

  # regime (high d1=1, low d1=1): tie at 1 iff t1 equal -> group by t1
  Hi <- which(sh1 == 1L)
  Lj <- which(sl1 == 1L)
  if (length(Hi) && length(Lj)) {
    for (v in intersect(unique(th1[Hi]), unique(tl1[Lj]))) {
      his <- Hi[th1[Hi] == v]
      ljs <- Lj[tl1[Lj] == v]
      if (length(his) && length(ljs)) {
        s <- .mrwin_fast_1d(th2[his], sh2[his], wh[his], tl2[ljs], sl2[ljs], wl[ljs])
        W <- W + s[1L]; Lo <- Lo + s[2L]
      }
    }
  }

  # regime (high d1=0, low d1=1): tie at 1 iff th1 <= tl1
  H <- which(sh1 == 0L)
  Bw <- which(sl1 == 1L & sl2 == 1L)
  W <- W + .mrwin_dom2d_count(th1[H], th2[H], wh[H], tl1[Bw], tl2[Bw], wl[Bw],
                              x_le = TRUE, y_a_gt_b = TRUE)
  Hl <- which(sh1 == 0L & sh2 == 1L)
  Bl <- which(sl1 == 1L)
  Lo <- Lo + .mrwin_dom2d_count(th1[Hl], th2[Hl], wh[Hl], tl1[Bl], tl2[Bl], wl[Bl],
                                x_le = TRUE, y_a_gt_b = FALSE)

  # regime (high d1=1, low d1=0): tie at 1 iff th1 >= tl1
  H <- which(sh1 == 1L)
  Bw <- which(sl1 == 0L & sl2 == 1L)
  W <- W + .mrwin_dom2d_count(th1[H], th2[H], wh[H], tl1[Bw], tl2[Bw], wl[Bw],
                              x_le = FALSE, y_a_gt_b = TRUE)
  Hl <- which(sh1 == 1L & sh2 == 1L)
  Bl <- which(sl1 == 0L)
  Lo <- Lo + .mrwin_dom2d_count(th1[Hl], th2[Hl], wh[Hl], tl1[Bl], tl2[Bl], wl[Bl],
                                x_le = FALSE, y_a_gt_b = FALSE)

  c(wins = as.numeric(W), losses = as.numeric(Lo), total = as.numeric(sum(wh) * sum(wl)))
}

# ---- K=3 fast path -------------------------------------------------------
# 1D strict weighted dominance: sum_{i,j} aw_i bw_j with a.x > b.x (gt) or < (lt).
.mrwin_dom1d_strict <- function(ax, aw, bx, bw, gt) {
  if (length(ax) == 0L || length(bx) == 0L) {
    return(0)
  }
  o <- order(ax)
  xs <- ax[o]
  cumw <- c(0, cumsum(aw[o]))
  total <- cumw[length(cumw)]
  if (gt) {
    k <- findInterval(bx, xs)                      # count xs <= bx
    s <- total - cumw[k + 1L]
  } else {
    k <- findInterval(bx, xs, left.open = TRUE)    # count xs < bx
    s <- cumw[k + 1L]
  }
  sum(bw * s)
}

# 3D weighted count via CDQ on dim1, reducing each merge to the 2D counter.
#   dim1: a1 <= b1 (r1='le') or a1 >= b1 (r1='ge')
#   dim2: a2 <= b2 (r2='le') or a2 >= b2 (r2='ge')
#   dim3: a3 >  b3 (s3='gt') or a3 <  b3 (s3='lt')  [strict]
.mrwin_dom3d_count <- function(a1, a2, a3, aw, b1, b2, b3, bw, r1, r2, s3) {
  na <- length(a1)
  nb <- length(b1)
  if (na == 0L || nb == 0L) {
    return(0)
  }
  key1 <- c(a1, b1)
  is_b <- c(rep(0L, na), rep(1L, nb))
  x2 <- c(a2, b2)
  x3 <- c(a3, b3)
  w <- c(aw, bw)
  ord <- if (r1 == "le") order(key1, is_b) else order(-key1, is_b)
  is_b <- is_b[ord]; x2 <- x2[ord]; x3 <- x3[ord]; w <- w[ord]

  x_le <- (r2 == "le")
  y_gt <- (s3 == "gt")
  res <- 0
  cdq <- function(lo, hi) {
    if (lo >= hi) {
      return(invisible())
    }
    mid <- (lo + hi) %/% 2L
    cdq(lo, mid)
    cdq(mid + 1L, hi)
    left <- lo:mid
    la <- left[is_b[left] == 0L]
    if (length(la) == 0L) {
      return(invisible())
    }
    right <- (mid + 1L):hi
    rb <- right[is_b[right] == 1L]
    if (length(rb) == 0L) {
      return(invisible())
    }
    res <<- res + .mrwin_dom2d_count(x2[la], x3[la], w[la], x2[rb], x3[rb], w[rb],
                                     x_le = x_le, y_a_gt_b = y_gt)
  }
  cdq(1L, na + nb)
  res
}

.mrwin_tie_kind <- function(di, dj) {
  if (di == 0L && dj == 0L) {
    "any"
  } else if (di == 0L && dj == 1L) {
    "le"
  } else if (di == 1L && dj == 0L) {
    "ge"
  } else {
    "eq"
  }
}

# Count over A x B with tie-constraints c1,c2 in {any,le,ge,eq} on dims 1,2 and
# strict s3 in {gt,lt} on dim 3. Reduces 'any'/'eq', then dispatches 1D/2D/3D.
.mrwin_level3_count <- function(ax1, ax2, ax3, aw, bx1, bx2, bx3, bw, c1, c2, s3) {
  if (length(ax1) == 0L || length(bx1) == 0L) {
    return(0)
  }
  if (c1 == "eq") {
    va <- split(seq_along(ax1), ax1)
    vb <- split(seq_along(bx1), bx1)
    tot <- 0
    for (nm in intersect(names(va), names(vb))) {
      ia <- va[[nm]]; jb <- vb[[nm]]
      tot <- tot + .mrwin_level3_count(
        ax1[ia], ax2[ia], ax3[ia], aw[ia], bx1[jb], bx2[jb], bx3[jb], bw[jb],
        "any", c2, s3)
    }
    return(tot)
  }
  if (c2 == "eq") {
    va <- split(seq_along(ax2), ax2)
    vb <- split(seq_along(bx2), bx2)
    tot <- 0
    for (nm in intersect(names(va), names(vb))) {
      ia <- va[[nm]]; jb <- vb[[nm]]
      tot <- tot + .mrwin_level3_count(
        ax1[ia], ax2[ia], ax3[ia], aw[ia], bx1[jb], bx2[jb], bx3[jb], bw[jb],
        c1, "any", s3)
    }
    return(tot)
  }
  dims <- list()
  if (c1 %in% c("le", "ge")) dims[[length(dims) + 1L]] <- list(d = 1L, c = c1)
  if (c2 %in% c("le", "ge")) dims[[length(dims) + 1L]] <- list(d = 2L, c = c2)
  if (length(dims) == 0L) {
    return(.mrwin_dom1d_strict(ax3, aw, bx3, bw, gt = (s3 == "gt")))
  }
  if (length(dims) == 1L) {
    d <- dims[[1L]]$d
    cc <- dims[[1L]]$c
    ax <- if (d == 1L) ax1 else ax2
    bx <- if (d == 1L) bx1 else bx2
    return(.mrwin_dom2d_count(ax, ax3, aw, bx, bx3, bw,
                              x_le = (cc == "le"), y_a_gt_b = (s3 == "gt")))
  }
  .mrwin_dom3d_count(ax1, ax2, ax3, aw, bx1, bx2, bx3, bw, c1, c2, s3)
}

# K=3 assembly. th/tl are N x 3 matrices; sh/sl are N x 3 status matrices.
.mrwin_fast_pair_3d <- function(th, sh, tl, sl, wh, wl) {
  base <- .mrwin_fast_pair_2d(th[, 1:2, drop = FALSE], sh[, 1:2, drop = FALSE],
                              tl[, 1:2, drop = FALSE], sl[, 1:2, drop = FALSE], wh, wl)
  W <- base[["wins"]]
  Lo <- base[["losses"]]
  tot <- base[["total"]]

  th1 <- th[, 1L]; th2 <- th[, 2L]; th3 <- th[, 3L]
  sh1 <- as.integer(sh[, 1L]); sh2 <- as.integer(sh[, 2L]); sh3 <- as.integer(sh[, 3L])
  tl1 <- tl[, 1L]; tl2 <- tl[, 2L]; tl3 <- tl[, 3L]
  sl1 <- as.integer(sl[, 1L]); sl2 <- as.integer(sl[, 2L]); sl3 <- as.integer(sl[, 3L])

  for (di1 in c(0L, 1L)) {
    for (di2 in c(0L, 1L)) {
      Hbase <- which(sh1 == di1 & sh2 == di2)
      if (length(Hbase) == 0L) next
      for (dj1 in c(0L, 1L)) {
        for (dj2 in c(0L, 1L)) {
          c1 <- .mrwin_tie_kind(di1, dj1)
          c2 <- .mrwin_tie_kind(di2, dj2)
          Lbase <- which(sl1 == dj1 & sl2 == dj2)
          if (length(Lbase) == 0L) next
          Lw <- Lbase[sl3[Lbase] == 1L]
          if (length(Lw)) {
            W <- W + .mrwin_level3_count(
              th1[Hbase], th2[Hbase], th3[Hbase], wh[Hbase],
              tl1[Lw], tl2[Lw], tl3[Lw], wl[Lw], c1, c2, "gt")
          }
          Hl <- Hbase[sh3[Hbase] == 1L]
          if (length(Hl)) {
            Lo <- Lo + .mrwin_level3_count(
              th1[Hl], th2[Hl], th3[Hl], wh[Hl],
              tl1[Lbase], tl2[Lbase], tl3[Lbase], wl[Lbase], c1, c2, "lt")
          }
        }
      }
    }
  }
  c(wins = as.numeric(W), losses = as.numeric(Lo), total = as.numeric(tot))
}

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
  k <- ncol(time_high)
  if (k > 3L || ncol(time_low) > 3L) {
    stop(
      "`mrwin_fast_pair_win_loss` supports K in {1, 2, 3}; use backend = \"sparse\" ",
      "for K>3 until the hierarchical fast path is extended further.",
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

  # Prefer the compiled C++ kernel; fall back to the pure-R reference if the
  # package was built without compilation.
  if (exists("mrwin_fast_pair_cpp", mode = "function")) {
    sh_i <- status_high; storage.mode(sh_i) <- "integer"
    sl_i <- status_low; storage.mode(sl_i) <- "integer"
    storage.mode(time_high) <- "double"
    storage.mode(time_low) <- "double"
    return(mrwin_fast_pair_cpp(time_high, sh_i, time_low, sl_i, wh, wl))
  }
  if (k == 1L) {
    s <- .mrwin_fast_1d(
      as.numeric(time_high[, 1L]), as.integer(status_high[, 1L]), wh,
      as.numeric(time_low[, 1L]), as.integer(status_low[, 1L]), wl
    )
    return(c(
      wins = as.numeric(s[1L]),
      losses = as.numeric(s[2L]),
      total = as.numeric(sum(wh) * sum(wl))
    ))
  }
  if (k == 2L) {
    return(.mrwin_fast_pair_2d(time_high, status_high, time_low, status_low, wh, wl))
  }
  .mrwin_fast_pair_3d(time_high, status_high, time_low, status_low, wh, wl)
}

mrwin_fast_adjacent_win_loss <- function(time, status, strata, weights = NULL) {
  time <- as.matrix(time)
  status <- as.matrix(status)
  if (ncol(time) > 3L) {
    stop(
      "`mrwin_fast_adjacent_win_loss` supports K in {1, 2, 3}; use ",
      "`mrwin_sparse_adjacent_win_loss` for K>3.",
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
# endpoint has K in {1, 2}, otherwise the dense pair kernel. Both return an
# identical c(wins, losses, total), so swapping is parity-preserving by contract.
.mrwin_pair_win_loss_backend <- function(
    time_high, status_high, time_low, status_low,
    weights_high = NULL, weights_low = NULL, fast = FALSE
) {
  use_fast <- isTRUE(fast) && ncol(as.matrix(time_high)) <= 3L
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

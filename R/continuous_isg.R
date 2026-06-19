# WP15 / M1: continuous kernel-smoothed instrument-standardised gradient (ISG).
#
# Replaces the arbitrary stratum count D with a bandwidth-controlled local
# gradient. At a centre rank c (in (0,1)) with bandwidth b and kernel K, define
# smooth upper/lower membership weights and the local win/loss + exposure
# contrast; pool the local contrasts with the SAME analytic influence-function
# covariance + Fieller machinery validated in R2 (the calibrated inference).
#
#   r_i = rank(PRS)/N
#   a_i(c) = K((r_i-c)/b) * 1{c <  r_i <= c+b}     (upper, half-open)
#   b_i(c) = K((r_i-c)/b) * 1{c-b < r_i <= c}      (lower)
#   W(c) = sum_ij a_i b_j 1{win};  L(c) = sum_ij a_i b_j 1{loss}
#   log_theta(c) = log(W/L);  Delta_X(c) = wmean(X|a) - wmean(X|b)
#   delta_isg(c) = log_theta(c) / Delta_X(c)
#
# Boxcar special case: K = boxcar, b = 1/D, centres c = k/D (k=1..D-1) reproduce
# the D-1 adjacent decile contrasts EXACTLY (the internal-consistency check).
#
# Influence function (multiplier xi_i = 1+e_i; each subject is upper OR lower at
# a given centre): C[t, c] = a_t (P+_t/W0 - P-_t/L0) + b_t (Q+_t/W0 - Q-_t/L0),
# C[t, m+c] = a_t (X_t - Xbar_a)/sum(a) - b_t (X_t - Xbar_b)/sum(b);
# cov_u = C^T C captures cross-centre correlation via shared subjects.

.mrwin_isg_kernel_fn <- function(name) {
  switch(name,
    boxcar = function(u) as.numeric(abs(u) <= 1),
    epanechnikov = function(u) pmax(0, 0.75 * (1 - u^2)),
    triangular = function(u) pmax(0, 1 - abs(u)),
    stop("unknown kernel: ", name, call. = FALSE)
  )
}

mrwin_continuous_isg <- function(kernel, score, X, weights = NULL,
                                 bandwidth = 0.1, n_centers = 25L,
                                 kernel_fn = c("epanechnikov", "boxcar", "triangular"),
                                 centers = NULL, floor = 1e-12) {
  kernel_fn <- match.arg(kernel_fn)
  kernel <- as.matrix(kernel)
  score <- as.numeric(score)
  X <- as.numeric(X)
  n <- length(score)
  if (nrow(kernel) != n || ncol(kernel) != n) {
    stop("`kernel` must be an N x N matrix matching `score`.", call. = FALSE)
  }
  if (length(X) != n) stop("`X` length must match `score`.", call. = FALSE)
  if (is.null(weights)) weights <- rep(1, n) else weights <- as.numeric(weights)
  b <- as.numeric(bandwidth)
  Kf <- .mrwin_isg_kernel_fn(kernel_fn)

  # Work in integer rank-position space to avoid floating-point boundary errors
  # (so the boxcar/decile reduction is exact). Centres/bandwidth are quantiles in
  # (0,1); convert to rank units by * n. Half-open bins: strict lower, closed upper.
  R <- rank(score, ties.method = "first")                # integer 1..n
  if (is.null(centers)) {
    centers <- seq(b, 1 - b, length.out = n_centers)     # keep windows interior
  }
  m <- length(centers)
  h <- b * n
  pc <- centers * n
  win <- (kernel == 1L) + 0
  loss <- (kernel == -1L) + 0

  log_theta <- delta_x <- rep(NA_real_, m)
  C <- matrix(0, nrow = n, ncol = 2L * m)
  for (k in seq_len(m)) {
    p0 <- pc[k]
    kw <- Kf((R - p0) / h) * weights
    a <- kw * (R > p0 & R <= p0 + h)        # upper (half-open, closed at top)
    bl <- kw * (R > p0 - h & R <= p0)        # lower (closed at p0)
    sa <- sum(a); sbl <- sum(bl)
    if (sa <= 0 || sbl <= 0) next
    W0 <- drop(crossprod(a, win %*% bl))
    L0 <- drop(crossprod(a, loss %*% bl))
    if (W0 <= 0 || L0 <= 0) next
    log_theta[k] <- log(W0 / L0)
    mxa <- sum(a * X) / sa
    mxb <- sum(bl * X) / sbl
    delta_x[k] <- mxa - mxb
    Pw_plus <- as.numeric(win %*% bl)              # subject as upper member (i)
    Pw_minus <- as.numeric(loss %*% bl)
    Qw_plus <- as.numeric(crossprod(win, a))       # subject as lower member (j)
    Qw_minus <- as.numeric(crossprod(loss, a))
    C[, k] <- a * (Pw_plus / W0 - Pw_minus / L0) + bl * (Qw_plus / W0 - Qw_minus / L0)
    C[, m + k] <- a * (X - mxa) / sa - bl * (X - mxb) / sbl
  }

  ok <- which(!is.na(log_theta) & abs(delta_x) > floor)
  if (length(ok) < 1L) {
    stop("No usable centres (all degenerate or Delta_X ~ 0); widen the bandwidth.", call. = FALSE)
  }
  cols <- c(ok, m + ok)
  cov_u <- crossprod(C[, cols, drop = FALSE])
  cov_u <- 0.5 * (cov_u + t(cov_u))
  lt <- log_theta[ok]; dx <- delta_x[ok]
  nm <- paste0("c", round(centers[ok], 4))
  names(lt) <- names(dx) <- nm
  dimnames(cov_u) <- list(c(paste0("log_theta:", nm), paste0("delta_x:", nm)),
                          c(paste0("log_theta:", nm), paste0("delta_x:", nm)))

  delta_isg <- lt / dx
  sigma_isg <- .mrwin_isg_covariance(cov_u, lt, dx, floor = floor)
  pooled <- mrwin_gls_pool(delta_isg, sigma_isg, shrink = TRUE)
  fieller <- .mrwin_fieller_ci(cov_u, sigma_isg, lt, dx, z = 1.96)

  list(
    bandwidth = b, kernel_fn = kernel_fn, n_centers = length(ok),
    centers = centers[ok], log_theta = lt, delta_x = dx, delta_isg = delta_isg,
    delta_gls = pooled$delta_gls, se_delta_gls = pooled$se_delta_gls,
    dscwr = pooled$dscwr, ci95_delta = pooled$ci95_delta,
    ci95_delta_fieller = fieller$delta, ci95_dscwr_fieller = exp(pmax(pmin(fieller$delta, 50), -50)),
    fieller_unbounded = fieller$unbounded, fieller_p_value = fieller$p_value,
    q = pooled$q, q_df = pooled$q_df, q_p_value = pooled$q_p_value,
    cov_u = cov_u, sigma_isg = sigma_isg
  )
}

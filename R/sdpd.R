mrwin_aalen_per_snp <- function(G, time, status) {
  G <- as.matrix(G)
  time <- as.numeric(time)
  status <- as.numeric(status)
  m <- ncol(G)
  beta <- se <- rep(NA_real_, m)

  for (j in seq_len(m)) {
    design <- cbind(time, G[, j] * time)
    xtx <- crossprod(design)
    inv <- tryCatch(solve(xtx), error = function(e) NULL)
    if (is.null(inv)) {
      beta[j] <- 0
      se[j] <- Inf
      next
    }
    coef <- inv %*% crossprod(design, status)
    resid <- status - drop(design %*% coef)
    meat <- crossprod(design * resid, design * resid)
    vcov <- inv %*% meat %*% inv
    beta[j] <- coef[2L, 1L]
    se[j] <- sqrt(max(vcov[2L, 2L], 1e-20))
  }

  list(beta = beta, se = se)
}

mrwin_cox_per_snp <- function(G, time, status, max_iter = 25L, tol = 1e-7) {
  G <- as.matrix(G)
  time <- as.numeric(time)
  status <- as.integer(status)
  ord <- order(time)
  time_s <- time[ord]
  status_s <- status[ord]
  G_s <- G[ord, , drop = FALSE]
  event_idx <- which(status_s == 1L)
  m <- ncol(G_s)

  if (length(event_idx) == 0L) {
    return(list(beta = rep(0, m), se = rep(NA_real_, m)))
  }

  beta <- rep(0, m)
  for (iter in seq_len(max_iter)) {
    lp <- sweep(G_s, 2L, beta, "*")
    risk <- exp(lp)
    cum_w <- .mrwin_tail_cumsum(risk)
    cum_gw <- .mrwin_tail_cumsum(G_s * risk)
    cum_g2w <- .mrwin_tail_cumsum(G_s * G_s * risk)

    mean_g <- cum_gw[event_idx, , drop = FALSE] /
      pmax(cum_w[event_idx, , drop = FALSE], 1e-12)
    score <- colSums(G_s[event_idx, , drop = FALSE] - mean_g)
    var_g <- cum_g2w[event_idx, , drop = FALSE] /
      pmax(cum_w[event_idx, , drop = FALSE], 1e-12) - mean_g^2
    info <- colSums(var_g)
    step <- pmax(pmin(score / pmax(info, 1e-10), 0.5), -0.5)
    beta_new <- beta + step
    if (max(abs(step)) < tol) {
      beta <- beta_new
      break
    }
    beta <- beta_new
  }

  lp <- sweep(G_s, 2L, beta, "*")
  risk <- exp(lp)
  cum_w <- .mrwin_tail_cumsum(risk)
  cum_gw <- .mrwin_tail_cumsum(G_s * risk)
  cum_g2w <- .mrwin_tail_cumsum(G_s * G_s * risk)
  mean_g <- cum_gw[event_idx, , drop = FALSE] /
    pmax(cum_w[event_idx, , drop = FALSE], 1e-12)
  var_g <- cum_g2w[event_idx, , drop = FALSE] /
    pmax(cum_w[event_idx, , drop = FALSE], 1e-12) - mean_g^2
  info <- colSums(var_g)
  list(beta = beta, se = 1 / sqrt(pmax(info, 1e-10)))
}

mrwin_mr_egger <- function(beta_x, beta_y, se_y) {
  beta_x <- as.numeric(beta_x)
  beta_y <- as.numeric(beta_y)
  se_y <- as.numeric(se_y)
  keep <- is.finite(beta_x) & is.finite(beta_y) & is.finite(se_y) & se_y > 0
  beta_x <- beta_x[keep]
  beta_y <- beta_y[keep]
  se_y <- se_y[keep]
  if (length(beta_x) < 3L) {
    stop("MR-Egger requires at least 3 valid SNPs.", call. = FALSE)
  }

  w <- 1 / pmax(se_y^2, 1e-12)
  design <- cbind(intercept = 1, beta_x = beta_x)
  xw <- design * sqrt(w)
  yw <- beta_y * sqrt(w)
  xtx <- crossprod(xw)
  coef <- solve(xtx, crossprod(xw, yw))
  resid <- beta_y - drop(design %*% coef)
  q <- sum(w * resid^2)
  phi <- max(1, q / max(length(beta_x) - 2L, 1L))
  vcov <- phi * solve(xtx)
  se <- sqrt(diag(vcov))

  list(
    intercept = coef[1L, 1L],
    intercept_se = se[1L],
    intercept_z = coef[1L, 1L] / se[1L],
    intercept_p_value = 2 * stats::pnorm(-abs(coef[1L, 1L] / se[1L])),
    slope = coef[2L, 1L],
    slope_se = se[2L],
    q = q,
    phi = phi,
    m = length(beta_x)
  )
}

mrwin_pleiotropy_bounded_ci <- function(delta_hat, se_delta, bias_radius, level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  sampling <- delta_hat + c(-z, z) * se_delta
  bounded <- c(sampling[1L] - abs(bias_radius), sampling[2L] + abs(bias_radius))
  list(
    delta_hat = delta_hat,
    se_delta = se_delta,
    bias_radius = abs(bias_radius),
    ci_delta_sampling = sampling,
    ci_delta_pleiotropy_bounded = bounded,
    dscwr = exp(delta_hat),
    ci_dscwr_sampling = exp(sampling),
    ci_dscwr_pleiotropy_bounded = exp(bounded),
    null_crossing_sampling = sampling[1L] <= 0 && sampling[2L] >= 0,
    null_crossing_pleiotropy_bounded = bounded[1L] <= 0 && bounded[2L] >= 0
  )
}

.mrwin_tail_cumsum <- function(x) {
  x <- as.matrix(x)
  out <- apply(x, 2L, function(col) rev(cumsum(rev(col))))
  matrix(out, nrow = nrow(x), ncol = ncol(x))
}

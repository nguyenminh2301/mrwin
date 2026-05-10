mrwin_sdpd <- function(
    G,
    X,
    time,
    status,
    scale = c("aalen", "cox", "both"),
    alpha = 0.05,
    min_snps = 10L
) {
  scale <- match.arg(scale)
  checked <- .mrwin_validate_sdpd_inputs(G, time, status, X = X)
  alpha <- .mrwin_validate_sdpd_alpha(alpha)
  min_snps <- suppressWarnings(as.integer(min_snps))
  if (length(min_snps) != 1L || is.na(min_snps) || min_snps < 3L) {
    stop("`min_snps` must be a single integer at least 3.", call. = FALSE)
  }

  G <- checked$G
  X <- checked$X
  time <- checked$time
  status <- checked$status
  beta_x <- .mrwin_lm_per_snp(G, X)
  scales <- if (scale == "both") c("aalen", "cox") else scale

  results <- lapply(scales, function(one_scale) {
    outcome <- if (one_scale == "aalen") {
      mrwin_aalen_per_snp(G, time, status)
    } else {
      mrwin_cox_per_snp(G, time, status)
    }
    egger <- tryCatch(
      mrwin_mr_egger(beta_x$beta, outcome$beta, outcome$se),
      error = function(e) list(error = conditionMessage(e))
    )
    if (is.null(egger$error)) {
      egger$rejected <- isTRUE(egger$intercept_p_value < alpha)
      egger$underpowered <- isTRUE(egger$m < min_snps)
    } else {
      egger$rejected <- NA
      egger$underpowered <- TRUE
    }
    egger$scale <- one_scale
    egger$outcome_summary <- outcome
    egger
  })
  names(results) <- scales

  rejected <- any(vapply(results, function(x) isTRUE(x$rejected), logical(1)))
  underpowered <- any(vapply(results, function(x) isTRUE(x$underpowered), logical(1)))

  structure(
    list(
      alpha = alpha,
      scale = scale,
      min_snps = min_snps,
      n = nrow(G),
      m_snps = ncol(G),
      n_events = sum(status == 1L),
      exposure_summary = beta_x,
      results = results,
      rejected = rejected,
      underpowered = underpowered
    ),
    class = c("mrwin_sdpd", "list")
  )
}

mrwin_aalen_per_snp <- function(G, time, status) {
  checked <- .mrwin_validate_sdpd_inputs(G, time, status)
  G <- checked$G
  time <- checked$time
  status <- as.numeric(checked$status)
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

  names(beta) <- names(se) <- colnames(G)
  structure(list(beta = beta, se = se, scale = "aalen"), class = c("mrwin_snp_summary", "list"))
}

mrwin_cox_per_snp <- function(G, time, status, max_iter = 25L, tol = 1e-7) {
  checked <- .mrwin_validate_sdpd_inputs(G, time, status)
  G <- checked$G
  time <- checked$time
  status <- checked$status
  max_iter <- suppressWarnings(as.integer(max_iter))
  if (length(max_iter) != 1L || is.na(max_iter) || max_iter < 1L) {
    stop("`max_iter` must be a positive integer.", call. = FALSE)
  }
  tol <- as.numeric(tol)
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("`tol` must be a positive finite number.", call. = FALSE)
  }
  ord <- order(time)
  time_s <- time[ord]
  status_s <- status[ord]
  G_s <- G[ord, , drop = FALSE]
  event_idx <- which(status_s == 1L)
  m <- ncol(G_s)

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
  se <- 1 / sqrt(pmax(info, 1e-10))
  names(beta) <- names(se) <- colnames(G_s)
  structure(list(beta = beta, se = se, scale = "cox"), class = c("mrwin_snp_summary", "list"))
}

mrwin_mr_egger <- function(beta_x, beta_y, se_y) {
  beta_x <- as.numeric(beta_x)
  beta_y <- as.numeric(beta_y)
  se_y <- as.numeric(se_y)
  if (length(beta_x) != length(beta_y) || length(beta_y) != length(se_y)) {
    stop("`beta_x`, `beta_y`, and `se_y` must have the same length.", call. = FALSE)
  }
  keep <- is.finite(beta_x) & is.finite(beta_y) & is.finite(se_y) & se_y > 0
  beta_x <- beta_x[keep]
  beta_y <- beta_y[keep]
  se_y <- se_y[keep]
  if (length(beta_x) < 3L) {
    stop("MR-Egger requires at least 3 valid SNPs.", call. = FALSE)
  }
  if (!is.finite(stats::var(beta_x)) || stats::var(beta_x) <= 0) {
    stop("MR-Egger requires positive variance in `beta_x`.", call. = FALSE)
  }

  w <- 1 / pmax(se_y^2, 1e-12)
  design <- cbind(intercept = 1, beta_x = beta_x)
  xw <- design * sqrt(w)
  yw <- beta_y * sqrt(w)
  xtx <- crossprod(xw)
  coef <- tryCatch(solve(xtx, crossprod(xw, yw)), error = function(e) NULL)
  if (is.null(coef)) {
    stop("MR-Egger design matrix is singular.", call. = FALSE)
  }
  resid <- beta_y - drop(design %*% coef)
  q <- sum(w * resid^2)
  phi <- max(1, q / max(length(beta_x) - 2L, 1L))
  xtx_inv <- tryCatch(solve(xtx), error = function(e) NULL)
  if (is.null(xtx_inv)) {
    stop("MR-Egger design matrix is singular.", call. = FALSE)
  }
  vcov <- phi * xtx_inv
  se <- sqrt(diag(vcov))

  structure(list(
    intercept = coef[1L, 1L],
    intercept_se = se[1L],
    intercept_z = coef[1L, 1L] / se[1L],
    intercept_p_value = 2 * stats::pnorm(-abs(coef[1L, 1L] / se[1L])),
    slope = coef[2L, 1L],
    slope_se = se[2L],
    q = q,
    phi = phi,
    m = length(beta_x)
  ), class = c("mrwin_mr_egger", "list"))
}

mrwin_pleiotropy_bounded_ci <- function(delta_hat, se_delta, bias_radius, level = 0.95) {
  delta_hat <- as.numeric(delta_hat)
  se_delta <- as.numeric(se_delta)
  bias_radius <- as.numeric(bias_radius)
  level <- as.numeric(level)
  if (length(delta_hat) != 1L || !is.finite(delta_hat)) {
    stop("`delta_hat` must be a single finite number.", call. = FALSE)
  }
  if (length(se_delta) != 1L || !is.finite(se_delta) || se_delta < 0) {
    stop("`se_delta` must be a single non-negative finite number.", call. = FALSE)
  }
  if (length(bias_radius) != 1L || !is.finite(bias_radius) || bias_radius < 0) {
    stop("`bias_radius` must be a single non-negative finite number.", call. = FALSE)
  }
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    stop("`level` must be between 0 and 1.", call. = FALSE)
  }
  z <- stats::qnorm(1 - (1 - level) / 2)
  sampling <- delta_hat + c(-z, z) * se_delta
  bounded <- c(sampling[1L] - abs(bias_radius), sampling[2L] + abs(bias_radius))
  names(sampling) <- names(bounded) <- c("low", "high")
  structure(list(
    delta_hat = delta_hat,
    se_delta = se_delta,
    bias_radius = abs(bias_radius),
    level = level,
    ci_delta_sampling = sampling,
    ci_delta_pleiotropy_bounded = bounded,
    dscwr = exp(delta_hat),
    ci_dscwr_sampling = exp(sampling),
    ci_dscwr_pleiotropy_bounded = exp(bounded),
    null_crossing_sampling = sampling[1L] <= 0 && sampling[2L] >= 0,
    null_crossing_pleiotropy_bounded = bounded[1L] <= 0 && bounded[2L] >= 0
  ), class = c("mrwin_pleiotropy_bounded_ci", "list"))
}

.mrwin_validate_sdpd_inputs <- function(G, time, status, X = NULL) {
  G <- as.matrix(G)
  storage.mode(G) <- "double"
  if (nrow(G) < 3L || ncol(G) < 1L) {
    stop("`G` must have at least 3 rows and 1 SNP column.", call. = FALSE)
  }
  if (any(!is.finite(G))) {
    stop("`G` must contain finite numeric values.", call. = FALSE)
  }
  zero_var <- vapply(seq_len(ncol(G)), function(j) stats::var(G[, j]) <= 0, logical(1))
  if (any(zero_var)) {
    stop("`G` contains zero-variance SNP columns.", call. = FALSE)
  }
  if (is.null(colnames(G))) {
    colnames(G) <- paste0("snp_", seq_len(ncol(G)))
  }

  time <- as.numeric(time)
  status <- as.numeric(status)
  if (length(time) != nrow(G) || length(status) != nrow(G)) {
    stop("`time`, `status`, and `G` rows must have the same length.", call. = FALSE)
  }
  if (any(!is.finite(time)) || any(time < 0)) {
    stop("`time` must contain finite non-negative values.", call. = FALSE)
  }
  if (any(!is.finite(status)) || any(!status %in% c(0, 1))) {
    stop("`status` must be binary with values 0 or 1.", call. = FALSE)
  }
  if (!any(status == 1)) {
    stop("SDPD requires at least one observed event.", call. = FALSE)
  }

  out <- list(G = G, time = time, status = as.integer(status))
  if (!is.null(X)) {
    X <- as.numeric(X)
    if (length(X) != nrow(G)) {
      stop("`X` must have one value per row of `G`.", call. = FALSE)
    }
    if (any(!is.finite(X))) {
      stop("`X` must contain finite numeric values.", call. = FALSE)
    }
    out$X <- X
  }
  out
}

.mrwin_validate_sdpd_alpha <- function(alpha) {
  alpha <- as.numeric(alpha)
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be between 0 and 1.", call. = FALSE)
  }
  alpha
}

.mrwin_tail_cumsum <- function(x) {
  x <- as.matrix(x)
  out <- apply(x, 2L, function(col) rev(cumsum(rev(col))))
  matrix(out, nrow = nrow(x), ncol = ncol(x))
}

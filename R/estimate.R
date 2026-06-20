mrwin_estimate <- function(
    time,
    status,
    G,
    X,
    beta_hat,
    n_strata = 10L,
    kernel = NULL,
    weights = NULL,
    block_size = 4000L,
    floor = 1e-12,
    active_strata = NULL,
    stratification = c("prs_rank", "doubly_ranked")
) {
  stratification <- match.arg(stratification)
  checked <- .mrwin_validate_estimate_inputs(
    time = time,
    status = status,
    G = G,
    X = X,
    beta_hat = beta_hat,
    n_strata = n_strata,
    kernel = kernel,
    weights = weights,
    floor = floor,
    active_strata = active_strata
  )
  time <- checked$time
  status <- checked$status
  G <- checked$G
  X <- checked$X
  beta_hat <- checked$beta_hat
  n_strata <- checked$n_strata
  weights <- checked$weights
  floor <- checked$floor
  active_strata <- checked$active_strata

  if (is.null(kernel)) {
    kernel <- mrwin_kernel(time, status, block_size = block_size)
  } else {
    kernel <- checked$kernel
  }

  strata_obj <- .mrwin_assign_strata(G, beta_hat, X, n_strata, stratification)
  strata <- strata_obj$strata
  x <- as.numeric(X)

  contrast_plan <- .mrwin_make_contrast_plan(active_strata, n_strata = n_strata)
  n_contrasts <- nrow(contrast_plan)
  log_theta <- cwr <- delta_x <- rep(NA_real_, n_contrasts)
  wins <- losses <- total <- rep(NA_real_, n_contrasts)
  contrast_names <- contrast_plan$label

  for (row in seq_len(n_contrasts)) {
    idx_high <- which(strata == contrast_plan$high[row])
    idx_low <- which(strata == contrast_plan$low[row])
    sums <- mrwin_stratum_win_loss(kernel, idx_high, idx_low, weights = weights)
    wins[row] <- sums[["wins"]]
    losses[row] <- sums[["losses"]]
    total[row] <- sums[["total"]]
    log_theta[row] <- log(max(wins[row], floor) / max(losses[row], floor))
    cwr[row] <- exp(log_theta[row])

    if (is.null(weights)) {
      delta_x[row] <- mean(x[idx_high]) - mean(x[idx_low])
    } else {
      wh <- weights[idx_high]
      wl <- weights[idx_low]
      delta_x[row] <- .mrwin_weighted_mean(x[idx_high], wh, "high stratum") -
        .mrwin_weighted_mean(x[idx_low], wl, "low stratum")
    }
  }

  delta_isg <- log_theta / delta_x
  names(delta_isg) <- names(log_theta) <- names(cwr) <- names(delta_x) <-
    names(wins) <- names(losses) <- names(total) <- contrast_names

  structure(list(
    strata = strata,
    score = strata_obj$score,
    kernel = kernel,
    wins = wins,
    losses = losses,
    total = total,
    cwr = cwr,
    log_theta = log_theta,
    delta_x = delta_x,
    delta_isg = delta_isg,
    active_strata = active_strata,
    contrast_plan = contrast_plan,
    weak_delta_x = !is.finite(delta_x) | abs(delta_x) <= floor
  ), class = c("mrwin_estimate", "list"))
}

mrwin_gls_pool <- function(delta_isg, sigma, shrink = TRUE) {
  delta_names <- names(delta_isg)
  delta_isg <- as.numeric(delta_isg)
  if (is.null(delta_names)) {
    delta_names <- paste0("contrast_", seq_along(delta_isg))
  }
  sigma <- as.matrix(sigma)
  if (length(delta_isg) == 0L) {
    stop("`delta_isg` must contain at least one contrast.", call. = FALSE)
  }
  if (any(!is.finite(delta_isg))) {
    stop("`delta_isg` must contain finite values.", call. = FALSE)
  }
  if (nrow(sigma) != length(delta_isg) || ncol(sigma) != length(delta_isg)) {
    stop("`sigma` must be a square matrix with length(delta_isg) rows.", call. = FALSE)
  }
  if (any(!is.finite(sigma))) {
    stop("`sigma` must contain finite values.", call. = FALSE)
  }
  if (max(abs(sigma - t(sigma))) > sqrt(.Machine$double.eps)) {
    stop("`sigma` must be symmetric.", call. = FALSE)
  }
  sigma_eigenvalues <- eigen(0.5 * (sigma + t(sigma)), symmetric = TRUE, only.values = TRUE)$values
  if (min(sigma_eigenvalues) < -sqrt(.Machine$double.eps)) {
    stop("`sigma` must be positive semi-definite.", call. = FALSE)
  }

  sigma_used <- if (shrink) .mrwin_ledoit_wolf(sigma) else list(sigma = sigma, rho = 0)
  inv_sigma <- .mrwin_pseudo_inverse(sigma_used$sigma)
  ones <- rep(1, length(delta_isg))
  precision <- drop(t(ones) %*% inv_sigma %*% ones)
  if (!is.finite(precision) || precision <= 0) {
    stop("`sigma` must imply positive finite GLS precision.", call. = FALSE)
  }
  var_gls <- 1 / precision
  delta_gls <- var_gls * drop(t(ones) %*% inv_sigma %*% delta_isg)
  se_gls <- sqrt(max(var_gls, 0))
  q <- max(drop(t(delta_isg - delta_gls) %*% inv_sigma %*% (delta_isg - delta_gls)), 0)
  df <- length(delta_isg) - 1L
  q_p_value <- if (df > 0L) stats::pchisq(q, df = df, lower.tail = FALSE) else NA_real_
  ci95_delta <- delta_gls + c(-1.96, 1.96) * se_gls
  names(ci95_delta) <- c("low", "high")
  ci95_dscwr <- exp(ci95_delta)
  names(ci95_dscwr) <- c("low", "high")
  gls_weights <- as.numeric(var_gls * drop(t(ones) %*% inv_sigma))
  names(gls_weights) <- delta_names

  list(
    delta_gls = delta_gls,
    se_delta_gls = se_gls,
    dscwr = exp(delta_gls),
    ci95_delta = ci95_delta,
    ci95_dscwr = ci95_dscwr,
    q = q,
    q_df = df,
    q_p_value = q_p_value,
    gls_weights = gls_weights,
    sigma = sigma_used$sigma,
    ledoit_wolf_rho = sigma_used$rho
  )
}

.mrwin_validate_estimate_inputs <- function(
    time,
    status,
    G,
    X,
    beta_hat,
    n_strata,
    kernel,
    weights,
    floor,
    active_strata = NULL
) {
  time <- as.matrix(time)
  status <- as.matrix(status)
  G <- as.matrix(G)
  X <- as.numeric(X)
  beta_hat <- as.numeric(beta_hat)
  n_strata_raw <- suppressWarnings(as.numeric(n_strata))
  floor <- as.numeric(floor)

  if (!all(dim(time) == dim(status))) {
    stop("`time` and `status` must have the same dimensions.", call. = FALSE)
  }
  if (nrow(time) != nrow(G) || length(X) != nrow(G)) {
    stop("`time`, `status`, `G`, and `X` must describe the same number of rows.", call. = FALSE)
  }
  if (ncol(G) != length(beta_hat)) {
    stop("`G` columns must match length of `beta_hat`.", call. = FALSE)
  }
  if (length(n_strata_raw) != 1L || is.na(n_strata_raw) || n_strata_raw < 2L ||
      n_strata_raw != base::floor(n_strata_raw)) {
    stop("`n_strata` must be a single integer at least 2.", call. = FALSE)
  }
  n_strata <- as.integer(n_strata_raw)
  if (n_strata > nrow(G)) {
    stop("`n_strata` cannot exceed the number of rows.", call. = FALSE)
  }
  if (length(floor) != 1L || !is.finite(floor) || floor <= 0) {
    stop("`floor` must be a single positive finite value.", call. = FALSE)
  }
  if (any(!is.finite(time))) {
    stop("`time` must contain finite observed times.", call. = FALSE)
  }
  if (!all(status %in% c(0, 1))) {
    stop("`status` values must be binary 0/1.", call. = FALSE)
  }
  if (any(!is.finite(G)) || any(!is.finite(X)) || any(!is.finite(beta_hat))) {
    stop("`G`, `X`, and `beta_hat` must contain finite values.", call. = FALSE)
  }

  if (!is.null(kernel)) {
    kernel <- as.matrix(kernel)
    if (!all(dim(kernel) == c(nrow(G), nrow(G)))) {
      stop("`kernel` must be an N x N matrix matching the analysis rows.", call. = FALSE)
    }
    if (!all(kernel %in% c(-1, 0, 1))) {
      stop("`kernel` values must be -1, 0, or 1.", call. = FALSE)
    }
  }

  if (!is.null(weights)) {
    weights <- as.numeric(weights)
    if (length(weights) != nrow(G)) {
      stop("`weights` length must match the number of rows.", call. = FALSE)
    }
    if (any(!is.finite(weights)) || any(weights < 0)) {
      stop("`weights` must be finite and non-negative.", call. = FALSE)
    }
  }
  if (is.null(active_strata)) {
    active_strata <- seq_len(n_strata)
  } else {
    active_strata <- suppressWarnings(as.integer(active_strata))
    if (length(active_strata) < 2L || any(is.na(active_strata)) ||
        any(active_strata < 1L | active_strata > n_strata)) {
      stop("`active_strata` must contain at least two valid stratum labels.", call. = FALSE)
    }
    active_strata <- sort(unique(active_strata))
  }

  list(
    time = time,
    status = status,
    G = G,
    X = X,
    beta_hat = beta_hat,
    n_strata = n_strata,
    kernel = kernel,
    weights = weights,
    floor = floor,
    active_strata = active_strata
  )
}

.mrwin_contrast_names <- function(n_strata) {
  paste0("s", 2:n_strata, "_vs_s", 1:(n_strata - 1L))
}

.mrwin_weighted_mean <- function(x, weights, label) {
  total_weight <- sum(weights)
  if (!is.finite(total_weight) || total_weight <= 0) {
    stop("Weights for the ", label, " must have positive total weight.", call. = FALSE)
  }
  sum(weights * x) / total_weight
}

.mrwin_ledoit_wolf <- function(sigma) {
  sigma <- 0.5 * (sigma + t(sigma))
  p <- nrow(sigma)
  mu <- sum(diag(sigma)) / p
  target <- diag(mu, p)
  denom <- sum(sigma * sigma)
  rho <- if (denom <= 1e-12) 1 else min(1, max(0, sum((sigma - target)^2) / max(p, 1) / denom))
  list(sigma = (1 - rho) * sigma + rho * target, rho = rho)
}

.mrwin_pseudo_inverse <- function(x, tolerance = sqrt(.Machine$double.eps)) {
  x <- 0.5 * (x + t(x))
  eig <- eigen(x, symmetric = TRUE)
  cutoff <- max(abs(eig$values), 1) * tolerance
  inv_values <- ifelse(abs(eig$values) > cutoff, 1 / eig$values, 0)
  eig$vectors %*% (inv_values * t(eig$vectors))
}

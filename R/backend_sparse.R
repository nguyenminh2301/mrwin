mrwin_precompute_pair_kernels <- function(time, status, strata, contrast_plan) {
  time <- as.matrix(time)
  status <- as.matrix(status)
  strata <- as.integer(strata)
  n_contrasts <- nrow(contrast_plan)
  pair_kernels <- vector("list", n_contrasts)

  for (row in seq_len(n_contrasts)) {
    idx_high <- which(strata == contrast_plan$high[row])
    idx_low <- which(strata == contrast_plan$low[row])
    if (length(idx_high) == 0L || length(idx_low) == 0L) {
      pair_kernels[[row]] <- NULL
      next
    }
    pair_kernels[[row]] <- list(
      high = contrast_plan$high[row],
      low = contrast_plan$low[row],
      idx_high = idx_high,
      idx_low = idx_low,
      kernel = mrwin_pair_kernel(
        time[idx_high, , drop = FALSE],
        status[idx_high, , drop = FALSE],
        time[idx_low, , drop = FALSE],
        status[idx_low, , drop = FALSE]
      )
    )
  }
  pair_kernels
}

mrwin_sparse_estimate <- function(
    time,
    status,
    G,
    X,
    beta_hat,
    n_strata = 10L,
    weights = NULL,
    floor = 1e-12,
    active_strata = NULL,
    pair_kernels = NULL,
    fast = FALSE,
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
    kernel = NULL,
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

  strata_obj <- .mrwin_assign_strata(G, beta_hat, X, n_strata, stratification)
  strata <- strata_obj$strata
  x <- as.numeric(X)

  contrast_plan <- .mrwin_make_contrast_plan(active_strata, n_strata = n_strata)
  n_contrasts <- nrow(contrast_plan)
  log_theta <- cwr <- delta_x <- rep(NA_real_, n_contrasts)
  wins <- losses <- total <- rep(NA_real_, n_contrasts)
  contrast_names <- contrast_plan$label

  precomputed <- !is.null(pair_kernels)

  for (row in seq_len(n_contrasts)) {
    if (precomputed && !is.null(pair_kernels[[row]])) {
      pk <- pair_kernels[[row]]
      idx_high <- pk$idx_high
      idx_low <- pk$idx_low
      h_block <- pk$kernel
      if (is.null(weights)) {
        w_sum <- sum(h_block == 1L)
        l_sum <- sum(h_block == -1L)
        t_sum <- length(h_block)
      } else {
        pair_w <- outer(weights[idx_high], weights[idx_low], "*")
        w_sum <- sum(pair_w * (h_block == 1L))
        l_sum <- sum(pair_w * (h_block == -1L))
        t_sum <- sum(pair_w)
      }
      wins[row] <- w_sum
      losses[row] <- l_sum
      total[row] <- t_sum
    } else {
      idx_high <- which(strata == contrast_plan$high[row])
      idx_low <- which(strata == contrast_plan$low[row])
      sums <- .mrwin_pair_win_loss_backend(
        time[idx_high, , drop = FALSE],
        status[idx_high, , drop = FALSE],
        time[idx_low, , drop = FALSE],
        status[idx_low, , drop = FALSE],
        weights_high = if (is.null(weights)) NULL else weights[idx_high],
        weights_low = if (is.null(weights)) NULL else weights[idx_low],
        fast = fast
      )
      wins[row] <- sums[["wins"]]
      losses[row] <- sums[["losses"]]
      total[row] <- sums[["total"]]
      idx_high <- which(strata == contrast_plan$high[row])
      idx_low <- which(strata == contrast_plan$low[row])
    }

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
    kernel = NULL,
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

mrwin_sparse_bootstrap <- function(
    time,
    status,
    G,
    X,
    beta_hat,
    sigma_beta = 0,
    n_strata = 10L,
    B = 200L,
    seed = NULL,
    beta_draws = NULL,
    multiplier_weights = NULL,
    floor = 1e-12,
    covariates = NULL,
    adjustment = c("none", "ordinal_iptw"),
    iptw_truncation = c(0.01, 0.99),
    ess_fraction = 0.5,
    fast = FALSE,
    stratification = c("prs_rank", "doubly_ranked")
) {
  adjustment <- match.arg(adjustment)
  stratification <- match.arg(stratification)
  checked <- .mrwin_validate_bootstrap_inputs(
    time = time,
    status = status,
    G = G,
    X = X,
    beta_hat = beta_hat,
    sigma_beta = sigma_beta,
    n_strata = n_strata,
    B = B,
    kernel = NULL,
    block_size = 1L,
    beta_draws = beta_draws,
    multiplier_weights = multiplier_weights,
    floor = floor,
    covariates = covariates,
    adjustment = adjustment,
    iptw_truncation = iptw_truncation,
    ess_fraction = ess_fraction
  )
  time <- checked$time
  status <- checked$status
  G <- checked$G
  X <- checked$X
  beta_hat <- checked$beta_hat
  sigma_beta <- checked$sigma_beta
  n_strata <- checked$n_strata
  B <- checked$B
  beta_draws <- checked$beta_draws
  multiplier_weights <- checked$multiplier_weights
  floor <- checked$floor
  covariates <- checked$covariates
  adjustment <- checked$adjustment
  iptw_truncation <- checked$iptw_truncation
  ess_fraction <- checked$ess_fraction

  if (!is.null(seed)) {
    set.seed(seed)
  }

  point_strata_obj <- .mrwin_assign_strata(G, beta_hat, X, n_strata, stratification)
  point_strata <- point_strata_obj$strata
  point_adjustment <- .mrwin_point_adjustment(
    strata = point_strata,
    covariates = covariates,
    adjustment = adjustment,
    iptw_truncation = iptw_truncation,
    ess_fraction = ess_fraction,
    n_strata = n_strata
  )

  contrast_plan <- .mrwin_make_contrast_plan(point_adjustment$active_strata, n_strata = n_strata)

  point_pair_kernels <- if (isTRUE(fast)) {
    NULL
  } else {
    mrwin_precompute_pair_kernels(time, status, point_strata, contrast_plan)
  }

  point <- mrwin_sparse_estimate(
    time = time,
    status = status,
    G = G,
    X = X,
    beta_hat = beta_hat,
    n_strata = n_strata,
    weights = point_adjustment$weights,
    floor = floor,
    active_strata = point_adjustment$active_strata,
    pair_kernels = point_pair_kernels,
    fast = fast,
    stratification = stratification
  )

  n_contrasts <- nrow(contrast_plan)
  lt <- matrix(NA_real_, nrow = B, ncol = n_contrasts)
  dx <- matrix(NA_real_, nrow = B, ncol = n_contrasts)
  colnames(lt) <- colnames(dx) <- names(point$log_theta)
  bootstrap_ess <- matrix(NA_real_, nrow = B, ncol = n_strata)
  colnames(bootstrap_ess) <- paste0("s", seq_len(n_strata))
  bootstrap_dropped <- matrix(FALSE, nrow = B, ncol = n_strata)
  colnames(bootstrap_dropped) <- paste0("s", seq_len(n_strata))

  for (b in seq_len(B)) {
    beta_star <- if (is.null(beta_draws)) {
      beta_hat + stats::rnorm(length(beta_hat)) * sigma_beta
    } else {
      beta_draws[b, ]
    }
    strata <- .mrwin_assign_strata(G, beta_star, X, n_strata, stratification)$strata
    xi <- if (is.null(multiplier_weights)) {
      stats::rexp(nrow(G), rate = 1)
    } else {
      multiplier_weights[b, ]
    }
    iter_adjustment <- .mrwin_iteration_adjustment(
      strata = strata,
      covariates = covariates,
      adjustment = adjustment,
      xi = xi,
      iptw_truncation = iptw_truncation,
      ess_fraction = ess_fraction,
      n_strata = n_strata
    )
    if (!is.null(iter_adjustment)) {
      bootstrap_ess[b, ] <- iter_adjustment$ess
      bootstrap_dropped[b, iter_adjustment$dropped_strata] <- TRUE
      if (any(point_adjustment$active_strata %in% iter_adjustment$dropped_strata)) {
        next
      }
      xi <- xi * iter_adjustment$weights
    }

    for (row in seq_len(n_contrasts)) {
      idx_high <- which(strata == contrast_plan$high[row])
      idx_low <- which(strata == contrast_plan$low[row])
      if (length(idx_high) == 0L || length(idx_low) == 0L) {
        next
      }
      wh <- xi[idx_high]
      wl <- xi[idx_low]
      if (sum(wh) <= 0 || sum(wl) <= 0) {
        next
      }

      sums <- .mrwin_pair_win_loss_backend(
        time[idx_high, , drop = FALSE],
        status[idx_high, , drop = FALSE],
        time[idx_low, , drop = FALSE],
        status[idx_low, , drop = FALSE],
        weights_high = wh,
        weights_low = wl,
        fast = fast
      )
      lt[b, row] <- log(max(sums[["wins"]], floor) / max(sums[["losses"]], floor))

      dx[b, row] <- .mrwin_weighted_mean(X[idx_high], wh, "bootstrap high stratum") -
        .mrwin_weighted_mean(X[idx_low], wl, "bootstrap low stratum")
    }
  }

  u <- cbind(lt, dx)
  colnames(u) <- c(paste0("log_theta:", colnames(lt)), paste0("delta_x:", colnames(dx)))
  valid <- stats::complete.cases(u)
  if (sum(valid) < 2L) {
    stop("Too few valid bootstrap iterations.", call. = FALSE)
  }
  cov_u <- .mrwin_bootstrap_covariance(u, valid)
  sigma_isg <- .mrwin_isg_covariance(
    cov_u = cov_u,
    point_log_theta = point$log_theta,
    point_delta_x = point$delta_x,
    floor = floor
  )

  pooled <- mrwin_gls_pool(point$delta_isg, sigma_isg, shrink = TRUE)
  fieller <- .mrwin_fieller_ci(
    cov_u = cov_u,
    sigma_isg = sigma_isg,
    point_log_theta = point$log_theta,
    point_delta_x = point$delta_x,
    z = 1.96
  )

  structure(list(
    B = as.integer(B),
    n_valid = as.integer(sum(valid)),
    n_invalid = as.integer(sum(!valid)),
    n_strata = as.integer(n_strata),
    active_strata = point_adjustment$active_strata,
    dropped_strata = point_adjustment$dropped_strata,
    contrast_plan = contrast_plan,
    point_log_theta = point$log_theta,
    point_cwr = point$cwr,
    point_delta_x = point$delta_x,
    delta_isg = point$delta_isg,
    delta_gls = pooled$delta_gls,
    se_delta_gls = pooled$se_delta_gls,
    dscwr = pooled$dscwr,
    ci95_delta = pooled$ci95_delta,
    ci95_dscwr = pooled$ci95_dscwr,
    ci95_delta_bivariate_delta = pooled$ci95_delta,
    ci95_dscwr_bivariate_delta = pooled$ci95_dscwr,
    ci95_delta_fieller = fieller$delta,
    ci95_dscwr_fieller = exp(pmax(pmin(fieller$delta, 50), -50)),
    fieller_unbounded = fieller$unbounded,
    fieller_coefficients = fieller$coefficients,
    fieller_p_value = fieller$p_value,
    q = pooled$q,
    q_df = pooled$q_df,
    q_p_value = pooled$q_p_value,
    gls_weights = pooled$gls_weights,
    ledoit_wolf_rho = pooled$ledoit_wolf_rho,
    cov_u = cov_u,
    sigma_isg = sigma_isg,
    sigma_lw = pooled$sigma,
    kurtosis_log_theta = apply(lt, 2, .mrwin_kurtosis),
    kurtosis_delta_x = apply(dx, 2, .mrwin_kurtosis),
    skew_log_theta = apply(lt, 2, .mrwin_skewness),
    skew_delta_x = apply(dx, 2, .mrwin_skewness),
    bootstrap_log_theta = lt,
    bootstrap_delta_x = dx,
    bootstrap_ess = bootstrap_ess,
    bootstrap_dropped_strata = bootstrap_dropped,
    valid_bootstrap = valid,
    covariance_method = "bivariate_delta_point_gradient",
    adjustment = point_adjustment,
    kernel = NULL,
    strata = point$strata
  ), class = c("mrwin_bootstrap", "list"))
}

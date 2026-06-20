mrwin_multiplier_bootstrap <- function(
    time,
    status,
    G,
    X,
    beta_hat,
    sigma_beta = 0,
    n_strata = 10L,
    B = 200L,
    seed = NULL,
    kernel = NULL,
    block_size = 4000L,
    beta_draws = NULL,
    multiplier_weights = NULL,
    floor = 1e-12,
    covariates = NULL,
    adjustment = c("none", "ordinal_iptw"),
    iptw_truncation = c(0.01, 0.99),
    ess_fraction = 0.5,
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
    kernel = kernel,
    block_size = block_size,
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
  kernel <- checked$kernel
  block_size <- checked$block_size
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

  if (is.null(kernel)) {
    kernel <- mrwin_kernel(time, status, block_size = block_size)
  }

  point_strata <- .mrwin_assign_strata(G, beta_hat, X, n_strata, stratification)$strata
  point_adjustment <- .mrwin_point_adjustment(
    strata = point_strata,
    covariates = covariates,
    adjustment = adjustment,
    iptw_truncation = iptw_truncation,
    ess_fraction = ess_fraction,
    n_strata = n_strata
  )

  point <- mrwin_estimate(
    time = time,
    status = status,
    G = G,
    X = X,
    beta_hat = beta_hat,
    n_strata = n_strata,
    kernel = kernel,
    weights = point_adjustment$weights,
    floor = floor,
    active_strata = point_adjustment$active_strata,
    stratification = stratification
  )

  contrast_plan <- point$contrast_plan
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
      sums <- mrwin_stratum_win_loss(kernel, idx_high, idx_low, weights = xi)
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
    kernel = kernel,
    strata = point$strata
  ), class = c("mrwin_bootstrap", "list"))
}

.mrwin_validate_bootstrap_inputs <- function(
    time,
    status,
    G,
    X,
    beta_hat,
    sigma_beta,
    n_strata,
    B,
    kernel,
    block_size,
    beta_draws,
    multiplier_weights,
    floor,
    covariates,
    adjustment,
    iptw_truncation,
    ess_fraction
) {
  checked <- .mrwin_validate_estimate_inputs(
    time = time,
    status = status,
    G = G,
    X = X,
    beta_hat = beta_hat,
    n_strata = n_strata,
    kernel = kernel,
    weights = NULL,
    floor = floor
  )
  if (adjustment != "none") {
    if (is.null(covariates)) {
      stop("`covariates` are required when adjustment is not 'none'.", call. = FALSE)
    }
    covariates <- as.matrix(covariates)
    if (nrow(covariates) != nrow(checked$G)) {
      stop("`covariates` must have one row per analysis row.", call. = FALSE)
    }
    if (any(!is.finite(covariates))) {
      stop("`covariates` must contain finite numeric values.", call. = FALSE)
    }
    iptw_truncation <- .mrwin_validate_truncation(iptw_truncation)
    ess_fraction <- .mrwin_validate_ess_fraction(ess_fraction)
  } else {
    covariates <- NULL
    iptw_truncation <- .mrwin_validate_truncation(iptw_truncation)
    ess_fraction <- .mrwin_validate_ess_fraction(ess_fraction)
  }
  B_raw <- suppressWarnings(as.numeric(B))
  if (length(B_raw) != 1L || is.na(B_raw) || B_raw < 2L || B_raw != base::floor(B_raw)) {
    stop("`B` must be a single integer at least 2.", call. = FALSE)
  }
  B <- as.integer(B_raw)

  block_size_raw <- suppressWarnings(as.numeric(block_size))
  if (length(block_size_raw) != 1L || is.na(block_size_raw) || block_size_raw < 1L ||
      block_size_raw != base::floor(block_size_raw)) {
    stop("`block_size` must be a positive integer.", call. = FALSE)
  }
  block_size <- as.integer(block_size_raw)

  sigma_beta <- as.numeric(sigma_beta)
  if (length(sigma_beta) == 1L) {
    sigma_beta <- rep(sigma_beta, length(checked$beta_hat))
  }
  if (length(sigma_beta) != length(checked$beta_hat)) {
    stop("`sigma_beta` must have length 1 or length(beta_hat).", call. = FALSE)
  }
  if (any(!is.finite(sigma_beta)) || any(sigma_beta < 0)) {
    stop("`sigma_beta` must be finite and non-negative.", call. = FALSE)
  }

  if (!is.null(beta_draws)) {
    beta_draws <- as.matrix(beta_draws)
    if (!all(dim(beta_draws) == c(B, length(checked$beta_hat)))) {
      stop("`beta_draws` must be a B x M matrix matching bootstrap iterations and beta length.", call. = FALSE)
    }
    if (any(!is.finite(beta_draws))) {
      stop("`beta_draws` must contain finite values.", call. = FALSE)
    }
  }

  if (!is.null(multiplier_weights)) {
    multiplier_weights <- as.matrix(multiplier_weights)
    if (!all(dim(multiplier_weights) == c(B, nrow(checked$G)))) {
      stop("`multiplier_weights` must be a B x N matrix matching bootstrap iterations and rows.", call. = FALSE)
    }
    if (any(!is.finite(multiplier_weights)) || any(multiplier_weights < 0)) {
      stop("`multiplier_weights` must be finite and non-negative.", call. = FALSE)
    }
  }

  c(
    checked,
    list(
      sigma_beta = sigma_beta,
      B = B,
      block_size = block_size,
      beta_draws = beta_draws,
      multiplier_weights = multiplier_weights,
      covariates = covariates,
      adjustment = adjustment,
      iptw_truncation = iptw_truncation,
      ess_fraction = ess_fraction
    )
  )
}

.mrwin_point_adjustment <- function(strata, covariates, adjustment, iptw_truncation, ess_fraction, n_strata = max(strata)) {
  if (adjustment == "none") {
    return(list(
      method = "none",
      weights = NULL,
      raw_weights = NULL,
      ess = rep(NA_real_, n_strata),
      ess_threshold = rep(NA_real_, n_strata),
      dropped_strata = integer(),
      active_strata = seq_len(n_strata),
      truncation = NULL,
      balance = NULL,
      has_positivity_failure = FALSE,
      has_bridging = FALSE
    ))
  }
  prop <- mrwin_propensity_weights(
    strata = strata,
    covariates = covariates,
    method = "ordinal_iptw",
    truncate = iptw_truncation,
    ess_fraction = ess_fraction,
    n_strata = n_strata
  )
  active <- prop$active_strata
  if (length(active) < 2L) {
    stop("Too few ESS-valid strata remain after IPTW positivity filtering.", call. = FALSE)
  }
  plan <- .mrwin_make_contrast_plan(active, n_strata = n_strata)
  prop$has_positivity_failure <- length(prop$dropped_strata) > 0L
  prop$has_bridging <- any(plan$bridged)
  prop
}

.mrwin_iteration_adjustment <- function(strata, covariates, adjustment, xi, iptw_truncation, ess_fraction, n_strata) {
  if (adjustment == "none") {
    return(NULL)
  }
  mrwin_propensity_weights(
    strata = strata,
    covariates = covariates,
    method = "ordinal_iptw",
    base_weights = xi,
    truncate = iptw_truncation,
    ess_fraction = ess_fraction,
    fail_on_empty = FALSE,
    n_strata = n_strata
  )
}

.mrwin_bootstrap_covariance <- function(u, valid) {
  cov_u <- stats::cov(u[valid, , drop = FALSE])
  if (any(!is.finite(cov_u))) {
    stop("Bootstrap covariance contains non-finite values.", call. = FALSE)
  }
  cov_u <- 0.5 * (cov_u + t(cov_u))
  cov_u
}

.mrwin_isg_covariance <- function(cov_u, point_log_theta, point_delta_x, floor = 1e-12) {
  dm1 <- length(point_log_theta)
  cov_u <- as.matrix(cov_u)
  if (!all(dim(cov_u) == c(2L * dm1, 2L * dm1))) {
    stop("`cov_u` must be a 2(D-1) x 2(D-1) covariance matrix.", call. = FALSE)
  }
  if (any(!is.finite(cov_u))) {
    stop("`cov_u` must contain finite values.", call. = FALSE)
  }
  if (any(!is.finite(point_log_theta)) || any(!is.finite(point_delta_x))) {
    stop("Point log-CWR and Delta-X values must be finite.", call. = FALSE)
  }
  if (any(abs(point_delta_x) <= floor)) {
    stop("Point Delta-X is too close to zero for bivariate Delta covariance.", call. = FALSE)
  }

  sigma <- matrix(0, nrow = dm1, ncol = dm1)
  colnames(sigma) <- rownames(sigma) <- names(point_log_theta)
  for (a in seq_len(dm1)) {
    for (c in seq_len(dm1)) {
      ga <- c(1 / point_delta_x[a], -point_log_theta[a] / point_delta_x[a]^2)
      gc <- c(1 / point_delta_x[c], -point_log_theta[c] / point_delta_x[c]^2)
      block <- matrix(
        c(
          cov_u[a, c],
          cov_u[a, c + dm1],
          cov_u[a + dm1, c],
          cov_u[a + dm1, c + dm1]
        ),
        nrow = 2,
        byrow = TRUE
      )
      sigma[a, c] <- drop(t(ga) %*% block %*% gc)
    }
  }
  sigma <- 0.5 * (sigma + t(sigma))
  if (any(!is.finite(sigma))) {
    stop("Bivariate Delta covariance contains non-finite values.", call. = FALSE)
  }
  sigma
}

.mrwin_fieller_ci <- function(cov_u, sigma_isg, point_log_theta, point_delta_x, z = 1.96) {
  dm1 <- length(point_log_theta)
  weights <- 1 / pmax(diag(sigma_isg), 1e-12)
  weights <- weights / sum(weights)

  u1 <- sum(weights * point_log_theta)
  u2 <- sum(weights * point_delta_x)
  cov_lt <- cov_u[seq_len(dm1), seq_len(dm1), drop = FALSE]
  cov_dx <- cov_u[dm1 + seq_len(dm1), dm1 + seq_len(dm1), drop = FALSE]
  cov_cross <- cov_u[seq_len(dm1), dm1 + seq_len(dm1), drop = FALSE]

  var_u1 <- drop(t(weights) %*% cov_lt %*% weights)
  var_u2 <- drop(t(weights) %*% cov_dx %*% weights)
  cov_u12 <- drop(t(weights) %*% cov_cross %*% weights)

  a <- u2^2 - z^2 * var_u2
  b <- -2 * (u1 * u2 - z^2 * cov_u12)
  cc <- u1^2 - z^2 * var_u1
  disc <- b^2 - 4 * a * cc

  if (is.finite(a) && a > 0 && is.finite(disc) && disc > 0) {
    roots <- sort(c((-b - sqrt(disc)) / (2 * a), (-b + sqrt(disc)) / (2 * a)))
  } else {
    roots <- c(-Inf, Inf)
  }
  names(roots) <- c("low", "high")
  # Fieller-consistent test of H0: delta = 0  <=>  pooled numerator u1 = 0.
  # Substituting delta = 0 into the Fieller inequality gives |u1|/sqrt(Var(u1)).
  z_null <- if (is.finite(var_u1) && var_u1 > 0) u1 / sqrt(var_u1) else NA_real_
  p_null <- if (is.finite(z_null)) 2 * stats::pnorm(-abs(z_null)) else NA_real_
  list(
    delta = roots,
    unbounded = any(!is.finite(roots)),
    coefficients = c(a = a, b = b, c = cc, discriminant = disc),
    z_null = z_null,
    p_value = p_null
  )
}

.mrwin_kurtosis <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) {
    return(NA_real_)
  }
  centered <- x - mean(x)
  v <- mean(centered^2)
  if (v <= 0) {
    return(NA_real_)
  }
  mean(centered^4) / v^2
}

.mrwin_skewness <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) {
    return(NA_real_)
  }
  centered <- x - mean(x)
  s <- sqrt(mean(centered^2))
  if (s <= 0) {
    return(NA_real_)
  }
  mean(centered^3) / s^3
}

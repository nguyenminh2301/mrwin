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
    floor = 1e-12
) {
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
    floor = floor
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

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (is.null(kernel)) {
    kernel <- mrwin_kernel(time, status, block_size = block_size)
  }

  point <- mrwin_estimate(
    time = time,
    status = status,
    G = G,
    X = X,
    beta_hat = beta_hat,
    n_strata = n_strata,
    kernel = kernel,
    floor = floor
  )

  dm1 <- n_strata - 1L
  lt <- matrix(NA_real_, nrow = B, ncol = dm1)
  dx <- matrix(NA_real_, nrow = B, ncol = dm1)
  colnames(lt) <- colnames(dx) <- names(point$log_theta)

  for (b in seq_len(B)) {
    beta_star <- if (is.null(beta_draws)) {
      beta_hat + stats::rnorm(length(beta_hat)) * sigma_beta
    } else {
      beta_draws[b, ]
    }
    strata <- mrwin_prs_strata(G, beta_star, n_strata = n_strata)$strata
    xi <- if (is.null(multiplier_weights)) {
      stats::rexp(nrow(G), rate = 1)
    } else {
      multiplier_weights[b, ]
    }

    for (d in 2:n_strata) {
      idx_high <- which(strata == d)
      idx_low <- which(strata == d - 1L)
      if (length(idx_high) == 0L || length(idx_low) == 0L) {
        next
      }
      wh <- xi[idx_high]
      wl <- xi[idx_low]
      if (sum(wh) <= 0 || sum(wl) <= 0) {
        next
      }
      sums <- mrwin_stratum_win_loss(kernel, idx_high, idx_low, weights = xi)
      lt[b, d - 1L] <- log(max(sums[["wins"]], floor) / max(sums[["losses"]], floor))

      dx[b, d - 1L] <- .mrwin_weighted_mean(X[idx_high], wh, "bootstrap high stratum") -
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
    valid_bootstrap = valid,
    covariance_method = "bivariate_delta_point_gradient",
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
    floor
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
      multiplier_weights = multiplier_weights
    )
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
  list(
    delta = roots,
    unbounded = any(!is.finite(roots)),
    coefficients = c(a = a, b = b, c = cc, discriminant = disc)
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

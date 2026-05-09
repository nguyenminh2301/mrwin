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
    block_size = 4000L
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  G <- as.matrix(G)
  beta_hat <- as.numeric(beta_hat)
  sigma_beta <- rep(as.numeric(sigma_beta), length.out = length(beta_hat))
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
    kernel = kernel
  )

  dm1 <- n_strata - 1L
  lt <- matrix(NA_real_, nrow = B, ncol = dm1)
  dx <- matrix(NA_real_, nrow = B, ncol = dm1)

  for (b in seq_len(B)) {
    beta_star <- beta_hat + stats::rnorm(length(beta_hat)) * sigma_beta
    strata <- mrwin_prs_strata(G, beta_star, n_strata = n_strata)$strata
    xi <- stats::rexp(nrow(G), rate = 1)

    for (d in 2:n_strata) {
      idx_high <- which(strata == d)
      idx_low <- which(strata == d - 1L)
      if (length(idx_high) == 0L || length(idx_low) == 0L) {
        next
      }
      sums <- mrwin_stratum_win_loss(kernel, idx_high, idx_low, weights = xi)
      lt[b, d - 1L] <- log(max(sums[["wins"]], 1e-12) / max(sums[["losses"]], 1e-12))

      wh <- xi[idx_high]
      wl <- xi[idx_low]
      dx[b, d - 1L] <- sum(wh * X[idx_high]) / sum(wh) -
        sum(wl * X[idx_low]) / sum(wl)
    }
  }

  u <- cbind(lt, dx)
  valid <- stats::complete.cases(u)
  if (sum(valid) < 2L) {
    stop("Too few valid bootstrap iterations.", call. = FALSE)
  }
  cov_u <- stats::cov(u[valid, , drop = FALSE])
  sigma_isg <- .mrwin_isg_covariance(
    cov_u = cov_u,
    point_log_theta = point$log_theta,
    point_delta_x = point$delta_x
  )

  pooled <- mrwin_gls_pool(point$delta_isg, sigma_isg, shrink = TRUE)
  fieller <- .mrwin_fieller_ci(
    cov_u = cov_u,
    sigma_isg = sigma_isg,
    point_log_theta = point$log_theta,
    point_delta_x = point$delta_x
  )

  list(
    B = as.integer(B),
    n_valid = as.integer(sum(valid)),
    n_strata = as.integer(n_strata),
    point_log_theta = point$log_theta,
    point_delta_x = point$delta_x,
    delta_isg = point$delta_isg,
    delta_gls = pooled$delta_gls,
    se_delta_gls = pooled$se_delta_gls,
    dscwr = pooled$dscwr,
    ci95_delta = pooled$ci95_delta,
    ci95_dscwr = pooled$ci95_dscwr,
    ci95_delta_fieller = fieller$delta,
    ci95_dscwr_fieller = exp(pmax(pmin(fieller$delta, 50), -50)),
    q = pooled$q,
    q_df = pooled$q_df,
    q_p_value = pooled$q_p_value,
    ledoit_wolf_rho = pooled$ledoit_wolf_rho,
    sigma_isg = sigma_isg,
    sigma_lw = pooled$sigma,
    kurtosis_log_theta = apply(lt, 2, .mrwin_kurtosis),
    kurtosis_delta_x = apply(dx, 2, .mrwin_kurtosis),
    skew_log_theta = apply(lt, 2, .mrwin_skewness),
    skew_delta_x = apply(dx, 2, .mrwin_skewness),
    bootstrap_log_theta = lt,
    bootstrap_delta_x = dx,
    kernel = kernel,
    strata = point$strata
  )
}

.mrwin_isg_covariance <- function(cov_u, point_log_theta, point_delta_x) {
  dm1 <- length(point_log_theta)
  sigma <- matrix(0, nrow = dm1, ncol = dm1)
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
  0.5 * (sigma + t(sigma))
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
  list(delta = roots)
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

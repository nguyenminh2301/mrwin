mrwin_estimate <- function(
    time,
    status,
    G,
    X,
    beta_hat,
    n_strata = 10L,
    kernel = NULL,
    weights = NULL,
    block_size = 4000L
) {
  if (is.null(kernel)) {
    kernel <- mrwin_kernel(time, status, block_size = block_size)
  }

  strata_obj <- mrwin_prs_strata(G, beta_hat, n_strata = n_strata)
  strata <- strata_obj$strata
  x <- as.numeric(X)

  log_theta <- delta_x <- rep(NA_real_, n_strata - 1L)
  wins <- losses <- total <- rep(NA_real_, n_strata - 1L)

  for (d in 2:n_strata) {
    idx_high <- which(strata == d)
    idx_low <- which(strata == d - 1L)
    sums <- mrwin_stratum_win_loss(kernel, idx_high, idx_low, weights = weights)
    wins[d - 1L] <- sums[["wins"]]
    losses[d - 1L] <- sums[["losses"]]
    total[d - 1L] <- sums[["total"]]
    log_theta[d - 1L] <- log(max(wins[d - 1L], 1e-12) / max(losses[d - 1L], 1e-12))

    if (is.null(weights)) {
      delta_x[d - 1L] <- mean(x[idx_high]) - mean(x[idx_low])
    } else {
      wh <- weights[idx_high]
      wl <- weights[idx_low]
      delta_x[d - 1L] <- sum(wh * x[idx_high]) / sum(wh) -
        sum(wl * x[idx_low]) / sum(wl)
    }
  }

  delta_isg <- log_theta / delta_x
  names(delta_isg) <- names(log_theta) <- names(delta_x) <-
    paste0("s", 2:n_strata, "_vs_s", 1:(n_strata - 1L))

  list(
    strata = strata,
    score = strata_obj$score,
    kernel = kernel,
    wins = wins,
    losses = losses,
    total = total,
    log_theta = log_theta,
    delta_x = delta_x,
    delta_isg = delta_isg
  )
}

mrwin_gls_pool <- function(delta_isg, sigma, shrink = TRUE) {
  delta_isg <- as.numeric(delta_isg)
  sigma <- as.matrix(sigma)
  if (nrow(sigma) != length(delta_isg) || ncol(sigma) != length(delta_isg)) {
    stop("`sigma` must be a square matrix with length(delta_isg) rows.", call. = FALSE)
  }

  sigma_used <- if (shrink) .mrwin_ledoit_wolf(sigma) else list(sigma = sigma, rho = 0)
  inv_sigma <- .mrwin_pseudo_inverse(sigma_used$sigma)
  ones <- rep(1, length(delta_isg))
  var_gls <- 1 / drop(t(ones) %*% inv_sigma %*% ones)
  delta_gls <- var_gls * drop(t(ones) %*% inv_sigma %*% delta_isg)
  se_gls <- sqrt(max(var_gls, 0))
  q <- drop(t(delta_isg - delta_gls) %*% inv_sigma %*% (delta_isg - delta_gls))
  df <- max(length(delta_isg) - 1L, 1L)

  list(
    delta_gls = delta_gls,
    se_delta_gls = se_gls,
    dscwr = exp(delta_gls),
    ci95_delta = delta_gls + c(-1.96, 1.96) * se_gls,
    ci95_dscwr = exp(delta_gls + c(-1.96, 1.96) * se_gls),
    q = q,
    q_df = df,
    q_p_value = stats::pchisq(q, df = df, lower.tail = FALSE),
    sigma = sigma_used$sigma,
    ledoit_wolf_rho = sigma_used$rho
  )
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

mrwin_simulate <- function(config = mrwin_config(), seed = config$seed) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  n <- config$n_outcome
  m <- config$m_snps
  k <- 3L

  mafs <- stats::runif(m, config$maf_low, config$maf_high)
  g <- .mrwin_simulate_genotypes(n, mafs)

  true_betas <- stats::rnorm(m, 0, config$sigma_beta)
  s_true <- drop(g %*% true_betas)
  u <- stats::rnorm(n)

  if (config$theta_f > 0) {
    w <- stats::rgamma(n, shape = 1 / config$theta_f, scale = config$theta_f)
  } else {
    w <- rep(1, n)
  }

  x <- config$alpha_s * s_true + config$alpha_u * u + stats::rnorm(n)
  log_w <- log(w)

  event_time <- matrix(NA_real_, nrow = n, ncol = k)
  for (j in seq_len(k)) {
    lp <- config$alpha_x[j] * x +
      config$nu_u[j] * u +
      config$gamma_direct[j] * s_true +
      log_w
    event_time[, j] <- .mrwin_weibull_inv(
      rate = config$baseline_haz[j],
      shape = config$shape_weibull,
      lp = lp
    )
  }

  censor_time <- pmin(
    stats::rexp(n, rate = config$censoring_rate),
    config$max_follow_up
  )

  t_death <- event_time[, 1L]
  t_hf <- event_time[, 2L]
  t_renal <- event_time[, 3L]

  time <- cbind(
    pmin(t_death, censor_time),
    pmin(t_hf, censor_time, t_death),
    pmin(t_renal, censor_time, t_death)
  )
  status <- cbind(
    as.integer(t_death <= censor_time),
    as.integer(t_hf <= censor_time & t_hf <= t_death),
    as.integer(t_renal <= censor_time & t_renal <= t_death)
  )
  colnames(time) <- colnames(status) <- c("death", "hf", "renal")

  list(
    G = g,
    mafs = mafs,
    true_betas = true_betas,
    S_true = s_true,
    U = u,
    W = w,
    X = x,
    time = time,
    status = status,
    event_time = event_time,
    censor_time = censor_time,
    config = config
  )
}

.mrwin_simulate_genotypes <- function(n, mafs, max_attempts = 100L) {
  m <- length(mafs)
  g_raw <- matrix(NA_real_, nrow = n, ncol = m)
  for (j in seq_len(m)) {
    for (attempt in seq_len(max_attempts)) {
      candidate <- stats::rbinom(n, 2L, mafs[j])
      if (stats::var(candidate) > 0) {
        g_raw[, j] <- candidate
        break
      }
    }
    if (anyNA(g_raw[, j])) {
      g_raw[, j] <- rep(c(0, 1, 2), length.out = n)
    }
  }
  g <- scale(g_raw)
  storage.mode(g) <- "double"
  colnames(g) <- paste0("snp_", seq_len(m))
  g
}

.mrwin_weibull_inv <- function(rate, shape, lp) {
  (-log(stats::runif(length(lp))) / (rate * exp(lp)))^(1 / shape)
}

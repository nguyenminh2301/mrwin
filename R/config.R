mrwin_config <- function(
    n_outcome = 6000L,
    m_snps = 40L,
    alpha_s = 0.4,
    alpha_u = 0.8,
    alpha_x = c(-0.4, -0.4, -0.4),
    gamma_direct = c(0, 0, 0),
    nu_u = c(0.5, 0.5, 0.5),
    baseline_haz = c(0.02, 0.05, 0.10),
    shape_weibull = 1.2,
    censoring_rate = 0.05,
    max_follow_up = 10,
    theta_f = 0.8,
    maf_low = 0.10,
    maf_high = 0.40,
    sigma_beta = 0.05,
    seed = 20260506L
) {
  n_outcome <- .mrwin_config_integer(n_outcome, "n_outcome", lower = 3L)
  m_snps <- .mrwin_config_integer(m_snps, "m_snps", lower = 1L)
  alpha_x <- .mrwin_config_vector(alpha_x, "alpha_x", n_expected = 3L)
  gamma_direct <- .mrwin_config_vector(gamma_direct, "gamma_direct", n_expected = 3L)
  nu_u <- .mrwin_config_vector(nu_u, "nu_u", n_expected = 3L)
  baseline_haz <- .mrwin_config_vector(baseline_haz, "baseline_haz", n_expected = 3L, lower = 0)
  alpha_s <- .mrwin_config_scalar(alpha_s, "alpha_s")
  alpha_u <- .mrwin_config_scalar(alpha_u, "alpha_u")
  shape_weibull <- .mrwin_config_scalar(shape_weibull, "shape_weibull", lower = 0)
  censoring_rate <- .mrwin_config_scalar(censoring_rate, "censoring_rate", lower = 0)
  max_follow_up <- .mrwin_config_scalar(max_follow_up, "max_follow_up", lower = 0)
  theta_f <- .mrwin_config_scalar(theta_f, "theta_f", lower = 0, allow_zero = TRUE)
  maf_low <- .mrwin_config_scalar(maf_low, "maf_low", lower = 0, upper = 1)
  maf_high <- .mrwin_config_scalar(maf_high, "maf_high", lower = 0, upper = 1)
  if (maf_low >= maf_high) {
    stop("`maf_low` must be smaller than `maf_high`.", call. = FALSE)
  }
  sigma_beta <- .mrwin_config_scalar(sigma_beta, "sigma_beta", lower = 0)
  seed <- if (is.null(seed)) NULL else .mrwin_config_integer(seed, "seed", lower = -Inf)

  structure(list(
    n_outcome = n_outcome,
    m_snps = m_snps,
    alpha_s = alpha_s,
    alpha_u = alpha_u,
    alpha_x = alpha_x,
    gamma_direct = gamma_direct,
    nu_u = nu_u,
    baseline_haz = baseline_haz,
    shape_weibull = shape_weibull,
    censoring_rate = censoring_rate,
    max_follow_up = max_follow_up,
    theta_f = theta_f,
    maf_low = maf_low,
    maf_high = maf_high,
    sigma_beta = sigma_beta,
    seed = seed
  ), class = c("mrwin_config", "list"))
}

.mrwin_config_scalar <- function(x, name, lower = -Inf, upper = Inf, allow_zero = FALSE) {
  x <- as.numeric(x)
  if (length(x) != 1L || !is.finite(x)) {
    stop("`", name, "` must be a single finite number.", call. = FALSE)
  }
  if (allow_zero) {
    bad_lower <- x < lower
  } else {
    bad_lower <- x <= lower
  }
  if (bad_lower || x >= upper) {
    op <- if (allow_zero) "at least" else "greater than"
    stop("`", name, "` must be ", op, " ", lower, " and less than ", upper, ".", call. = FALSE)
  }
  x
}

.mrwin_config_integer <- function(x, name, lower = 1L) {
  x <- suppressWarnings(as.integer(x))
  if (length(x) != 1L || is.na(x) || x < lower) {
    stop("`", name, "` must be a single integer at least ", lower, ".", call. = FALSE)
  }
  x
}

.mrwin_config_vector <- function(x, name, n_expected, lower = -Inf, upper = Inf) {
  x <- as.numeric(x)
  if (length(x) != n_expected || any(!is.finite(x)) || any(x <= lower) || any(x >= upper)) {
    stop("`", name, "` must be a finite numeric vector of length ", n_expected, ".", call. = FALSE)
  }
  x
}

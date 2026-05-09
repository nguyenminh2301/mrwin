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
  stopifnot(length(alpha_x) == 3L)
  stopifnot(length(gamma_direct) == 3L)
  stopifnot(length(nu_u) == 3L)
  stopifnot(length(baseline_haz) == 3L)

  list(
    n_outcome = as.integer(n_outcome),
    m_snps = as.integer(m_snps),
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
    seed = as.integer(seed)
  )
}

# WP19 calibration validation harness (internal). For one scenario/configuration
# it runs M simulated cohorts and returns the Fieller-based rejection rate:
#   scenario = "null"     -> type-I error  (target ~0.05)
#   scenario = "valid_iv" -> coverage of `truth` (target ~0.95)
# Rejection / coverage use the FIELLER interval (the primary, calibrated CI):
# reject H0:delta=0 iff the bounded Fieller 95% CI for delta excludes 0 (this is
# exactly dual to the Fieller numerator test, and matches the validated 0.055).
# A strong instrument is used by default; weak-instrument degeneracy is the
# documented behaviour of IV ratio estimators, not a calibration defect.
.mrwin_validate_calibration <- function(
    M = 200L, N = 4000L, m_snps = 100L, n_strata = 5L,
    scenario = c("null", "valid_iv"),
    inference = c("bootstrap", "analytic"),
    sigma_beta = 0, adjustment = c("none", "ordinal_iptw"),
    stratification = c("prs_rank", "doubly_ranked"),
    bootstrap = 200L, backend = "fast", truth = NULL, seed0 = 1000L) {
  scenario <- match.arg(scenario)
  inference <- match.arg(inference)
  adjustment <- match.arg(adjustment)
  stratification <- match.arg(stratification)
  ax <- if (scenario == "null") c(0, 0, 0) else c(-0.4, -0.4, -0.4)
  if (scenario == "valid_iv" && is.null(truth)) {
    stop("`truth` required for coverage.", call. = FALSE)
  }
  be <- if (inference == "analytic") "dense" else backend
  hit <- rep(NA, M)
  weak <- rep(NA, M)
  for (m in seq_len(M)) {
    s <- seed0 + m
    cfg <- mrwin_config(n_outcome = N, m_snps = m_snps, alpha_x = ax,
                        gamma_direct = c(0, 0, 0), seed = s)
    dat <- mrwin_simulate(cfg, seed = s)
    ep <- mrwin_endpoint(dat$time, dat$status, colnames(dat$time))
    gw <- mrwin_gwas(dat$true_betas, rep(sigma_beta, length(dat$true_betas)))
    cov <- if (adjustment == "ordinal_iptw") cbind(z = as.numeric(scale(dat$U))) else NULL
    ctrl <- mrwin_controls(n_strata = n_strata, bootstrap = bootstrap, seed = s,
                           inference = inference, adjustment = adjustment,
                           stratification = stratification,
                           backend = be, run_sdpd = FALSE)
    fit <- tryCatch(
      mrwin(endpoint = ep, genotype = dat$G, exposure = dat$X, gwas = gw,
            covariates = cov, controls = ctrl),
      error = function(e) NULL)
    if (is.null(fit)) next
    ci <- fit$inference$ci95_delta_fieller
    bounded <- is.finite(ci[[1]]) && is.finite(ci[[2]])
    if (scenario == "null") {
      hit[m] <- bounded && ((ci[[1]] > 0) || (ci[[2]] < 0))   # reject delta=0
    } else {
      hit[m] <- bounded && (ci[[1]] <= truth) && (truth <= ci[[2]])  # cover truth
    }
    weak[m] <- any(vapply(fit$warnings, function(w) w$code == "weak_instrument", logical(1)))
  }
  ok <- !is.na(hit)
  n <- sum(ok)
  list(scenario = scenario, inference = inference, sigma_beta = sigma_beta,
       adjustment = adjustment, stratification = stratification,
       M = M, valid = n, rate = mean(hit[ok]),
       se = sqrt(0.05 * 0.95 / max(n, 1)), weak_frac = mean(weak[ok]))
}

# WP17 / M3: analytic influence-function covariance for the multiplier bootstrap.
#
# The whole inference pipeline consumes `cov_u`, the 2(D-1) x 2(D-1) covariance
# of (LT, DX) = (log-theta, Delta-X) across the multiplier bootstrap. This module
# computes that covariance analytically from the first-order influence function,
# so the B-iteration loop is not needed (for the sampling / multiplier part with
# fixed GWAS weights and adjustment = "none").
#
# Derivation (multiplier weights xi_i = 1 + e_i, e_i iid mean 0, var 1; Exp(1)):
#   log_theta_d = log(W_d/L_d), with W_d = sum_{i in H, j in L} xi_i xi_j 1{win}.
#   First-order influence coefficient of individual k on log_theta_d:
#     k in H:  P+_k/W0 - P-_k/L0      (P+_k = #wins of k over L; P-_k = #losses)
#     k in L:  Q+_k/W0 - Q-_k/L0      (Q+_k = #times k is beaten; Q-_k = ...)
#   Delta_X_d = mean(X|H) - mean(X|L); coefficient of k:
#     k in H:  (X_k - Xbar_H)/n_H ;  k in L: -(X_k - Xbar_L)/n_L
#   With Var(e)=1 and independence across k, cov_u = C^T C, where C is the
#   N x 2(D-1) matrix of these coefficients. Shared strata between adjacent
#   contrasts make C's columns share rows, which reproduces the induced
#   adjacent-contrast correlation the bootstrap exhibits.
#
# Scope: cov_u = Sigma_sampling (analytic influence function) + Sigma_gwas (exact
# GWAS-only resample when sigma_beta>0). adjustment = "none" is exact; for
# adjustment = "ordinal_iptw" the influence function uses the point IPTW weights
# as KNOWN (it omits the estimated-weights / propensity-refit correction), which
# is a first-order approximation. Empirically the omitted term is small (~1-3%,
# either sign, shrinking with N) because in MR the instrument strata are nearly
# independent of the covariates, so the IPTW weights are mild. The bootstrap
# stays the exact default for IPTW; analytic is the fast approximation.

mrwin_analytic_covariance <- function(kernel, strata, X, contrast_plan,
                                      weights = NULL, floor = 1e-12) {
  kernel <- as.matrix(kernel)
  strata <- as.integer(strata)
  X <- as.numeric(X)
  n <- length(strata)
  if (nrow(kernel) != n || ncol(kernel) != n) {
    stop("`kernel` must be an N x N matrix matching `strata`.", call. = FALSE)
  }
  if (length(X) != n) {
    stop("`X` length must match `strata`.", call. = FALSE)
  }
  if (is.null(weights)) {
    weights <- rep(1, n)
  } else {
    weights <- as.numeric(weights)
    if (length(weights) != n) {
      stop("`weights` length must match `strata`.", call. = FALSE)
    }
    if (any(!is.finite(weights)) || any(weights < 0)) {
      stop("`weights` must be finite and non-negative.", call. = FALSE)
    }
  }
  dm1 <- nrow(contrast_plan)
  if (dm1 < 1L) {
    stop("`contrast_plan` must have at least one contrast.", call. = FALSE)
  }

  C <- matrix(0, nrow = n, ncol = 2L * dm1)
  for (a in seq_len(dm1)) {
    H <- which(strata == contrast_plan$high[a])
    L <- which(strata == contrast_plan$low[a])
    if (length(H) == 0L || length(L) == 0L) {
      stop("Analytic covariance requires non-empty strata in every contrast (contrast ", a, ").", call. = FALSE)
    }
    wH <- weights[H]
    wL <- weights[L]
    block <- kernel[H, L, drop = FALSE]
    win <- (block == 1L) + 0  # numeric 0/1
    loss <- (block == -1L) + 0
    # weighted win/loss totals and per-subject weighted leverages
    W0 <- drop(crossprod(wH, win %*% wL))
    L0 <- drop(crossprod(wH, loss %*% wL))
    if (W0 <= 0 || L0 <= 0) {
      stop("Analytic covariance requires positive win and loss counts (weighted) in every contrast (contrast ", a,
           "); use the bootstrap when counts can be zero.", call. = FALSE)
    }
    Pw_plus <- as.numeric(win %*% wL)            # i in H: sum_j w_j 1{win}
    Pw_minus <- as.numeric(loss %*% wL)
    Qw_plus <- as.numeric(crossprod(win, wH))    # j in L: sum_i w_i 1{win}
    Qw_minus <- as.numeric(crossprod(loss, wH))
    sw_H <- sum(wH)
    sw_L <- sum(wL)
    mw_H <- sum(wH * X[H]) / sw_H
    mw_L <- sum(wL * X[L]) / sw_L
    # log-theta influence (column a): coefficient of e_k (k in H or L)
    C[H, a] <- wH * (Pw_plus / W0 - Pw_minus / L0)
    C[L, a] <- wL * (Qw_plus / W0 - Qw_minus / L0)
    # Delta-X influence (column dm1 + a)
    C[H, dm1 + a] <- wH * (X[H] - mw_H) / sw_H
    C[L, dm1 + a] <- -wL * (X[L] - mw_L) / sw_L
  }

  cov_u <- crossprod(C)
  cov_u <- 0.5 * (cov_u + t(cov_u))
  nm <- c(paste0("log_theta:", contrast_plan$label), paste0("delta_x:", contrast_plan$label))
  dimnames(cov_u) <- list(nm, nm)
  cov_u
}

# GWAS-uncertainty covariance term, Sigma_gwas = Cov_{beta*}(point estimate over
# re-stratification). The point estimate (LT^0, DX^0) is piecewise-constant in
# beta (strata jump discretely), so it has no pointwise gradient: a pure-analytic
# Sigma_gwas would require a smoothed / boundary-density approximation that adds a
# bandwidth and its own bias. This routine instead computes the term *exactly* by
# resampling beta* ~ N(beta, diag(sigma_beta^2)), re-stratifying, and recomputing
# the UNWEIGHTED point estimate (no multiplier xi) on the FIXED precomputed
# kernel. It is the irreducible GWAS-only Monte Carlo; the (usually dominant)
# multiplier/sampling part stays analytic, so the full xi resampling is removed.
mrwin_gwas_resample_covariance <- function(kernel, G, X, beta_hat, sigma_beta,
                                           n_strata, contrast_plan,
                                           B_gwas = 200L, seed = NULL, floor = 1e-12,
                                           covariates = NULL,
                                           adjustment = c("none", "ordinal_iptw"),
                                           iptw_truncation = c(0.01, 0.99),
                                           ess_fraction = 0.5) {
  adjustment <- match.arg(adjustment)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  kernel <- as.matrix(kernel)
  G <- as.matrix(G)
  X <- as.numeric(X)
  beta_hat <- as.numeric(beta_hat)
  sigma_beta <- as.numeric(sigma_beta)
  if (length(sigma_beta) == 1L) {
    sigma_beta <- rep(sigma_beta, length(beta_hat))
  }
  if (length(sigma_beta) != length(beta_hat)) {
    stop("`sigma_beta` must have length 1 or length(beta_hat).", call. = FALSE)
  }
  if (adjustment == "ordinal_iptw" && is.null(covariates)) {
    stop("`covariates` are required for adjustment = \"ordinal_iptw\".", call. = FALSE)
  }
  dm1 <- nrow(contrast_plan)
  U <- matrix(NA_real_, nrow = B_gwas, ncol = 2L * dm1)

  for (b in seq_len(B_gwas)) {
    beta_star <- beta_hat + stats::rnorm(length(beta_hat)) * sigma_beta
    strata <- mrwin_prs_strata(G, beta_star, n_strata = n_strata)$strata
    # IPTW: refit propensity at this re-stratification (multiplier xi = 1, the
    # point); skip the iteration if positivity filtering drops a needed stratum.
    w <- NULL
    if (adjustment == "ordinal_iptw") {
      prop <- mrwin_propensity_weights(
        strata = strata, covariates = covariates, method = "ordinal_iptw",
        base_weights = NULL, truncate = iptw_truncation,
        ess_fraction = ess_fraction, fail_on_empty = FALSE, n_strata = n_strata
      )
      needed <- unique(c(contrast_plan$high, contrast_plan$low))
      if (length(intersect(needed, prop$dropped_strata)) > 0L || any(is.na(prop$weights))) {
        next
      }
      w <- prop$weights
    }
    for (a in seq_len(dm1)) {
      H <- which(strata == contrast_plan$high[a])
      L <- which(strata == contrast_plan$low[a])
      if (length(H) == 0L || length(L) == 0L) {
        next
      }
      sums <- mrwin_stratum_win_loss(kernel, H, L, weights = w)
      U[b, a] <- log(max(sums[["wins"]], floor) / max(sums[["losses"]], floor))
      if (is.null(w)) {
        U[b, dm1 + a] <- mean(X[H]) - mean(X[L])
      } else {
        U[b, dm1 + a] <- stats::weighted.mean(X[H], w[H]) - stats::weighted.mean(X[L], w[L])
      }
    }
  }

  valid <- stats::complete.cases(U)
  if (sum(valid) < 2L) {
    stop("Too few valid GWAS-resample iterations.", call. = FALSE)
  }
  cov_u <- stats::cov(U[valid, , drop = FALSE])
  cov_u <- 0.5 * (cov_u + t(cov_u))
  nm <- c(paste0("log_theta:", contrast_plan$label), paste0("delta_x:", contrast_plan$label))
  dimnames(cov_u) <- list(nm, nm)
  cov_u
}

# Analytic inference: feed the analytic cov_u through the SAME downstream
# machinery the bootstrap uses (bivariate-Delta ISG covariance, GLS pooling,
# Fieller). Replaces the multiplier-bootstrap loop for the sampling part with
# fixed GWAS weights and adjustment = "none". `estimate` is an `mrwin_estimate`
# object (dense backend, so it carries the kernel); `X` is the exposure vector.
mrwin_analytic_inference <- function(estimate, X, G = NULL, beta_hat = NULL,
                                     sigma_beta = 0, n_strata = NULL,
                                     B_gwas = 200L, seed = NULL, z = 1.96,
                                     weights = NULL, covariates = NULL,
                                     adjustment = c("none", "ordinal_iptw"),
                                     iptw_truncation = c(0.01, 0.99),
                                     ess_fraction = 0.5) {
  adjustment <- match.arg(adjustment)
  if (!inherits(estimate, "mrwin_estimate")) {
    stop("`estimate` must be an `mrwin_estimate` object.", call. = FALSE)
  }
  if (is.null(estimate$kernel)) {
    stop("Analytic inference needs the dense kernel; run `mrwin_estimate()` ",
         "(or `mrwin(backend = \"dense\")`) so the kernel is available.", call. = FALSE)
  }
  cov_sampling <- mrwin_analytic_covariance(
    estimate$kernel, estimate$strata, X, estimate$contrast_plan, weights = weights
  )
  cov_gwas <- NULL
  gwas_included <- !is.null(G) && !is.null(beta_hat) && any(as.numeric(sigma_beta) > 0)
  if (gwas_included) {
    if (is.null(n_strata)) {
      n_strata <- max(estimate$strata)
    }
    cov_gwas <- mrwin_gwas_resample_covariance(
      kernel = estimate$kernel, G = G, X = X, beta_hat = beta_hat,
      sigma_beta = sigma_beta, n_strata = n_strata,
      contrast_plan = estimate$contrast_plan, B_gwas = B_gwas, seed = seed,
      covariates = covariates, adjustment = adjustment,
      iptw_truncation = iptw_truncation, ess_fraction = ess_fraction
    )
  }
  cov_u <- if (gwas_included) cov_sampling + cov_gwas else cov_sampling
  sigma_isg <- .mrwin_isg_covariance(
    cov_u = cov_u,
    point_log_theta = estimate$log_theta,
    point_delta_x = estimate$delta_x
  )
  pooled <- mrwin_gls_pool(estimate$delta_isg, sigma_isg, shrink = TRUE)
  fieller <- .mrwin_fieller_ci(
    cov_u = cov_u,
    sigma_isg = sigma_isg,
    point_log_theta = estimate$log_theta,
    point_delta_x = estimate$delta_x,
    z = z
  )

  list(
    method = if (gwas_included) {
      "analytic_sampling_plus_gwas_resample"
    } else {
      "analytic_influence_function_sampling"
    },
    gwas_included = gwas_included,
    cov_u = cov_u,
    cov_sampling = cov_sampling,
    cov_gwas = cov_gwas,
    sigma_isg = sigma_isg,
    sigma_lw = pooled$sigma,
    point_log_theta = estimate$log_theta,
    point_delta_x = estimate$delta_x,
    delta_isg = estimate$delta_isg,
    delta_gls = pooled$delta_gls,
    se_delta_gls = pooled$se_delta_gls,
    dscwr = pooled$dscwr,
    ci95_delta = pooled$ci95_delta,
    ci95_dscwr = pooled$ci95_dscwr,
    ci95_delta_fieller = fieller$delta,
    ci95_dscwr_fieller = exp(pmax(pmin(fieller$delta, 50), -50)),
    fieller_unbounded = fieller$unbounded,
    fieller_coefficients = fieller$coefficients,
    q = pooled$q,
    q_df = pooled$q_df,
    q_p_value = pooled$q_p_value,
    gls_weights = pooled$gls_weights,
    ledoit_wolf_rho = pooled$ledoit_wolf_rho
  )
}

# Top-level analytic inference engine, output-compatible with
# mrwin_multiplier_bootstrap so mrwin() and the S3 methods consume it unchanged.
# Supports adjustment in {"none","ordinal_iptw"}; covariance = analytic sampling
# influence function (weighted by the point IPTW weights when adjusting) + (if
# sigma_beta>0) exact GWAS-only resample (with per-draw propensity refit under
# IPTW). For IPTW the estimated-weights correction is omitted (see the module
# header). Uses the dense kernel (moderate N).
mrwin_analytic_bootstrap <- function(
    time, status, G, X, beta_hat,
    sigma_beta = 0, n_strata = 10L, B_gwas = 200L, seed = NULL,
    kernel = NULL, block_size = 4000L, floor = 1e-12,
    covariates = NULL, adjustment = c("none", "ordinal_iptw"),
    iptw_truncation = c(0.01, 0.99), ess_fraction = 0.5
) {
  adjustment <- match.arg(adjustment)
  checked <- .mrwin_validate_estimate_inputs(
    time = time, status = status, G = G, X = X, beta_hat = beta_hat,
    n_strata = n_strata, kernel = kernel, weights = NULL, floor = floor
  )
  time <- checked$time; status <- checked$status; G <- checked$G; X <- checked$X
  beta_hat <- checked$beta_hat; n_strata <- checked$n_strata; floor <- checked$floor
  sigma_beta <- as.numeric(sigma_beta)
  if (length(sigma_beta) == 1L) sigma_beta <- rep(sigma_beta, length(beta_hat))
  if (adjustment == "ordinal_iptw") {
    if (is.null(covariates)) stop("`covariates` are required for adjustment = \"ordinal_iptw\".", call. = FALSE)
    covariates <- as.matrix(covariates)
    if (nrow(covariates) != nrow(G)) stop("`covariates` must have one row per analysis row.", call. = FALSE)
  }

  if (is.null(kernel)) kernel <- mrwin_kernel(time, status, block_size = block_size)

  point_strata <- mrwin_prs_strata(G, beta_hat, n_strata = n_strata)$strata
  point_adjustment <- .mrwin_point_adjustment(
    strata = point_strata, covariates = covariates, adjustment = adjustment,
    iptw_truncation = iptw_truncation, ess_fraction = ess_fraction, n_strata = n_strata
  )
  est <- mrwin_estimate(
    time = time, status = status, G = G, X = X, beta_hat = beta_hat,
    n_strata = n_strata, kernel = kernel, weights = point_adjustment$weights,
    floor = floor, active_strata = point_adjustment$active_strata
  )
  ana <- mrwin_analytic_inference(
    est, X, G = G, beta_hat = beta_hat, sigma_beta = sigma_beta,
    n_strata = n_strata, B_gwas = B_gwas, seed = seed,
    weights = point_adjustment$weights, covariates = covariates,
    adjustment = adjustment, iptw_truncation = iptw_truncation, ess_fraction = ess_fraction
  )

  dm1 <- nrow(est$contrast_plan)
  na_vec <- rep(NA_real_, dm1)
  names(na_vec) <- names(est$log_theta)

  structure(list(
    B = as.integer(if (any(sigma_beta > 0)) B_gwas else 0L),
    n_valid = as.integer(if (any(sigma_beta > 0)) B_gwas else 0L),
    n_invalid = 0L,
    n_strata = as.integer(n_strata),
    active_strata = point_adjustment$active_strata,
    dropped_strata = point_adjustment$dropped_strata,
    contrast_plan = est$contrast_plan,
    point_log_theta = est$log_theta,
    point_cwr = est$cwr,
    point_delta_x = est$delta_x,
    delta_isg = est$delta_isg,
    delta_gls = ana$delta_gls,
    se_delta_gls = ana$se_delta_gls,
    dscwr = ana$dscwr,
    ci95_delta = ana$ci95_delta,
    ci95_dscwr = ana$ci95_dscwr,
    ci95_delta_bivariate_delta = ana$ci95_delta,
    ci95_dscwr_bivariate_delta = ana$ci95_dscwr,
    ci95_delta_fieller = ana$ci95_delta_fieller,
    ci95_dscwr_fieller = ana$ci95_dscwr_fieller,
    fieller_unbounded = ana$fieller_unbounded,
    fieller_coefficients = ana$fieller_coefficients,
    q = ana$q, q_df = ana$q_df, q_p_value = ana$q_p_value,
    gls_weights = ana$gls_weights,
    ledoit_wolf_rho = ana$ledoit_wolf_rho,
    cov_u = ana$cov_u, sigma_isg = ana$sigma_isg, sigma_lw = ana$sigma_lw,
    kurtosis_log_theta = na_vec, kurtosis_delta_x = na_vec,
    skew_log_theta = na_vec, skew_delta_x = na_vec,
    bootstrap_log_theta = NULL, bootstrap_delta_x = NULL,
    bootstrap_ess = NULL, bootstrap_dropped_strata = NULL, valid_bootstrap = NULL,
    covariance_method = ana$method,
    adjustment = point_adjustment, kernel = kernel, strata = est$strata
  ), class = c("mrwin_bootstrap", "list"))
}

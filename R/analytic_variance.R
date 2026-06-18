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
# Scope (this step): sampling/multiplier part, adjustment = "none". GWAS-draw
# (re-stratification) and IPTW terms are documented follow-ups (WP17 step B).

mrwin_analytic_covariance <- function(kernel, strata, X, contrast_plan, floor = 1e-12) {
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
    block <- kernel[H, L, drop = FALSE]
    win <- block == 1L
    loss <- block == -1L
    W0 <- sum(win)
    L0 <- sum(loss)
    if (W0 <= 0 || L0 <= 0) {
      stop("Analytic covariance requires positive win and loss counts in every contrast (contrast ", a,
           "); use the bootstrap when counts can be zero.", call. = FALSE)
    }
    nH <- length(H)
    nL <- length(L)
    # log-theta influence (column a)
    C[H, a] <- rowSums(win) / W0 - rowSums(loss) / L0
    C[L, a] <- colSums(win) / W0 - colSums(loss) / L0
    # Delta-X influence (column dm1 + a)
    C[H, dm1 + a] <- (X[H] - mean(X[H])) / nH
    C[L, dm1 + a] <- -(X[L] - mean(X[L])) / nL
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
                                           B_gwas = 200L, seed = NULL, floor = 1e-12) {
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
  dm1 <- nrow(contrast_plan)
  U <- matrix(NA_real_, nrow = B_gwas, ncol = 2L * dm1)

  for (b in seq_len(B_gwas)) {
    beta_star <- beta_hat + stats::rnorm(length(beta_hat)) * sigma_beta
    strata <- mrwin_prs_strata(G, beta_star, n_strata = n_strata)$strata
    for (a in seq_len(dm1)) {
      H <- which(strata == contrast_plan$high[a])
      L <- which(strata == contrast_plan$low[a])
      if (length(H) == 0L || length(L) == 0L) {
        next
      }
      sums <- mrwin_stratum_win_loss(kernel, H, L)
      U[b, a] <- log(max(sums[["wins"]], floor) / max(sums[["losses"]], floor))
      U[b, dm1 + a] <- mean(X[H]) - mean(X[L])
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
                                     B_gwas = 200L, seed = NULL, z = 1.96) {
  if (!inherits(estimate, "mrwin_estimate")) {
    stop("`estimate` must be an `mrwin_estimate` object.", call. = FALSE)
  }
  if (is.null(estimate$kernel)) {
    stop("Analytic inference needs the dense kernel; run `mrwin_estimate()` ",
         "(or `mrwin(backend = \"dense\")`) so the kernel is available.", call. = FALSE)
  }
  cov_sampling <- mrwin_analytic_covariance(
    estimate$kernel, estimate$strata, X, estimate$contrast_plan
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
      contrast_plan = estimate$contrast_plan, B_gwas = B_gwas, seed = seed
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
# Supports adjustment = "none"; covariance = analytic sampling + (if sigma_beta>0)
# exact GWAS-only resample. Uses the dense kernel (moderate N).
mrwin_analytic_bootstrap <- function(
    time, status, G, X, beta_hat,
    sigma_beta = 0, n_strata = 10L, B_gwas = 200L, seed = NULL,
    kernel = NULL, block_size = 4000L, floor = 1e-12
) {
  checked <- .mrwin_validate_estimate_inputs(
    time = time, status = status, G = G, X = X, beta_hat = beta_hat,
    n_strata = n_strata, kernel = kernel, weights = NULL, floor = floor
  )
  time <- checked$time; status <- checked$status; G <- checked$G; X <- checked$X
  beta_hat <- checked$beta_hat; n_strata <- checked$n_strata; floor <- checked$floor
  sigma_beta <- as.numeric(sigma_beta)
  if (length(sigma_beta) == 1L) sigma_beta <- rep(sigma_beta, length(beta_hat))

  if (is.null(kernel)) kernel <- mrwin_kernel(time, status, block_size = block_size)

  point_strata <- mrwin_prs_strata(G, beta_hat, n_strata = n_strata)$strata
  point_adjustment <- .mrwin_point_adjustment(
    strata = point_strata, covariates = NULL, adjustment = "none",
    iptw_truncation = c(0.01, 0.99), ess_fraction = 0.5, n_strata = n_strata
  )
  est <- mrwin_estimate(
    time = time, status = status, G = G, X = X, beta_hat = beta_hat,
    n_strata = n_strata, kernel = kernel, weights = point_adjustment$weights,
    floor = floor, active_strata = point_adjustment$active_strata
  )
  ana <- mrwin_analytic_inference(
    est, X, G = G, beta_hat = beta_hat, sigma_beta = sigma_beta,
    n_strata = n_strata, B_gwas = B_gwas, seed = seed
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

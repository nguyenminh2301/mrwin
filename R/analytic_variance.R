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

# Analytic inference: feed the analytic cov_u through the SAME downstream
# machinery the bootstrap uses (bivariate-Delta ISG covariance, GLS pooling,
# Fieller). Replaces the multiplier-bootstrap loop for the sampling part with
# fixed GWAS weights and adjustment = "none". `estimate` is an `mrwin_estimate`
# object (dense backend, so it carries the kernel); `X` is the exposure vector.
mrwin_analytic_inference <- function(estimate, X, z = 1.96) {
  if (!inherits(estimate, "mrwin_estimate")) {
    stop("`estimate` must be an `mrwin_estimate` object.", call. = FALSE)
  }
  if (is.null(estimate$kernel)) {
    stop("Analytic inference needs the dense kernel; run `mrwin_estimate()` ",
         "(or `mrwin(backend = \"dense\")`) so the kernel is available.", call. = FALSE)
  }
  cov_u <- mrwin_analytic_covariance(
    estimate$kernel, estimate$strata, X, estimate$contrast_plan
  )
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
    method = "analytic_influence_function_sampling",
    cov_u = cov_u,
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
    q = pooled$q,
    q_df = pooled$q_df,
    q_p_value = pooled$q_p_value,
    gls_weights = pooled$gls_weights,
    ledoit_wolf_rho = pooled$ledoit_wolf_rho
  )
}

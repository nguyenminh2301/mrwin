# Two-sample / summary-data win-ratio MR.
#
# The outcome-side summary statistic for a priority-ranked (hierarchical) outcome
# is a per-SNP WIN-ODDS coefficient: the log-win-odds slope of the outcome on
# allele dosage. mrwin_win_snp() estimates it with a closed-form standard error,
# mrwin_win_gwas() runs it across a panel of SNPs, and mrwin_twosample_ivw() pools
# the resulting (beta_GX, delta_winodds) pairs by inverse-variance weighting --
# the win-ratio analogue of ordinary two-sample MR. The estimator is consistent
# for the causal win-odds gradient (dev/findings-two-sample-winratio.md); its SE
# is the influence-function variance of log(W/L) used by mrwin_analytic_covariance,
# evaluated subquadratically (point estimate on the fast kernel, variance from
# per-subject counts on a capped subject subsample).

# log-win-odds influence-function variance for one high-vs-low contrast.
# coef_k = w_k (P+_k / W0 - P-_k / L0); Var(log theta) = sum_k coef_k^2, estimated
# as (n/|S|) sum over a uniform subject subsample S of each group (exact when the
# group is not subsampled). Returns Var(log theta).
.mrwin_logwinodds_var <- function(Ht, Hs, wH, Lt, Ls, wL, W0, L0,
                                  max_subjects, rng_seed) {
  nh <- nrow(Ht); nl <- nrow(Lt)
  subj <- getFromNamespace("mrwin_subject_win_loss_cpp", "mrwin")
  if (!is.null(rng_seed)) {
    old <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
    set.seed(rng_seed)
    on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
  }
  iH <- if (nh > max_subjects) sort(sample.int(nh, max_subjects)) else seq_len(nh)
  iL <- if (nl > max_subjects) sort(sample.int(nl, max_subjects)) else seq_len(nl)
  # high subjects vs the FULL low group: (P+_k, P-_k) weighted by wL
  cH <- subj(Ht[iH, , drop = FALSE], Hs[iH, , drop = FALSE], Lt, Ls, wL)
  coefH <- wH[iH] * (cH[, 1] / W0 - cH[, 2] / L0)
  # low subjects vs the FULL high group: out[,1] = low beats high (= high loss, Q-_j),
  # out[,2] = high beats low (= high win, Q+_j); weighted by wH
  cL <- subj(Lt[iL, , drop = FALSE], Ls[iL, , drop = FALSE], Ht, Hs, wH)
  coefL <- wL[iL] * (cL[, 2] / W0 - cL[, 1] / L0)
  (nh / length(iH)) * sum(coefH^2) + (nl / length(iL)) * sum(coefL^2)
}

#' Per-SNP win-odds coefficient with a closed-form standard error
#'
#' Estimates the log-win-odds slope of a priority-ranked (hierarchical) outcome on
#' allele dosage for a single SNP -- the outcome-side summary statistic for
#' two-sample win-ratio Mendelian randomization. The point estimate uses the
#' subquadratic hierarchical kernel; the standard error is the influence-function
#' variance of the log-win-odds, evaluated on a capped subject subsample so the
#' whole routine stays near-linear in sample size.
#'
#' @param time,status N x K matrices of event times and 0/1 status, in priority
#'   order (column 1 = highest priority), as elsewhere in the package.
#' @param dosage Numeric allele dosage (length N); any scale is allowed.
#' @param weights Optional non-negative subject weights (length N).
#' @param split How to dichotomize dosage into the high/low groups whose win-odds
#'   is contrasted: `"median"` (default) or `"threshold"`.
#' @param threshold Cut value for `split = "threshold"` (high = dosage > threshold).
#' @param max_subjects Subject-subsample cap per group for the variance (default
#'   2000); the point estimate always uses all subjects.
#' @param seed Optional RNG seed for the variance subsample (reproducibility).
#' @param floor Minimum weighted win/loss count; below it the estimate is `NA`.
#' @return A list with `delta` (win-odds slope per dosage unit) and its `se`,
#'   plus `log_theta`, `se_log_theta`, `d_diff` (mean-dosage gap between groups),
#'   `wins`, `losses`, `n_high`, `n_low`.
#' @export
mrwin_win_snp <- function(time, status, dosage, weights = NULL,
                          split = c("median", "threshold"), threshold = NULL,
                          max_subjects = 2000L, seed = NULL, floor = 1e-9) {
  time <- as.matrix(time); status <- as.matrix(status)
  dosage <- as.numeric(dosage)
  n <- length(dosage)
  if (nrow(time) != n || nrow(status) != n) {
    stop("`time`/`status` rows must match length(dosage).", call. = FALSE)
  }
  if (is.null(weights)) {
    weights <- rep(1, n)
  } else {
    weights <- as.numeric(weights)
    if (length(weights) != n) stop("`weights` must have length N.", call. = FALSE)
    if (any(!is.finite(weights)) || any(weights < 0)) {
      stop("`weights` must be finite and non-negative.", call. = FALSE)
    }
  }
  split <- match.arg(split)
  cut <- if (split == "median") stats::median(dosage) else {
    if (is.null(threshold)) stop("`threshold` is required for split = \"threshold\".", call. = FALSE)
    as.numeric(threshold)
  }
  hi <- which(dosage > cut)
  lo <- which(dosage <= cut)
  na_out <- list(delta = NA_real_, se = NA_real_, log_theta = NA_real_,
                 se_log_theta = NA_real_, d_diff = NA_real_, wins = NA_real_,
                 losses = NA_real_, n_high = length(hi), n_low = length(lo))
  if (length(hi) < 2L || length(lo) < 2L) return(na_out)
  d_diff <- mean(dosage[hi]) - mean(dosage[lo])
  if (!is.finite(d_diff) || abs(d_diff) < .Machine$double.eps) return(na_out)

  Ht <- time[hi, , drop = FALSE]; Hs <- status[hi, , drop = FALSE]; wH <- weights[hi]
  Lt <- time[lo, , drop = FALSE]; Ls <- status[lo, , drop = FALSE]; wL <- weights[lo]
  wl <- mrwin_fast_pair_win_loss(Ht, Hs, Lt, Ls, weights_high = wH, weights_low = wL)
  W0 <- as.numeric(wl[["wins"]]); L0 <- as.numeric(wl[["losses"]])
  if (W0 <= floor || L0 <= floor) return(na_out)

  log_theta <- log(W0 / L0)
  var_lt <- .mrwin_logwinodds_var(Ht, Hs, wH, Lt, Ls, wL, W0, L0,
                                  max_subjects = max_subjects, rng_seed = seed)
  se_lt <- sqrt(var_lt)
  list(delta = log_theta / d_diff, se = se_lt / abs(d_diff),
       log_theta = log_theta, se_log_theta = se_lt, d_diff = d_diff,
       wins = W0, losses = L0, n_high = length(hi), n_low = length(lo))
}

#' Win-odds GWAS: per-SNP win-odds coefficients across a SNP panel
#'
#' Applies [mrwin_win_snp()] to each column of a genotype matrix, producing the
#' shareable outcome-side summary statistics for two-sample win-ratio MR.
#'
#' @inheritParams mrwin_win_snp
#' @param G N x M genotype/dosage matrix (one column per SNP).
#' @param ... Passed to [mrwin_win_snp()].
#' @return A data frame with one row per SNP: `snp`, `delta`, `se`, `log_theta`,
#'   `se_log_theta`, `d_diff`, `n_high`, `n_low`.
#' @export
mrwin_win_gwas <- function(time, status, G, weights = NULL, max_subjects = 2000L,
                           seed = NULL, ...) {
  G <- as.matrix(G); m <- ncol(G)
  res <- lapply(seq_len(m), function(k) {
    mrwin_win_snp(time, status, G[, k], weights = weights,
                  max_subjects = max_subjects, seed = seed, ...)
  })
  data.frame(
    snp = seq_len(m),
    delta = vapply(res, `[[`, numeric(1), "delta"),
    se = vapply(res, `[[`, numeric(1), "se"),
    log_theta = vapply(res, `[[`, numeric(1), "log_theta"),
    se_log_theta = vapply(res, `[[`, numeric(1), "se_log_theta"),
    d_diff = vapply(res, `[[`, numeric(1), "d_diff"),
    n_high = vapply(res, `[[`, numeric(1), "n_high"),
    n_low = vapply(res, `[[`, numeric(1), "n_low")
  )
}

#' Two-sample win-ratio MR by inverse-variance weighting
#'
#' Pools per-SNP exposure coefficients `beta_gx` (from an exposure GWAS) with
#' per-SNP win-odds coefficients `delta_gy` (from [mrwin_win_gwas()]) into a single
#' causal win-odds gradient, using fixed-effect IVW with real inverse-variance
#' weights. By default the exposure coefficients are treated as fixed (the usual
#' two-sample convention); supplying `se_gx` adds the second-order Wald-ratio
#' variance term.
#'
#' @param beta_gx Per-SNP exposure associations.
#' @param delta_gy Per-SNP win-odds coefficients (e.g. `mrwin_win_gwas()$delta`).
#' @param se_gy Standard errors of `delta_gy` (e.g. `mrwin_win_gwas()$se`).
#' @param se_gx Optional standard errors of `beta_gx`.
#' @return A list: `gamma` (pooled win-odds gradient per exposure unit), `se`,
#'   `ci95`, `n_snp`, and Cochran's `Q`, `Q_df`, `Q_p` heterogeneity test.
#' @export
mrwin_twosample_ivw <- function(beta_gx, delta_gy, se_gy, se_gx = NULL) {
  beta_gx <- as.numeric(beta_gx); delta_gy <- as.numeric(delta_gy); se_gy <- as.numeric(se_gy)
  ok <- is.finite(beta_gx) & is.finite(delta_gy) & is.finite(se_gy) & se_gy > 0 & beta_gx != 0
  if (!is.null(se_gx)) { se_gx <- as.numeric(se_gx); ok <- ok & is.finite(se_gx) }
  if (sum(ok) < 1L) stop("No usable SNPs (need finite beta_gx, delta_gy, se_gy > 0).", call. = FALSE)
  bx <- beta_gx[ok]; dy <- delta_gy[ok]; sy <- se_gy[ok]
  ratio <- dy / bx
  var_ratio <- sy^2 / bx^2
  if (!is.null(se_gx)) var_ratio <- var_ratio + (dy^2) * (se_gx[ok]^2) / (bx^4)  # delta method
  w <- 1 / var_ratio
  gamma <- sum(w * ratio) / sum(w)
  se <- sqrt(1 / sum(w))
  Q <- sum(w * (ratio - gamma)^2); df <- max(length(bx) - 1L, 0L)
  list(gamma = gamma, se = se, ci95 = gamma + c(-1.96, 1.96) * se, n_snp = length(bx),
       Q = Q, Q_df = df,
       Q_p = if (df > 0L) stats::pchisq(Q, df, lower.tail = FALSE) else NA_real_)
}

mrwin_propensity_weights <- function(
    strata,
    covariates,
    method = c("ordinal_iptw"),
    base_weights = NULL,
    truncate = c(0.01, 0.99),
    ess_fraction = 0.5,
    fail_on_empty = TRUE,
    n_strata = NULL
) {
  method <- match.arg(method)
  strata <- suppressWarnings(as.integer(strata))
  covariates <- as.matrix(covariates)
  if (length(strata) != nrow(covariates)) {
    stop("`strata` length must match the number of covariate rows.", call. = FALSE)
  }
  if (any(is.na(strata)) || any(strata < 1L)) {
    stop("`strata` must contain positive integer stratum labels.", call. = FALSE)
  }
  if (any(!is.finite(covariates))) {
    stop("`covariates` must contain finite numeric values.", call. = FALSE)
  }
  if (ncol(covariates) < 1L) {
    stop("`covariates` must contain at least one column.", call. = FALSE)
  }
  truncate <- .mrwin_validate_truncation(truncate)
  ess_fraction <- .mrwin_validate_ess_fraction(ess_fraction)
  if (is.null(base_weights)) {
    base_weights <- rep(1, length(strata))
  } else {
    base_weights <- as.numeric(base_weights)
    if (length(base_weights) != length(strata)) {
      stop("`base_weights` length must match `strata`.", call. = FALSE)
    }
    if (any(!is.finite(base_weights)) || any(base_weights < 0)) {
      stop("`base_weights` must be finite and non-negative.", call. = FALSE)
    }
  }

  if (is.null(n_strata)) {
    n_strata <- max(strata)
  } else {
    n_strata <- as.integer(n_strata)
    if (length(n_strata) != 1L || is.na(n_strata) || n_strata < max(strata)) {
      stop("`n_strata` must cover all observed stratum labels.", call. = FALSE)
    }
  }
  counts <- tabulate(strata, nbins = n_strata)
  empty <- which(counts == 0L)
  if (length(empty) > 0L) {
    if (fail_on_empty) {
      stop("Cannot fit propensity model with empty strata: ", paste(empty, collapse = ", "), call. = FALSE)
    }
    return(.mrwin_empty_propensity_result(strata, covariates, base_weights, empty, truncate, ess_fraction, method))
  }

  model <- .mrwin_fit_ordinal_propensity(strata, covariates, base_weights)
  probabilities <- model$probabilities
  marginal <- as.numeric(rowsum(base_weights, strata, reorder = FALSE))
  marginal <- marginal / sum(marginal)
  raw_weights <- marginal[strata] / pmax(probabilities[cbind(seq_along(strata), strata)], 1e-12)

  limits <- stats::quantile(raw_weights, probs = truncate, names = FALSE, type = 7, na.rm = TRUE)
  weights <- pmin(pmax(raw_weights, limits[1L]), limits[2L])
  ess <- .mrwin_stratum_ess(strata, weights, n_strata)
  counts <- tabulate(strata, nbins = n_strata)
  ess_threshold <- ess_fraction * counts
  names(ess_threshold) <- paste0("s", seq_len(n_strata))
  dropped <- which(ess < ess_threshold)
  balance <- .mrwin_balance_diagnostics(strata, covariates, weights)

  structure(
    list(
      method = method,
      weights = as.numeric(weights),
      raw_weights = as.numeric(raw_weights),
      probabilities = probabilities,
      marginal = marginal,
      ess = ess,
      ess_threshold = ess_threshold,
      dropped_strata = as.integer(dropped),
      active_strata = setdiff(seq_len(n_strata), dropped),
      truncation = limits,
      truncation_probs = truncate,
      balance = balance,
      model = model
    ),
    class = c("mrwin_propensity_weights", "list")
  )
}

.mrwin_validate_truncation <- function(truncate) {
  truncate <- as.numeric(truncate)
  if (length(truncate) != 2L || any(!is.finite(truncate)) ||
      truncate[1L] < 0 || truncate[2L] > 1 || truncate[1L] > truncate[2L]) {
    stop("`iptw_truncation` must be two probabilities with low <= high.", call. = FALSE)
  }
  truncate
}

.mrwin_validate_ess_fraction <- function(ess_fraction) {
  ess_fraction <- as.numeric(ess_fraction)
  if (length(ess_fraction) != 1L || !is.finite(ess_fraction) ||
      ess_fraction < 0 || ess_fraction > 1) {
    stop("`ess_fraction` must be between 0 and 1.", call. = FALSE)
  }
  ess_fraction
}

.mrwin_empty_propensity_result <- function(strata, covariates, base_weights, empty, truncate, ess_fraction, method) {
  n_strata <- length(tabulate(strata, nbins = max(c(strata, empty))))
  balance <- .mrwin_balance_diagnostics(strata, covariates, rep(1, length(strata)))
  ess_threshold <- ess_fraction * tabulate(strata, nbins = n_strata)
  names(ess_threshold) <- paste0("s", seq_len(n_strata))
  structure(
    list(
      method = method,
      weights = rep(NA_real_, length(strata)),
      raw_weights = rep(NA_real_, length(strata)),
      probabilities = matrix(NA_real_, nrow = length(strata), ncol = n_strata),
      marginal = rep(NA_real_, n_strata),
      ess = rep(NA_real_, n_strata),
      ess_threshold = ess_threshold,
      dropped_strata = as.integer(seq_len(n_strata)),
      active_strata = integer(),
      empty_strata = as.integer(empty),
      truncation = rep(NA_real_, 2L),
      truncation_probs = truncate,
      balance = balance,
      model = NULL
    ),
    class = c("mrwin_propensity_weights", "list")
  )
}

.mrwin_fit_ordinal_propensity <- function(strata, covariates, base_weights) {
  n_strata <- max(strata)
  z <- .mrwin_standardize_covariates(covariates, base_weights)
  p <- ncol(z$x)
  cumulative <- cumsum(tabulate(strata, nbins = n_strata)) / length(strata)
  cut_init <- stats::qlogis(pmin(pmax(cumulative[-n_strata], 1e-4), 1 - 1e-4))
  if (length(cut_init) > 1L) {
    cut_init <- sort(cut_init)
    theta_cut <- c(cut_init[1L], log(pmax(diff(cut_init), 0.05)))
  } else {
    theta_cut <- cut_init
  }
  start <- c(theta_cut, rep(0, p))
  objective <- function(par) {
    cuts <- .mrwin_unpack_cutpoints(par, n_strata)
    beta <- if (p > 0L) par[n_strata:length(par)] else numeric()
    eta <- if (p > 0L) drop(z$x %*% beta) else rep(0, length(strata))
    prob <- .mrwin_ordinal_probabilities(cuts, eta)
    obs <- pmax(prob[cbind(seq_along(strata), strata)], 1e-12)
    -sum(base_weights * log(obs)) + 1e-6 * sum(beta^2)
  }
  fit <- stats::optim(start, objective, method = "BFGS", control = list(maxit = 1000))
  if (!is.finite(fit$value) || fit$convergence != 0L) {
    stop("Ordinal IPTW propensity model failed to converge.", call. = FALSE)
  }
  cuts <- .mrwin_unpack_cutpoints(fit$par, n_strata)
  beta <- if (p > 0L) fit$par[n_strata:length(fit$par)] else numeric()
  eta <- if (p > 0L) drop(z$x %*% beta) else rep(0, length(strata))
  prob <- .mrwin_ordinal_probabilities(cuts, eta)
  colnames(prob) <- paste0("s", seq_len(n_strata))

  list(
    cutpoints = cuts,
    beta = beta,
    center = z$center,
    scale = z$scale,
    probabilities = prob,
    convergence = fit$convergence,
    objective = fit$value
  )
}

.mrwin_standardize_covariates <- function(covariates, base_weights) {
  covariates <- as.matrix(covariates)
  center <- scale <- numeric(ncol(covariates))
  out <- covariates
  total_weight <- sum(base_weights)
  if (total_weight <= 0) {
    stop("`base_weights` must have positive total weight.", call. = FALSE)
  }
  for (j in seq_len(ncol(covariates))) {
    center[j] <- sum(base_weights * covariates[, j]) / total_weight
    variance <- sum(base_weights * (covariates[, j] - center[j])^2) / total_weight
    scale[j] <- sqrt(variance)
    if (!is.finite(scale[j]) || scale[j] <= 0) {
      stop("Covariates used for IPTW must have positive variance.", call. = FALSE)
    }
    out[, j] <- (covariates[, j] - center[j]) / scale[j]
  }
  list(x = out, center = center, scale = scale)
}

.mrwin_unpack_cutpoints <- function(par, n_strata) {
  n_cut <- n_strata - 1L
  if (n_cut == 1L) {
    return(par[1L])
  }
  c(par[1L], par[1L] + cumsum(exp(par[2:n_cut])))
}

.mrwin_ordinal_probabilities <- function(cutpoints, eta) {
  n <- length(eta)
  n_strata <- length(cutpoints) + 1L
  cumulative <- matrix(NA_real_, nrow = n, ncol = length(cutpoints))
  for (k in seq_along(cutpoints)) {
    cumulative[, k] <- stats::plogis(cutpoints[k] - eta)
  }
  prob <- matrix(NA_real_, nrow = n, ncol = n_strata)
  prob[, 1L] <- cumulative[, 1L]
  if (n_strata > 2L) {
    for (k in 2:(n_strata - 1L)) {
      prob[, k] <- cumulative[, k] - cumulative[, k - 1L]
    }
  }
  prob[, n_strata] <- 1 - cumulative[, n_strata - 1L]
  prob <- pmax(prob, 1e-12)
  prob / rowSums(prob)
}

.mrwin_stratum_ess <- function(strata, weights, n_strata = max(strata)) {
  out <- rep(NA_real_, n_strata)
  for (d in seq_len(n_strata)) {
    w <- weights[strata == d]
    out[d] <- if (length(w) == 0L || sum(w^2) <= 0) NA_real_ else sum(w)^2 / sum(w^2)
  }
  names(out) <- paste0("s", seq_len(n_strata))
  out
}

.mrwin_balance_diagnostics <- function(strata, covariates, weights = NULL) {
  covariates <- as.matrix(covariates)
  if (is.null(colnames(covariates))) {
    colnames(covariates) <- paste0("covariate_", seq_len(ncol(covariates)))
  }
  before <- .mrwin_balance_smd(strata, covariates, rep(1, length(strata)))
  after <- if (is.null(weights)) before else .mrwin_balance_smd(strata, covariates, weights)
  data.frame(
    covariate = colnames(covariates),
    before_max_abs_smd = before,
    after_max_abs_smd = after,
    stringsAsFactors = FALSE
  )
}

.mrwin_balance_smd <- function(strata, covariates, weights) {
  out <- numeric(ncol(covariates))
  total_weight <- sum(weights)
  for (j in seq_len(ncol(covariates))) {
    x <- covariates[, j]
    overall <- sum(weights * x) / total_weight
    sd_all <- sqrt(sum(weights * (x - overall)^2) / total_weight)
    if (!is.finite(sd_all) || sd_all <= 0) {
      out[j] <- NA_real_
      next
    }
    smd <- numeric(max(strata))
    for (d in seq_len(max(strata))) {
      idx <- strata == d
      if (!any(idx) || sum(weights[idx]) <= 0) {
        smd[d] <- NA_real_
      } else {
        smd[d] <- (sum(weights[idx] * x[idx]) / sum(weights[idx]) - overall) / sd_all
      }
    }
    out[j] <- max(abs(smd), na.rm = TRUE)
  }
  names(out) <- colnames(covariates)
  out
}

.mrwin_make_contrast_plan <- function(active_strata, n_strata = max(active_strata)) {
  active_strata <- sort(unique(as.integer(active_strata)))
  if (length(active_strata) < 2L) {
    stop("At least two active strata are required.", call. = FALSE)
  }
  low <- active_strata[-length(active_strata)]
  high <- active_strata[-1L]
  data.frame(
    high = high,
    low = low,
    label = paste0("s", high, "_vs_s", low),
    bridged = high - low > 1L,
    stringsAsFactors = FALSE
  )
}

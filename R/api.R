mrwin_endpoint <- function(time, status, priority = NULL) {
  if (missing(time) || missing(status)) {
    stop("`time` and `status` are required.", call. = FALSE)
  }
  if (is.character(time) != is.character(status)) {
    stop("`time` and `status` must both be column names or both be matrices.", call. = FALSE)
  }

  n_priority <- if (is.character(time)) {
    length(time)
  } else {
    ncol(as.matrix(time))
  }
  if (length(status) != length(time) && is.character(time)) {
    stop("`time` and `status` column vectors must have the same length.", call. = FALSE)
  }
  if (!is.character(time) && !all(dim(as.matrix(time)) == dim(as.matrix(status)))) {
    stop("`time` and `status` matrices must have the same dimensions.", call. = FALSE)
  }

  if (is.null(priority)) {
    priority <- if (!is.character(time) && !is.null(colnames(as.matrix(time)))) {
      colnames(as.matrix(time))
    } else {
      paste0("priority_", seq_len(n_priority))
    }
  }
  if (length(priority) != n_priority) {
    stop("`priority` must match the number of endpoint columns.", call. = FALSE)
  }

  structure(
    list(time = time, status = status, priority = as.character(priority)),
    class = "mrwin_endpoint_spec"
  )
}

mrwin_gwas <- function(beta, se = NULL, snp = NULL, covariance = NULL) {
  if (missing(beta)) {
    stop("`beta` is required.", call. = FALSE)
  }
  beta <- as.numeric(beta)
  if (any(!is.finite(beta))) {
    stop("`beta` must contain finite values.", call. = FALSE)
  }

  if (!is.null(covariance)) {
    covariance <- as.matrix(covariance)
    if (!all(dim(covariance) == c(length(beta), length(beta)))) {
      stop("`covariance` must be an M x M matrix matching `beta`.", call. = FALSE)
    }
    if (is.null(se)) {
      se <- sqrt(pmax(diag(covariance), 0))
    }
  }
  if (is.null(se)) {
    se <- rep(0, length(beta))
  }
  se <- as.numeric(se)
  if (length(se) != length(beta)) {
    stop("`se` must have the same length as `beta`.", call. = FALSE)
  }
  if (any(!is.finite(se)) || any(se < 0)) {
    stop("`se` must be finite and non-negative.", call. = FALSE)
  }
  if (is.null(snp)) {
    snp <- paste0("snp_", seq_along(beta))
  }
  if (length(snp) != length(beta)) {
    stop("`snp` must have the same length as `beta`.", call. = FALSE)
  }

  structure(
    list(beta = beta, se = se, snp = as.character(snp), covariance = covariance),
    class = "mrwin_gwas_spec"
  )
}

mrwin_controls <- function(
    n_strata = 10L,
    bootstrap = 200L,
    seed = NULL,
    block_size = 4000L,
    run_sdpd = TRUE,
    sdpd_scale = c("aalen", "cox", "both", "none"),
    sdpd_alpha = 0.05,
    sdpd_min_snps = 10L,
    pleiotropy_bias_radius = NULL,
    adjustment = c("none", "ordinal_iptw", "gps"),
    backend = c("dense", "sparse", "fast", "rcpp"),
    inference = c("bootstrap", "analytic"),
    delta_x_tol = 1e-8,
    iptw_truncation = c(0.01, 0.99),
    ess_fraction = 0.5
) {
  sdpd_scale <- match.arg(sdpd_scale)
  adjustment <- match.arg(adjustment)
  backend <- match.arg(backend)
  inference <- match.arg(inference)
  if (!is.null(pleiotropy_bias_radius)) {
    pleiotropy_bias_radius <- as.numeric(pleiotropy_bias_radius)
  }

  out <- list(
    n_strata = as.integer(n_strata),
    bootstrap = as.integer(bootstrap),
    seed = seed,
    block_size = as.integer(block_size),
    run_sdpd = isTRUE(run_sdpd),
    sdpd_scale = sdpd_scale,
    sdpd_alpha = sdpd_alpha,
    sdpd_min_snps = as.integer(sdpd_min_snps),
    pleiotropy_bias_radius = pleiotropy_bias_radius,
    adjustment = adjustment,
    backend = backend,
    inference = inference,
    delta_x_tol = delta_x_tol,
    iptw_truncation = iptw_truncation,
    ess_fraction = ess_fraction
  )
  .mrwin_validate_controls(out)
  structure(out, class = "mrwin_controls")
}

mrwin <- function(
    data = NULL,
    endpoint,
    genotype,
    exposure,
    gwas = NULL,
    beta_gwas = NULL,
    se_gwas = NULL,
    covariates = NULL,
    controls = mrwin_controls(),
    n_strata = NULL,
    bootstrap = NULL,
    seed = NULL,
    run_sdpd = NULL,
    ...
) {
  call <- match.call()
  dots <- list(...)
  if (length(dots) > 0L) {
    stop("Unused arguments: ", paste(names(dots), collapse = ", "), call. = FALSE)
  }

  controls <- .mrwin_as_controls(controls)
  if (!is.null(n_strata)) controls$n_strata <- as.integer(n_strata)
  if (!is.null(bootstrap)) controls$bootstrap <- as.integer(bootstrap)
  if (!is.null(seed)) controls$seed <- seed
  if (!is.null(run_sdpd)) controls$run_sdpd <- isTRUE(run_sdpd)
  .mrwin_validate_controls(controls)

  warning_log <- list()
  if (controls$backend == "rcpp") {
    stop("`backend = \"rcpp\"` is planned for a later work package; use `dense` or `sparse`.", call. = FALSE)
  }
  if (controls$adjustment == "gps") {
    stop("`adjustment = \"gps\"` is planned for a later work package; use `ordinal_iptw` or `none`.", call. = FALSE)
  }
  if (controls$adjustment != "none" && is.null(covariates)) {
    stop("`covariates` are required when adjustment is not 'none'.", call. = FALSE)
  }
  if (controls$adjustment == "none" && !is.null(covariates)) {
    warning_log <- .mrwin_add_warning(
      warning_log,
      "covariates_ignored",
      "Covariates were supplied, but adjustment = 'none', so they were not used."
    )
  }

  validated <- mrwin_validate_data(
    data = data,
    endpoint = endpoint,
    genotype = genotype,
    exposure = exposure,
    gwas = gwas,
    beta_gwas = beta_gwas,
    se_gwas = se_gwas,
    covariates = covariates,
    controls = controls
  )
  endpoint <- validated$endpoint
  endpoint_data <- list(time = validated$time, status = validated$status)
  G <- validated$G
  X <- validated$X
  gwas <- validated$gwas
  warning_log <- c(warning_log, .mrwin_issue_warnings_to_fit_warnings(validated$issues))

  if (!is.null(gwas$covariance)) {
    warning_log <- .mrwin_add_warning(
      warning_log,
      "full_covariance_not_used",
      "Full GWAS covariance was supplied but WP1 bootstrap uses the diagonal standard errors."
    )
  }

  inference <- if (is.null(controls$inference)) "bootstrap" else controls$inference
  if (inference == "analytic") {
    if (controls$adjustment != "none") {
      stop("`inference = \"analytic\"` currently supports `adjustment = \"none\"` only; ",
           "use `inference = \"bootstrap\"` with IPTW adjustment.", call. = FALSE)
    }
    boot <- mrwin_analytic_bootstrap(
      time = endpoint_data$time,
      status = endpoint_data$status,
      G = G,
      X = X,
      beta_hat = gwas$beta,
      sigma_beta = gwas$se,
      n_strata = controls$n_strata,
      B_gwas = controls$bootstrap,
      seed = controls$seed,
      block_size = controls$block_size
    )
  } else if (controls$backend %in% c("sparse", "fast")) {
    boot <- mrwin_sparse_bootstrap(
      time = endpoint_data$time,
      status = endpoint_data$status,
      G = G,
      X = X,
      beta_hat = gwas$beta,
      sigma_beta = gwas$se,
      n_strata = controls$n_strata,
      B = controls$bootstrap,
      seed = controls$seed,
      covariates = validated$covariates,
      adjustment = controls$adjustment,
      iptw_truncation = controls$iptw_truncation,
      ess_fraction = controls$ess_fraction,
      fast = identical(controls$backend, "fast")
    )
  } else {
    boot <- mrwin_multiplier_bootstrap(
      time = endpoint_data$time,
      status = endpoint_data$status,
      G = G,
      X = X,
      beta_hat = gwas$beta,
      sigma_beta = gwas$se,
      n_strata = controls$n_strata,
      B = controls$bootstrap,
      seed = controls$seed,
      block_size = controls$block_size,
      covariates = validated$covariates,
      adjustment = controls$adjustment,
      iptw_truncation = controls$iptw_truncation,
      ess_fraction = controls$ess_fraction
    )
  }

  if (any(abs(boot$point_delta_x) <= controls$delta_x_tol, na.rm = TRUE) ||
      any(!is.finite(boot$ci95_delta_fieller))) {
    warning_log <- .mrwin_add_warning(
      warning_log,
      "weak_instrument",
      "At least one adjacent phenotypic shift is near zero or the Fieller interval is unbounded."
    )
  }
  if (isTRUE(boot$adjustment$has_positivity_failure)) {
    warning_log <- .mrwin_add_warning(
      warning_log,
      "positivity_failure",
      "At least one PRS stratum failed the IPTW effective-sample-size threshold."
    )
  }
  if (isTRUE(boot$adjustment$has_bridging)) {
    warning_log <- .mrwin_add_warning(
      warning_log,
      "bridged_strata",
      "Adjacent contrasts were bridged across one or more positivity-filtered strata."
    )
  }

  pleiotropy_bounded <- NULL
  if (!is.null(controls$pleiotropy_bias_radius)) {
    pleiotropy_bounded <- mrwin_pleiotropy_bounded_ci(
      delta_hat = boot$delta_gls,
      se_delta = boot$se_delta_gls,
      bias_radius = controls$pleiotropy_bias_radius
    )
  }

  sdpd <- NULL
  if (controls$run_sdpd && controls$sdpd_scale != "none") {
    sdpd <- mrwin_sdpd(
      G = G,
      X = X,
      time = endpoint_data$time[, 1L],
      status = endpoint_data$status[, 1L],
      scale = controls$sdpd_scale,
      alpha = controls$sdpd_alpha,
      min_snps = controls$sdpd_min_snps
    )
    if (isTRUE(sdpd$rejected)) {
      warning_log <- .mrwin_add_warning(
        warning_log,
        "sdpd_rejected",
        "SDPD MR-Egger intercept rejected on at least one scale; cCWR validity is questionable."
      )
    }
    if (isTRUE(sdpd$underpowered)) {
      warning_log <- .mrwin_add_warning(
        warning_log,
        "sdpd_underpowered",
        "SDPD used fewer SNPs than the configured minimum; non-rejection may have low power."
      )
    }
  }

  fit <- list(
    call = call,
    data_info = list(
      n = nrow(G),
      m_snps = ncol(G),
      n_priorities = ncol(endpoint_data$time),
      n_strata = controls$n_strata,
      bootstrap = controls$bootstrap,
      rows_used = validated$rows_used
    ),
    endpoint_info = list(priority = endpoint$priority),
    instrument_info = list(
      snp = gwas$snp,
      beta = gwas$beta,
      se = gwas$se,
      diagonal_gwas_covariance = is.null(gwas$covariance)
    ),
    point = list(
      log_theta = boot$point_log_theta,
      delta_x = boot$point_delta_x,
      delta_isg = boot$delta_isg,
      delta_gls = boot$delta_gls,
      dscwr = boot$dscwr
    ),
    inference = list(
      se_delta_gls = boot$se_delta_gls,
      ci95_delta = boot$ci95_delta,
      ci95_dscwr = boot$ci95_dscwr,
      ci95_delta_fieller = boot$ci95_delta_fieller,
      ci95_dscwr_fieller = boot$ci95_dscwr_fieller,
      pleiotropy_bounded = pleiotropy_bounded
    ),
    heterogeneity = list(
      q = boot$q,
      q_df = boot$q_df,
      q_p_value = boot$q_p_value
    ),
    diagnostics = list(
      validation = validated$issues,
      ledoit_wolf_rho = boot$ledoit_wolf_rho,
      kurtosis_log_theta = boot$kurtosis_log_theta,
      kurtosis_delta_x = boot$kurtosis_delta_x,
      skew_log_theta = boot$skew_log_theta,
      skew_delta_x = boot$skew_delta_x,
      n_valid_bootstrap = boot$n_valid,
      ess = boot$adjustment$ess,
      weights = boot$adjustment$weights,
      dropped_strata = boot$dropped_strata,
      active_strata = boot$active_strata,
      balance = boot$adjustment$balance
    ),
    sdpd = sdpd,
    warnings = warning_log,
    controls = controls,
    bootstrap = boot,
    session_info = list(r_version = R.version.string)
  )
  class(fit) <- c("mrwin_fit", "list")
  fit
}

.mrwin_as_controls <- function(controls) {
  if (inherits(controls, "mrwin_controls")) {
    return(controls)
  }
  if (is.list(controls)) {
    return(do.call(mrwin_controls, controls))
  }
  stop("`controls` must be an `mrwin_controls()` object or a list.", call. = FALSE)
}

.mrwin_validate_controls <- function(controls) {
  if (!is.numeric(controls$n_strata) || controls$n_strata < 2L) {
    stop("`n_strata` must be at least 2.", call. = FALSE)
  }
  if (!is.numeric(controls$bootstrap) || controls$bootstrap < 2L) {
    stop("`bootstrap` must be at least 2.", call. = FALSE)
  }
  if (!is.numeric(controls$block_size) || controls$block_size < 1L) {
    stop("`block_size` must be positive.", call. = FALSE)
  }
  if (!is.numeric(controls$sdpd_alpha) || controls$sdpd_alpha <= 0 || controls$sdpd_alpha >= 1) {
    stop("`sdpd_alpha` must be between 0 and 1.", call. = FALSE)
  }
  if (!is.numeric(controls$sdpd_min_snps) || controls$sdpd_min_snps < 3L) {
    stop("`sdpd_min_snps` must be at least 3.", call. = FALSE)
  }
  if (!is.null(controls$pleiotropy_bias_radius)) {
    radius <- as.numeric(controls$pleiotropy_bias_radius)
    if (length(radius) != 1L || !is.finite(radius) || radius < 0) {
      stop("`pleiotropy_bias_radius` must be NULL or a single non-negative finite number.", call. = FALSE)
    }
  }
  if (!is.numeric(controls$delta_x_tol) || controls$delta_x_tol < 0) {
    stop("`delta_x_tol` must be non-negative.", call. = FALSE)
  }
  .mrwin_validate_truncation(controls$iptw_truncation)
  .mrwin_validate_ess_fraction(controls$ess_fraction)
  invisible(TRUE)
}

.mrwin_resolve_endpoint <- function(endpoint, data) {
  time <- .mrwin_resolve_matrix(endpoint$time, data, "endpoint time")
  status <- .mrwin_resolve_matrix(endpoint$status, data, "endpoint status")
  storage.mode(status) <- "integer"
  list(time = time, status = status)
}

.mrwin_resolve_matrix <- function(x, data, label) {
  if (is.character(x)) {
    if (is.null(data)) {
      stop("`data` is required when ", label, " is specified by column names.", call. = FALSE)
    }
    missing_cols <- setdiff(x, names(data))
    if (length(missing_cols) > 0L) {
      stop("Missing ", label, " columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
    }
    out <- as.matrix(data[, x, drop = FALSE])
  } else {
    out <- as.matrix(x)
  }
  storage.mode(out) <- "double"
  out
}

.mrwin_resolve_vector <- function(x, data, label) {
  if (is.character(x) && length(x) == 1L) {
    if (is.null(data)) {
      stop("`data` is required when ", label, " is specified by a column name.", call. = FALSE)
    }
    if (!x %in% names(data)) {
      stop("Missing ", label, " column: ", x, call. = FALSE)
    }
    out <- data[[x]]
  } else {
    out <- x
  }
  as.numeric(out)
}

.mrwin_resolve_gwas <- function(gwas, beta_gwas, se_gwas) {
  if (inherits(gwas, "mrwin_gwas_spec")) {
    return(gwas)
  }
  if (is.list(gwas) && !is.null(gwas$beta)) {
    return(mrwin_gwas(gwas$beta, gwas$se, gwas$snp, gwas$covariance))
  }
  if (is.null(beta_gwas)) {
    stop("Supply `gwas = mrwin_gwas(...)` or `beta_gwas`.", call. = FALSE)
  }
  mrwin_gwas(beta_gwas, se_gwas)
}

.mrwin_lm_per_snp <- function(G, X) {
  G <- as.matrix(G)
  X <- as.numeric(X)
  m <- ncol(G)
  beta <- se <- rep(NA_real_, m)
  for (j in seq_len(m)) {
    design <- cbind(1, G[, j])
    xtx <- crossprod(design)
    inv <- tryCatch(solve(xtx), error = function(e) NULL)
    if (is.null(inv)) {
      beta[j] <- 0
      se[j] <- Inf
      next
    }
    coef <- inv %*% crossprod(design, X)
    resid <- X - drop(design %*% coef)
    sigma2 <- sum(resid^2) / max(length(X) - 2L, 1L)
    vcov <- sigma2 * inv
    beta[j] <- coef[2L, 1L]
    se[j] <- sqrt(max(vcov[2L, 2L], 1e-20))
  }
  list(beta = beta, se = se)
}

.mrwin_run_sdpd <- function(G, X, time, status, scale, alpha) {
  mrwin_sdpd(G = G, X = X, time = time, status = status, scale = scale, alpha = alpha)
}

.mrwin_add_warning <- function(warnings, code, message) {
  warnings[[length(warnings) + 1L]] <- list(code = code, message = message)
  warnings
}

.mrwin_warnings_df <- function(warnings) {
  if (length(warnings) == 0L) {
    return(data.frame(code = character(), message = character()))
  }
  data.frame(
    code = vapply(warnings, `[[`, character(1), "code"),
    message = vapply(warnings, `[[`, character(1), "message"),
    stringsAsFactors = FALSE
  )
}

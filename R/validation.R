mrwin_validate_data <- function(
    data = NULL,
    endpoint,
    genotype,
    exposure,
    gwas = NULL,
    beta_gwas = NULL,
    se_gwas = NULL,
    covariates = NULL,
    controls = mrwin_controls(),
    complete_cases = c("error", "drop"),
    allow_time_zero = TRUE,
    terminal_priority = 1L,
    terminal_tolerance = 0,
    check_terminal_order = TRUE,
    check_strata = TRUE
) {
  complete_cases <- match.arg(complete_cases)
  controls <- .mrwin_as_controls(controls)
  .mrwin_validate_controls(controls)

  endpoint <- if (inherits(endpoint, "mrwin_endpoint_spec")) {
    endpoint
  } else if (is.list(endpoint) && !is.null(endpoint$time) && !is.null(endpoint$status)) {
    mrwin_endpoint(endpoint$time, endpoint$status, endpoint$priority)
  } else {
    stop("`endpoint` must be an `mrwin_endpoint()` object or a list with `time` and `status`.", call. = FALSE)
  }

  time <- .mrwin_resolve_matrix(endpoint$time, data, "endpoint time")
  status <- .mrwin_resolve_matrix(endpoint$status, data, "endpoint status")
  G <- .mrwin_resolve_matrix(genotype, data, "genotype")
  X <- .mrwin_resolve_vector(exposure, data, "exposure")
  gwas <- .mrwin_resolve_gwas(gwas, beta_gwas, se_gwas)
  Z <- if (is.null(covariates)) NULL else .mrwin_resolve_matrix(covariates, data, "covariates")

  issues <- list(errors = list(), warnings = list())
  n <- nrow(G)

  if (!all(dim(time) == dim(status))) {
    issues <- .mrwin_issue(issues, "error", "endpoint_dim_mismatch", "`time` and `status` must have the same dimensions.")
  }
  if (nrow(time) != n || length(X) != n || (!is.null(Z) && nrow(Z) != n)) {
    issues <- .mrwin_issue(issues, "error", "row_count_mismatch", "Endpoint, genotype, exposure, and covariates must have the same row count.")
  }
  if (length(unique(endpoint$priority)) != length(endpoint$priority)) {
    issues <- .mrwin_issue(issues, "error", "duplicate_priority_labels", "`priority` labels must be unique.")
  }
  if (ncol(G) != length(gwas$beta)) {
    issues <- .mrwin_issue(issues, "error", "gwas_beta_mismatch", "Number of genotype columns must match GWAS beta length.")
  }
  if (length(gwas$se) != length(gwas$beta)) {
    issues <- .mrwin_issue(issues, "error", "gwas_se_mismatch", "GWAS standard errors must match beta length.")
  }

  same_rows <- nrow(time) == n && nrow(status) == n && length(X) == n &&
    (is.null(Z) || nrow(Z) == n)
  row_ok <- rep(TRUE, n)
  if (n > 0L && same_rows) {
    row_ok <- row_ok & stats::complete.cases(time) & stats::complete.cases(status) &
      stats::complete.cases(G) & is.finite(X)
    if (!is.null(Z)) {
      row_ok <- row_ok & stats::complete.cases(Z)
    }
  }
  if (any(!row_ok)) {
    n_bad <- sum(!row_ok)
    if (complete_cases == "drop") {
      time <- time[row_ok, , drop = FALSE]
      status <- status[row_ok, , drop = FALSE]
      G <- G[row_ok, , drop = FALSE]
      X <- X[row_ok]
      if (!is.null(Z)) Z <- Z[row_ok, , drop = FALSE]
      issues <- .mrwin_issue(
        issues,
        "warning",
        "rows_dropped",
        paste0(n_bad, " row(s) with missing/non-finite values were dropped.")
      )
    } else {
      issues <- .mrwin_issue(
        issues,
        "error",
        "missing_or_nonfinite",
        paste0(n_bad, " row(s) contain missing/non-finite endpoint, genotype, covariate, or exposure values.")
      )
    }
  }

  status_non_binary <- length(status) > 0L && !all(status %in% c(0, 1))
  if (status_non_binary) {
    issues <- .mrwin_issue(issues, "error", "status_not_binary", "Endpoint status values must be binary 0/1.")
  }
  if (length(time) > 0L && any(!is.finite(time))) {
    issues <- .mrwin_issue(issues, "error", "time_not_finite", "Endpoint times must be finite.")
  }
  if (length(time) > 0L && any(time < 0)) {
    issues <- .mrwin_issue(issues, "error", "time_negative", "Endpoint times must be non-negative.")
  }
  if (!allow_time_zero && length(time) > 0L && any(time == 0)) {
    issues <- .mrwin_issue(issues, "error", "time_zero", "Endpoint times equal to zero are not allowed with `allow_time_zero = FALSE`.")
  }
  if (any(!is.finite(G))) {
    issues <- .mrwin_issue(issues, "error", "genotype_not_finite", "Genotype matrix must contain finite values.")
  }
  if (any(!is.finite(X))) {
    issues <- .mrwin_issue(issues, "error", "exposure_not_finite", "Exposure vector must contain finite values.")
  }
  if (!is.null(Z) && any(!is.finite(Z))) {
    issues <- .mrwin_issue(issues, "error", "covariates_not_finite", "Covariates must contain finite values.")
  }
  if (ncol(G) > 0L) {
    zero_var <- which(apply(G, 2L, stats::sd) <= 0)
    if (length(zero_var) > 0L) {
      issues <- .mrwin_issue(
        issues,
        "error",
        "zero_variance_genotype",
        paste0("Genotype columns with zero variance: ", paste(zero_var, collapse = ", "), ".")
      )
    }
  }
  if (length(gwas$se) > 0L && all(gwas$se == 0)) {
    issues <- .mrwin_issue(
      issues,
      "warning",
      "no_gwas_uncertainty",
      "All GWAS standard errors are zero; bootstrap will not propagate external GWAS uncertainty."
    )
  }
  if (!is.null(gwas$covariance)) {
    eig <- tryCatch(eigen(0.5 * (gwas$covariance + t(gwas$covariance)), symmetric = TRUE, only.values = TRUE)$values,
                    error = function(e) NA_real_)
    if (any(!is.finite(eig)) || min(eig, na.rm = TRUE) < -1e-8) {
      issues <- .mrwin_issue(issues, "error", "gwas_covariance_not_psd", "GWAS covariance must be symmetric positive semi-definite.")
    }
  }

  if (check_terminal_order && same_rows && all(dim(time) == dim(status)) &&
      ncol(time) > 1L && terminal_priority >= 1L && terminal_priority <= ncol(time)) {
    terminal_event <- status[, terminal_priority] == 1
    for (k in setdiff(seq_len(ncol(time)), terminal_priority)) {
      bad <- terminal_event & status[, k] == 1 & time[, k] > time[, terminal_priority] + terminal_tolerance
      if (any(bad)) {
        issues <- .mrwin_issue(
          issues,
          "error",
          "terminal_event_order_violation",
          paste0("Priority column ", k, " has observed events after the terminal priority in ", sum(bad), " row(s).")
        )
      }
    }
  }

  strata <- NULL
  if (check_strata && .mrwin_has_no_errors(issues) && nrow(G) > 0L && ncol(G) == length(gwas$beta)) {
    if (controls$n_strata > nrow(G)) {
      issues <- .mrwin_issue(issues, "error", "too_many_strata", "`n_strata` cannot exceed the number of analysis rows.")
    } else {
      strata <- mrwin_prs_strata(G, gwas$beta, controls$n_strata)$strata
      counts <- tabulate(strata, nbins = controls$n_strata)
      empty <- which(counts == 0L)
      if (length(empty) > 0L) {
        issues <- .mrwin_issue(
          issues,
          "error",
          "empty_stratum",
          paste0("PRS strata with no observations: ", paste(empty, collapse = ", "), ".")
        )
      }
      small <- which(counts < 2L)
      if (length(small) > 0L) {
        issues <- .mrwin_issue(
          issues,
          "warning",
          "small_stratum",
          paste0("PRS strata with fewer than two observations: ", paste(small, collapse = ", "), ".")
        )
      }
    }
  }

  if (!.mrwin_has_no_errors(issues)) {
    .mrwin_stop_issues(issues)
  }

  out <- list(
    time = time,
    status = matrix(as.integer(status), nrow = nrow(status), ncol = ncol(status), dimnames = dimnames(status)),
    G = G,
    X = X,
    covariates = Z,
    gwas = gwas,
    endpoint = endpoint,
    controls = controls,
    strata = strata,
    rows_used = if (exists("row_ok")) which(row_ok) else seq_len(nrow(G)),
    issues = issues
  )
  class(out) <- c("mrwin_validated_data", "list")
  out
}

print.mrwin_validated_data <- function(x, ...) {
  cat("mrwin validated data\n")
  cat("  N:", nrow(x$G), "\n")
  cat("  SNPs:", ncol(x$G), "\n")
  cat("  Priorities:", ncol(x$time), "\n")
  cat("  Strata:", x$controls$n_strata, "\n")
  if (length(x$issues$warnings) > 0L) {
    cat("  Warnings:", length(x$issues$warnings), "\n")
  }
  invisible(x)
}

.mrwin_issue <- function(issues, severity, code, message) {
  issues[[paste0(severity, "s")]][[length(issues[[paste0(severity, "s")]]) + 1L]] <-
    list(code = code, message = message)
  issues
}

.mrwin_has_no_errors <- function(issues) {
  length(issues$errors) == 0L
}

.mrwin_stop_issues <- function(issues) {
  messages <- vapply(
    issues$errors,
    function(x) paste0("[", x$code, "] ", x$message),
    character(1)
  )
  stop(paste(messages, collapse = "\n"), call. = FALSE)
}

.mrwin_issue_warnings_to_fit_warnings <- function(issues) {
  lapply(issues$warnings, function(x) list(code = x$code, message = x$message))
}

print.mrwin_fit <- function(x, ...) {
  cat("mrwin fit\n")
  cat("  N:", x$data_info$n, "\n")
  cat("  SNPs:", x$data_info$m_snps, "\n")
  cat("  Priorities:", x$data_info$n_priorities, "\n")
  cat("  Strata:", x$data_info$n_strata, "\n")
  cat("  Bootstrap:", x$data_info$bootstrap, "(", x$diagnostics$n_valid_bootstrap, "valid )\n", sep = "")
  cat("  delta_GLS:", formatC(x$point$delta_gls, digits = 4, format = "f"), "\n")
  cat("  DS-CWR:", formatC(x$point$dscwr, digits = 4, format = "f"), "\n")
  cat(
    "  95% CI:",
    paste(formatC(x$inference$ci95_dscwr, digits = 4, format = "f"), collapse = " to "),
    "\n"
  )
  if (length(x$warnings) > 0L) {
    cat("  Warnings:", length(x$warnings), "\n")
  }
  invisible(x)
}

summary.mrwin_fit <- function(object, ...) {
  estimate <- data.frame(
    term = "DS-CWR",
    delta = object$point$delta_gls,
    se_delta = object$inference$se_delta_gls,
    dscwr = object$point$dscwr,
    ci_low = object$inference$ci95_dscwr[1L],
    ci_high = object$inference$ci95_dscwr[2L],
    pleiotropy_ci_low = .mrwin_bound_or_na(object$inference$pleiotropy_bounded, "ci_dscwr_pleiotropy_bounded", 1L),
    pleiotropy_ci_high = .mrwin_bound_or_na(object$inference$pleiotropy_bounded, "ci_dscwr_pleiotropy_bounded", 2L),
    stringsAsFactors = FALSE
  )
  heterogeneity <- data.frame(
    q = object$heterogeneity$q,
    df = object$heterogeneity$q_df,
    p_value = object$heterogeneity$q_p_value
  )
  sdpd <- .mrwin_sdpd_df(object$sdpd)
  out <- list(
    call = object$call,
    estimate = estimate,
    adjacent = data.frame(
      contrast = names(object$point$delta_isg),
      log_theta = as.numeric(object$point$log_theta),
      delta_x = as.numeric(object$point$delta_x),
      delta_isg = as.numeric(object$point$delta_isg),
      stringsAsFactors = FALSE
    ),
    heterogeneity = heterogeneity,
    sdpd = sdpd,
    warnings = .mrwin_warnings_df(object$warnings)
  )
  class(out) <- c("summary.mrwin_fit", "list")
  out
}

.mrwin_bound_or_na <- function(x, field, index) {
  if (is.null(x) || is.null(x[[field]])) {
    return(NA_real_)
  }
  x[[field]][index]
}

print.summary.mrwin_fit <- function(x, ...) {
  cat("mrwin summary\n\n")
  print(x$estimate, row.names = FALSE)
  cat("\nHeterogeneity\n")
  print(x$heterogeneity, row.names = FALSE)
  if (nrow(x$sdpd) > 0L) {
    cat("\nSDPD\n")
    print(x$sdpd, row.names = FALSE)
  }
  if (nrow(x$warnings) > 0L) {
    cat("\nWarnings\n")
    print(x$warnings, row.names = FALSE)
  }
  invisible(x)
}

plot.mrwin_fit <- function(x, type = c("isg"), ...) {
  type <- match.arg(type)
  y <- as.numeric(x$point$delta_isg)
  labels <- names(x$point$delta_isg)
  graphics::plot(
    seq_along(y),
    y,
    type = "b",
    xaxt = "n",
    xlab = "Adjacent PRS-stratum contrast",
    ylab = "Instrument-standardized gradient",
    main = "mrwin adjacent gradients",
    ...
  )
  graphics::axis(1, at = seq_along(y), labels = labels, las = 2)
  graphics::abline(h = x$point$delta_gls, lty = 2)
  invisible(x)
}

.mrwin_sdpd_df <- function(sdpd) {
  if (is.null(sdpd) || length(sdpd$results) == 0L) {
    return(data.frame())
  }
  rows <- lapply(sdpd$results, function(x) {
    if (!is.null(x$error)) {
      return(data.frame(
        scale = x$scale,
        intercept = NA_real_,
        intercept_se = NA_real_,
        p_value = NA_real_,
        rejected = NA,
        error = x$error
      ))
    }
    data.frame(
      scale = x$scale,
      intercept = x$intercept,
      intercept_se = x$intercept_se,
      p_value = x$intercept_p_value,
      rejected = x$rejected,
      error = NA_character_
    )
  })
  do.call(rbind, rows)
}

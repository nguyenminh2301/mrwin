print.mrwin_fit <- function(x, ...) {
  cat("mrwin fit\n")
  cat("  N:", x$data_info$n, "\n")
  cat("  SNPs:", x$data_info$m_snps, "\n")
  cat("  Priorities:", x$data_info$n_priorities, "\n")
  cat("  Strata:", x$data_info$n_strata, "\n")
  cat("  Bootstrap:", x$data_info$bootstrap, "(", x$diagnostics$n_valid_bootstrap, "valid )\n", sep = "")
  cat("  Backend:", x$controls$backend, "\n")
  cat("  Adjustment:", x$controls$adjustment, "\n")
  cat("  delta_GLS:", formatC(x$point$delta_gls, digits = 4, format = "f"), "\n")
  cat("  DS-CWR:", formatC(x$point$dscwr, digits = 4, format = "f"), "\n")
  cat(
    "  95% CI:",
    paste(formatC(x$inference$ci95_dscwr, digits = 4, format = "f"), collapse = " to "),
    "\n"
  )
  cat("  Q:", formatC(x$heterogeneity$q, digits = 4, format = "f"),
      "( df =", x$heterogeneity$q_df, ", p =", formatC(x$heterogeneity$q_p_value, digits = 4, format = "f"), ")\n")
  if (length(x$warnings) > 0L) {
    codes <- vapply(x$warnings, `[[`, character(1), "code")
    cat("  Warnings:", length(x$warnings), "(", paste(codes, collapse = ", "), ")\n")
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
    p_value = 2 * stats::pnorm(-abs(object$point$delta_gls / max(object$inference$se_delta_gls, 1e-12))),
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
  adjacent <- data.frame(
    contrast = names(object$point$delta_isg),
    log_theta = as.numeric(object$point$log_theta),
    delta_x = as.numeric(object$point$delta_x),
    delta_isg = as.numeric(object$point$delta_isg),
    cwr = exp(as.numeric(object$point$log_theta)),
    stringsAsFactors = FALSE
  )
  fieller <- data.frame(
    ci_low = object$inference$ci95_delta_fieller[1L],
    ci_high = object$inference$ci95_delta_fieller[2L],
    unbounded = isTRUE(object$bootstrap$fieller_unbounded),
    stringsAsFactors = FALSE
  )
  diagnostics <- data.frame(
    n_valid_bootstrap = object$diagnostics$n_valid_bootstrap,
    ledoit_wolf_rho = object$diagnostics$ledoit_wolf_rho,
    n_strata = object$data_info$n_strata,
    n_contrasts = length(object$point$delta_isg),
    stringsAsFactors = FALSE
  )
  out <- list(
    call = object$call,
    estimate = estimate,
    adjacent = adjacent,
    heterogeneity = heterogeneity,
    fieller = fieller,
    sdpd = sdpd,
    diagnostics = diagnostics,
    warnings = .mrwin_warnings_df(object$warnings),
    data_info = object$data_info,
    controls = object$controls,
    endpoint_info = object$endpoint_info
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
  cat("mrwin summary\n")
  cat(strrep("-", 60), "\n")
  cat("Data:", x$data_info$n, "observations,", x$data_info$m_snps, "SNPs,",
      x$data_info$n_priorities, "priorities\n")
  cat("Strata:", x$data_info$n_strata, " | Bootstrap:", x$data_info$bootstrap,
      "(", x$diagnostics$n_valid_bootstrap, "valid) | Backend:", x$controls$backend, "\n")
  cat(strrep("-", 60), "\n\n")

  cat("Main Estimate\n")
  est <- x$estimate
  cat(sprintf("  DS-CWR:     %.4f  (95%% CI: %.4f to %.4f)\n", est$dscwr, est$ci_low, est$ci_high))
  cat(sprintf("  delta_GLS:  %.4f  (SE: %.4f, p = %.4f)\n", est$delta, est$se_delta, est$p_value))
  if (!is.na(est$pleiotropy_ci_low)) {
    cat(sprintf("  Pleiotropy-bounded CI: %.4f to %.4f\n", est$pleiotropy_ci_low, est$pleiotropy_ci_high))
  }

  cat("\nAdjacent Stratum Gradients\n")
  adj <- x$adjacent
  for (i in seq_len(nrow(adj))) {
    cat(sprintf("  %-12s  log(theta)=%.4f  Delta_X=%.4f  ISG=%.4f  CWR=%.4f\n",
                adj$contrast[i], adj$log_theta[i], adj$delta_x[i], adj$delta_isg[i], adj$cwr[i]))
  }

  cat("\nHeterogeneity\n")
  het <- x$heterogeneity
  cat(sprintf("  Q = %.4f  (df = %d, p = %.4f)\n", het$q, het$df, het$p_value))

  cat("\nFieller CI\n")
  cat(sprintf("  delta: %.4f to %.4f  (unbounded: %s)\n",
              x$fieller$ci_low, x$fieller$ci_high, x$fieller$unbounded))

  if (nrow(x$sdpd) > 0L) {
    cat("\nSDPD Diagnostics\n")
    for (i in seq_len(nrow(x$sdpd))) {
      s <- x$sdpd[i, ]
      if (is.na(s$error)) {
        cat(sprintf("  %s: intercept=%.4f (p=%.4f) rejected=%s\n",
                    s$scale, s$intercept, s$p_value, s$rejected))
      } else {
        cat(sprintf("  %s: error=%s\n", s$scale, s$error))
      }
    }
  }

  if (nrow(x$warnings) > 0L) {
    cat("\nWarnings\n")
    for (i in seq_len(nrow(x$warnings))) {
      cat(sprintf("  [%s] %s\n", x$warnings$code[i], x$warnings$message[i]))
    }
  }
  cat(strrep("-", 60), "\n")
  invisible(x)
}

plot.mrwin_fit <- function(x, type = c("isg", "forest", "bootstrap"), ...) {
  type <- match.arg(type)
  switch(type,
    isg = .mrwin_plot_isg(x, ...),
    forest = .mrwin_plot_forest(x, ...),
    bootstrap = .mrwin_plot_bootstrap(x, ...)
  )
  invisible(x)
}

.mrwin_plot_isg <- function(x, ...) {
  y <- as.numeric(x$point$delta_isg)
  labels <- names(x$point$delta_isg)
  graphics::plot(
    seq_along(y), y,
    type = "b", pch = 19, lwd = 2,
    xaxt = "n", xlab = "Adjacent PRS-stratum contrast",
    ylab = "Instrument-standardized gradient",
    main = "Adjacent ISG Gradients",
    ...
  )
  graphics::axis(1, at = seq_along(y), labels = labels, las = 2)
  graphics::abline(h = x$point$delta_gls, lty = 2, col = "red", lwd = 1.5)
  graphics::abline(h = 0, lty = 3, col = "gray50")
  graphics::legend("topright",
    legend = c("ISG", "GLS pooled", "null"),
    lty = c(1, 2, 3), pch = c(19, NA, NA),
    col = c("black", "red", "gray50"), lwd = c(2, 1.5, 1),
    cex = 0.8
  )
}

.mrwin_plot_forest <- function(x, ...) {
  adj <- summary(x)$adjacent
  n <- nrow(adj)
  if (n == 0L) {
    message("No adjacent contrasts to plot.")
    return(invisible(x))
  }
  y_pos <- seq_len(n)
  graphics::plot(
    adj$cwr, y_pos,
    type = "n", pch = 19,
    xlim = range(c(adj$cwr, 1), na.rm = TRUE),
    ylim = c(0.5, n + 0.5),
    xlab = "CWR (exp(log-theta))",
    ylab = "",
    yaxt = "n",
    main = "Forest Plot: Adjacent CWR",
    ...
  )
  graphics::axis(2, at = y_pos, labels = adj$contrast, las = 2)
  graphics::abline(v = 1, lty = 2, col = "gray50")
  graphics::points(adj$cwr, y_pos, pch = 19, cex = 1.2)
  graphics::segments(
    exp(adj$log_theta - 1.96 * abs(adj$log_theta) / pmax(abs(adj$delta_isg), 1e-12)),
    y_pos,
    exp(adj$log_theta + 1.96 * abs(adj$log_theta) / pmax(abs(adj$delta_isg), 1e-12)),
    y_pos,
    lwd = 1.5
  )
}

.mrwin_plot_bootstrap <- function(x, ...) {
  lt <- x$bootstrap$bootstrap_log_theta
  if (is.null(lt)) {
    message("No bootstrap draws available.")
    return(invisible(x))
  }
  valid <- x$bootstrap$valid_bootstrap
  lt_valid <- lt[valid, , drop = FALSE]
  n_contrasts <- ncol(lt_valid)
  if (n_contrasts == 0L) {
    message("No valid bootstrap contrasts.")
    return(invisible(x))
  }
  graphics::boxplot(
    lt_valid,
    names = colnames(lt_valid),
    main = "Bootstrap Distribution: log-theta",
    ylab = "log(theta)",
    las = 2,
    ...
  )
  graphics::abline(h = 0, lty = 3, col = "gray50")
  points <- as.numeric(x$point$log_theta)
  graphics::points(seq_along(points), points, pch = 19, col = "red", cex = 1.2)
  graphics::legend("topright",
    legend = c("bootstrap draws", "point estimate"),
    pch = c(NA, 19), lty = c(1, NA), col = c("black", "red"),
    cex = 0.8
  )
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

tidy <- function(x, ...) {
  UseMethod("tidy")
}

tidy.mrwin_fit <- function(x, ...) {
  data.frame(
    term = "DS-CWR",
    estimate = x$point$dscwr,
    delta = x$point$delta_gls,
    se_delta = x$inference$se_delta_gls,
    statistic = x$point$delta_gls / max(x$inference$se_delta_gls, 1e-12),
    p_value = 2 * stats::pnorm(-abs(x$point$delta_gls / max(x$inference$se_delta_gls, 1e-12))),
    ci_low = x$inference$ci95_dscwr[1L],
    ci_high = x$inference$ci95_dscwr[2L],
    ci_delta_low = x$inference$ci95_delta[1L],
    ci_delta_high = x$inference$ci95_delta[2L],
    q_stat = x$heterogeneity$q,
    q_p_value = x$heterogeneity$q_p_value,
    n = x$data_info$n,
    m_snps = x$data_info$m_snps,
    n_strata = x$data_info$n_strata,
    n_valid_bootstrap = x$diagnostics$n_valid_bootstrap,
    stringsAsFactors = FALSE
  )
}

mrwin_report <- function(x, file = NULL, format = c("text", "markdown")) {
  if (!inherits(x, "mrwin_fit")) {
    stop("`x` must be an `mrwin_fit` object.", call. = FALSE)
  }
  format <- match.arg(format)
  if (format == "markdown") {
    lines <- .mrwin_report_markdown(x)
  } else {
    lines <- .mrwin_report_text(x)
  }
  text <- paste(lines, collapse = "\n")
  if (!is.null(file)) {
    writeLines(text, file)
    invisible(text)
  } else {
    cat(text)
    invisible(text)
  }
}

.mrwin_report_text <- function(x) {
  s <- summary(x)
  lines <- character()
  lines <- c(lines, "mrwin Analysis Report", strrep("=", 60), "")
  lines <- c(lines, sprintf("Date: %s", Sys.time()))
  lines <- c(lines, sprintf("R version: %s", x$session_info$r_version), "")

  lines <- c(lines, "Data Summary", strrep("-", 40))
  lines <- c(lines, sprintf("  Observations: %d", x$data_info$n))
  lines <- c(lines, sprintf("  SNPs: %d", x$data_info$m_snps))
  lines <- c(lines, sprintf("  Priorities: %d", x$data_info$n_priorities))
  lines <- c(lines, sprintf("  Strata: %d", x$data_info$n_strata))
  lines <- c(lines, sprintf("  Backend: %s", x$controls$backend))
  lines <- c(lines, sprintf("  Adjustment: %s", x$controls$adjustment), "")

  lines <- c(lines, "Main Estimate", strrep("-", 40))
  est <- s$estimate
  lines <- c(lines, sprintf("  DS-CWR: %.4f (95%% CI: %.4f to %.4f)", est$dscwr, est$ci_low, est$ci_high))
  lines <- c(lines, sprintf("  delta_GLS: %.4f (SE: %.4f, p = %.4f)", est$delta, est$se_delta, est$p_value))
  if (!is.na(est$pleiotropy_ci_low)) {
    lines <- c(lines, sprintf("  Pleiotropy-bounded CI: %.4f to %.4f", est$pleiotropy_ci_low, est$pleiotropy_ci_high))
  }
  lines <- c(lines, "")

  lines <- c(lines, "Adjacent Stratum Gradients", strrep("-", 40))
  adj <- s$adjacent
  for (i in seq_len(nrow(adj))) {
    lines <- c(lines, sprintf("  %s: log(theta)=%.4f, Delta_X=%.4f, ISG=%.4f, CWR=%.4f",
                              adj$contrast[i], adj$log_theta[i], adj$delta_x[i], adj$delta_isg[i], adj$cwr[i]))
  }
  lines <- c(lines, "")

  lines <- c(lines, "Heterogeneity", strrep("-", 40))
  het <- s$heterogeneity
  lines <- c(lines, sprintf("  Q = %.4f (df = %d, p = %.4f)", het$q, het$df, het$p_value))
  lines <- c(lines, "")

  lines <- c(lines, "Inference", strrep("-", 40))
  lines <- c(lines, sprintf("  Bootstrap: %d/%d valid", x$diagnostics$n_valid_bootstrap, x$data_info$bootstrap))
  lines <- c(lines, sprintf("  Fieller CI: %.4f to %.4f (unbounded: %s)",
                            s$fieller$ci_low, s$fieller$ci_high, s$fieller$unbounded))
  lines <- c(lines, sprintf("  Ledoit-Wolf rho: %.4f", x$diagnostics$ledoit_wolf_rho), "")

  if (nrow(s$sdpd) > 0L) {
    lines <- c(lines, "SDPD Diagnostics", strrep("-", 40))
    for (i in seq_len(nrow(s$sdpd))) {
      row <- s$sdpd[i, ]
      if (is.na(row$error)) {
        lines <- c(lines, sprintf("  %s: intercept=%.4f, p=%.4f, rejected=%s",
                                  row$scale, row$intercept, row$p_value, row$rejected))
      } else {
        lines <- c(lines, sprintf("  %s: %s", row$scale, row$error))
      }
    }
    lines <- c(lines, "")
  }

  if (nrow(s$warnings) > 0L) {
    lines <- c(lines, "Warnings", strrep("-", 40))
    for (i in seq_len(nrow(s$warnings))) {
      lines <- c(lines, sprintf("  [%s] %s", s$warnings$code[i], s$warnings$message[i]))
    }
    lines <- c(lines, "")
  }

  lines <- c(lines, strrep("=", 60))
  lines
}

.mrwin_report_markdown <- function(x) {
  s <- summary(x)
  lines <- character()
  lines <- c(lines, "# mrwin Analysis Report", "")
  lines <- c(lines, sprintf("**Date:** %s  ", Sys.time()))
  lines <- c(lines, sprintf("**R version:** %s  ", x$session_info$r_version), "")

  lines <- c(lines, "## Data Summary", "")
  lines <- c(lines, "| Parameter | Value |")
  lines <- c(lines, "|---|---|")
  lines <- c(lines, sprintf("| Observations | %d |", x$data_info$n))
  lines <- c(lines, sprintf("| SNPs | %d |", x$data_info$m_snps))
  lines <- c(lines, sprintf("| Priorities | %d |", x$data_info$n_priorities))
  lines <- c(lines, sprintf("| Strata | %d |", x$data_info$n_strata))
  lines <- c(lines, sprintf("| Backend | %s |", x$controls$backend))
  lines <- c(lines, sprintf("| Adjustment | %s |", x$controls$adjustment), "")

  lines <- c(lines, "## Main Estimate", "")
  est <- s$estimate
  lines <- c(lines, "| Metric | Value | 95% CI |")
  lines <- c(lines, "|---|---|---|")
  lines <- c(lines, sprintf("| DS-CWR | %.4f | %.4f to %.4f |", est$dscwr, est$ci_low, est$ci_high))
  lines <- c(lines, sprintf("| delta_GLS | %.4f | %.4f to %.4f |", est$delta, est$ci_delta_low, est$ci_delta_high))
  lines <- c(lines, sprintf("| p-value | %.4f | |", est$p_value), "")

  lines <- c(lines, "## Adjacent Stratum Gradients", "")
  lines <- c(lines, "| Contrast | log(theta) | Delta_X | ISG | CWR |")
  lines <- c(lines, "|---|---|---|---|---|")
  adj <- s$adjacent
  for (i in seq_len(nrow(adj))) {
    lines <- c(lines, sprintf("| %s | %.4f | %.4f | %.4f | %.4f |",
                              adj$contrast[i], adj$log_theta[i], adj$delta_x[i], adj$delta_isg[i], adj$cwr[i]))
  }
  lines <- c(lines, "")

  lines <- c(lines, "## Heterogeneity", "")
  het <- s$heterogeneity
  lines <- c(lines, sprintf("- **Q** = %.4f (df = %d, p = %.4f)", het$q, het$df, het$p_value), "")

  lines <- c(lines, "## Inference Details", "")
  lines <- c(lines, sprintf("- Bootstrap: %d/%d valid iterations", x$diagnostics$n_valid_bootstrap, x$data_info$bootstrap))
  lines <- c(lines, sprintf("- Fieller CI: %.4f to %.4f (unbounded: %s)",
                            s$fieller$ci_low, s$fieller$ci_high, s$fieller$unbounded))
  lines <- c(lines, sprintf("- Ledoit-Wolf rho: %.4f", x$diagnostics$ledoit_wolf_rho), "")

  if (nrow(s$sdpd) > 0L) {
    lines <- c(lines, "## SDPD Diagnostics", "")
    lines <- c(lines, "| Scale | Intercept | p-value | Rejected |")
    lines <- c(lines, "|---|---|---|---|")
    for (i in seq_len(nrow(s$sdpd))) {
      row <- s$sdpd[i, ]
      if (is.na(row$error)) {
        lines <- c(lines, sprintf("| %s | %.4f | %.4f | %s |", row$scale, row$intercept, row$p_value, row$rejected))
      } else {
        lines <- c(lines, sprintf("| %s | %s | | |", row$scale, row$error))
      }
    }
    lines <- c(lines, "")
  }

  if (nrow(s$warnings) > 0L) {
    lines <- c(lines, "## Warnings", "")
    for (i in seq_len(nrow(s$warnings))) {
      lines <- c(lines, sprintf("- **[%s]** %s", s$warnings$code[i], s$warnings$message[i]))
    }
    lines <- c(lines, "")
  }

  lines <- c(lines, "---", sprintf("*Report generated by mrwin v%s*", "0.0.0.9000"))
  lines
}

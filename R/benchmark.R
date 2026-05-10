mrwin_benchmark <- function(
    n_vec = c(200, 500, 1000),
    m_snps = 20L,
    n_strata = 5L,
    B = 50L,
    backends = c("dense", "sparse"),
    seed = 20260510L,
    n_iter = 1L
) {
  if (!is.numeric(n_vec) || any(n_vec < 10L)) {
    stop("`n_vec` must be a numeric vector of sample sizes at least 10.", call. = FALSE)
  }
  backends <- match.arg(backends, several.ok = TRUE)
  n_iter <- as.integer(n_iter)
  if (length(n_iter) != 1L || is.na(n_iter) || n_iter < 1L) {
    stop("`n_iter` must be a positive integer.", call. = FALSE)
  }

  rows <- list()
  row_id <- 1L
  for (n in n_vec) {
    cfg <- mrwin_config(n_outcome = n, m_snps = m_snps, seed = seed)
    dat <- mrwin_simulate(cfg, seed = seed)
    for (iter in seq_len(n_iter)) {
      for (backend in backends) {
        ctrl <- mrwin_controls(
          n_strata = n_strata,
          bootstrap = B,
          seed = seed + iter,
          backend = backend,
          run_sdpd = FALSE
        )
        t_start <- proc.time()
        fit <- tryCatch(
          mrwin(
            endpoint = mrwin_endpoint(dat$time, dat$status, colnames(dat$time)),
            genotype = dat$G,
            exposure = dat$X,
            gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
            controls = ctrl
          ),
          error = function(e) NULL
        )
        t_elapsed <- (proc.time() - t_start)[["elapsed"]]

        rows[[row_id]] <- data.frame(
          n = as.integer(n),
          m_snps = as.integer(m_snps),
          n_strata = as.integer(n_strata),
          B = as.integer(B),
          backend = backend,
          iteration = iter,
          success = !is.null(fit),
          time_seconds = t_elapsed,
          delta_gls = if (!is.null(fit)) fit$point$delta_gls else NA_real_,
          dscwr = if (!is.null(fit)) fit$point$dscwr else NA_real_,
          n_valid_boot = if (!is.null(fit)) fit$diagnostics$n_valid_bootstrap else NA_integer_,
          stringsAsFactors = FALSE
        )
        row_id <- row_id + 1L
      }
    }
  }

  results <- do.call(rbind, rows)
  rownames(results) <- NULL
  structure(
    list(
      results = results,
      n_vec = as.integer(n_vec),
      m_snps = as.integer(m_snps),
      n_strata = as.integer(n_strata),
      B = as.integer(B),
      backends = backends,
      seed = seed,
      n_iter = as.integer(n_iter)
    ),
    class = c("mrwin_benchmark", "list")
  )
}

print.mrwin_benchmark <- function(x, ...) {
  cat("mrwin benchmark\n")
  cat("  sample sizes:", paste(x$n_vec, collapse = ", "), "\n")
  cat("  SNPs:", x$m_snps, "\n")
  cat("  strata:", x$n_strata, "\n")
  cat("  bootstrap:", x$B, "\n")
  cat("  backends:", paste(x$backends, collapse = ", "), "\n")
  cat("  iterations per config:", x$n_iter, "\n")
  cat("  total runs:", nrow(x$results), "\n")
  if (nrow(x$results) > 0L) {
    cat("\n  Median runtime (seconds):\n")
    agg <- stats::aggregate(
      time_seconds ~ n + backend,
      data = x$results,
      FUN = stats::median
    )
    print(agg, row.names = FALSE)
  }
  invisible(x)
}

summary.mrwin_benchmark <- function(object, ...) {
  if (nrow(object$results) == 0L) {
    return(data.frame())
  }
  agg <- stats::aggregate(
    time_seconds ~ n + m_snps + n_strata + B + backend,
    data = object$results,
    FUN = function(x) c(
      median = stats::median(x),
      mean = mean(x),
      min = min(x),
      max = max(x),
      sd = stats::sd(x)
    )
  )
  cbind(
    agg[, 1:5, drop = FALSE],
    as.data.frame(agg$time_seconds)
  )
}

mrwin_verify_sparse_dense_parity <- function(
    n = 200L,
    m_snps = 10L,
    n_strata = 4L,
    B = 20L,
    seed = 20260510L,
    tol = 1e-10
) {
  cfg <- mrwin_config(n_outcome = n, m_snps = m_snps, seed = seed)
  dat <- mrwin_simulate(cfg, seed = seed)

  ctrl_dense <- mrwin_controls(
    n_strata = n_strata,
    bootstrap = B,
    seed = seed,
    backend = "dense",
    run_sdpd = FALSE
  )
  ctrl_sparse <- mrwin_controls(
    n_strata = n_strata,
    bootstrap = B,
    seed = seed,
    backend = "sparse",
    run_sdpd = FALSE
  )

  fit_dense <- mrwin(
    endpoint = mrwin_endpoint(dat$time, dat$status, colnames(dat$time)),
    genotype = dat$G,
    exposure = dat$X,
    gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
    controls = ctrl_dense
  )
  fit_sparse <- mrwin(
    endpoint = mrwin_endpoint(dat$time, dat$status, colnames(dat$time)),
    genotype = dat$G,
    exposure = dat$X,
    gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
    controls = ctrl_sparse
  )

  checks <- list(
    delta_gls = abs(fit_dense$point$delta_gls - fit_sparse$point$delta_gls),
    dscwr = abs(fit_dense$point$dscwr - fit_sparse$point$dscwr),
    se_delta = abs(fit_dense$inference$se_delta_gls - fit_sparse$inference$se_delta_gls),
    q = abs(fit_dense$heterogeneity$q - fit_sparse$heterogeneity$q),
    log_theta_max = max(abs(fit_dense$point$log_theta - fit_sparse$point$log_theta)),
    delta_x_max = max(abs(fit_dense$point$delta_x - fit_sparse$point$delta_x)),
    delta_isg_max = max(abs(fit_dense$point$delta_isg - fit_sparse$point$delta_isg))
  )

  all_close <- all(vapply(checks, function(x) is.finite(x) && x < tol, logical(1)))

  structure(
    list(
      n = as.integer(n),
      m_snps = as.integer(m_snps),
      n_strata = as.integer(n_strata),
      B = as.integer(B),
      seed = seed,
      tol = tol,
      checks = checks,
      all_close = all_close,
      dense = fit_dense,
      sparse = fit_sparse
    ),
    class = c("mrwin_parity_check", "list")
  )
}

print.mrwin_parity_check <- function(x, ...) {
  cat("mrwin sparse/dense parity check\n")
  cat("  N:", x$n, " SNPs:", x$m_snps, " strata:", x$n_strata, " B:", x$B, "\n")
  cat("  tolerance:", formatC(x$tol, format = "e"), "\n")
  cat("  all close:", x$all_close, "\n")
  cat("\n  Max absolute differences:\n")
  for (nm in names(x$checks)) {
    cat(sprintf("    %-20s %s\n", nm, formatC(x$checks[[nm]], format = "e", digits = 4)))
  }
  invisible(x)
}

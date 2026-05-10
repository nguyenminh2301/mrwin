mrwin_scenarios <- function(
    base_config = mrwin_config(),
    scenarios = c("A_null", "B_valid_IV", "C_pleiotropy", "D_hierarchy_discordant")
) {
  if (!inherits(base_config, "mrwin_config")) {
    base_config <- do.call(mrwin_config, base_config)
  }
  scenarios <- match.arg(scenarios, several.ok = TRUE)
  out <- lapply(scenarios, function(name) {
    cfg <- base_config
    meta <- switch(
      name,
      A_null = list(
        alpha_x = c(0, 0, 0),
        gamma_direct = c(0, 0, 0),
        true_effect = 0,
        expected_direction = "null",
        package = "A",
        description = "Null calibration scenario."
      ),
      B_valid_IV = list(
        alpha_x = c(-0.4, -0.4, -0.4),
        gamma_direct = c(0, 0, 0),
        true_effect = NA_real_,
        expected_direction = "protective",
        package = "A",
        description = "Concordant valid-IV causal effect scenario."
      ),
      C_pleiotropy = list(
        alpha_x = c(-0.4, -0.4, -0.4),
        gamma_direct = c(0.05, 0, 0),
        true_effect = NA_real_,
        expected_direction = "pleiotropy",
        package = "A",
        description = "Priority-1 direct pleiotropy scenario."
      ),
      D_hierarchy_discordant = list(
        alpha_x = c(0.4, 0, -0.4),
        gamma_direct = c(0, 0, 0),
        true_effect = NA_real_,
        expected_direction = "mixed",
        package = "D",
        description = "Discordant component effects across endpoint priorities."
      )
    )
    cfg$alpha_x <- meta$alpha_x
    cfg$gamma_direct <- meta$gamma_direct
    list(
      name = name,
      package = meta$package,
      description = meta$description,
      expected_direction = meta$expected_direction,
      true_effect = meta$true_effect,
      config = cfg
    )
  })
  names(out) <- scenarios
  structure(out, class = c("mrwin_scenario_set", "list"))
}

mrwin_run_simulation_grid <- function(
    scenarios = mrwin_scenarios(),
    n_iter = 1L,
    controls = mrwin_controls(n_strata = 4L, bootstrap = 10L, run_sdpd = FALSE),
    seed = 20260510L,
    run_sdpd = TRUE,
    sdpd_scale = c("aalen", "cox", "both"),
    run_benchmark = TRUE,
    benchmark_permutations = 50L
) {
  if (!inherits(scenarios, "mrwin_scenario_set")) {
    scenarios <- do.call(mrwin_scenarios, list(scenarios = scenarios))
  }
  controls <- .mrwin_as_controls(controls)
  n_iter <- .mrwin_grid_integer(n_iter, "n_iter", lower = 1L)
  seed <- .mrwin_grid_integer(seed, "seed", lower = -Inf)
  benchmark_permutations <- .mrwin_grid_integer(benchmark_permutations, "benchmark_permutations", lower = 1L)
  sdpd_scale <- match.arg(sdpd_scale)

  rows <- list()
  warnings <- list()
  row_id <- 1L
  for (s_idx in seq_along(scenarios)) {
    scen <- scenarios[[s_idx]]
    for (iter in seq_len(n_iter)) {
      iter_seed <- seed + s_idx * 10000L + iter
      cfg <- scen$config
      cfg$seed <- iter_seed
      one <- .mrwin_run_simulation_iteration(
        scenario = scen,
        iteration = iter,
        seed = iter_seed,
        config = cfg,
        controls = controls,
        run_sdpd = run_sdpd,
        sdpd_scale = sdpd_scale,
        run_benchmark = run_benchmark,
        benchmark_permutations = benchmark_permutations
      )
      rows[[row_id]] <- one$row
      if (length(one$warnings) > 0L) {
        warnings <- c(warnings, one$warnings)
      }
      row_id <- row_id + 1L
    }
  }

  results <- do.call(rbind, rows)
  rownames(results) <- NULL
  structure(
    list(
      results = results,
      warnings = .mrwin_simulation_warnings_df(warnings),
      scenarios = scenarios,
      n_iter = n_iter,
      controls = controls,
      seed = seed
    ),
    class = c("mrwin_simulation_grid", "list")
  )
}

mrwin_simulation_summary <- function(x) {
  if (!inherits(x, "mrwin_simulation_grid")) {
    stop("`x` must be an `mrwin_simulation_grid` object.", call. = FALSE)
  }
  dat <- x$results
  scenario_names <- unique(dat$scenario)
  rows <- lapply(scenario_names, function(name) {
    idx <- dat$scenario == name
    z <- dat[idx, , drop = FALSE]
    data.frame(
      scenario = name,
      n_iter = nrow(z),
      n_success = sum(z$success),
      mean_delta = mean(z$delta_gls, na.rm = TRUE),
      median_dscwr = stats::median(z$dscwr, na.rm = TRUE),
      rejection_rate = mean(z$p_delta < 0.05, na.rm = TRUE),
      sdpd_rejection_rate = mean(z$sdpd_rejected, na.rm = TRUE),
      discordant_warning_rate = mean(z$discordant_components, na.rm = TRUE),
      expected_direction = z$expected_direction[1L],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

mrwin_per_component_benchmark <- function(
    G,
    X,
    time,
    status,
    priorities = NULL,
    n_perm = 50L,
    seed = 1L
) {
  G <- as.matrix(G)
  X <- as.numeric(X)
  time <- as.matrix(time)
  status <- as.matrix(status)
  if (nrow(G) != length(X) || nrow(time) != nrow(G) || nrow(status) != nrow(G) ||
      ncol(time) != ncol(status)) {
    stop("`G`, `X`, `time`, and `status` dimensions are incompatible.", call. = FALSE)
  }
  n_perm <- .mrwin_grid_integer(n_perm, "n_perm", lower = 1L)
  seed <- .mrwin_grid_integer(seed, "seed", lower = -Inf)
  if (is.null(priorities)) {
    priorities <- colnames(time)
    if (is.null(priorities)) priorities <- paste0("priority_", seq_len(ncol(time)))
  }
  if (length(priorities) != ncol(time)) {
    stop("`priorities` must match the number of endpoint columns.", call. = FALSE)
  }

  beta_x <- .mrwin_lm_per_snp(G, X)
  component_rows <- list()
  row_id <- 1L
  for (k in seq_len(ncol(time))) {
    outcome <- tryCatch(
      mrwin_aalen_per_snp(G, time[, k], status[, k]),
      error = function(e) NULL
    )
    if (is.null(outcome)) next
    ivw <- .mrwin_summary_ivw(beta_x$beta, outcome$beta, outcome$se)
    egger <- .mrwin_summary_egger(beta_x$beta, outcome$beta, outcome$se)
    wm <- .mrwin_summary_weighted_median(beta_x$beta, outcome$beta, outcome$se)
    presso <- .mrwin_summary_presso(beta_x$beta, outcome$beta, outcome$se, n_perm = n_perm, seed = seed + k)
    for (entry in list(ivw, egger, wm, presso)) {
      component_rows[[row_id]] <- data.frame(
        priority = priorities[k],
        priority_index = k,
        method = entry$method,
        theta = entry$theta,
        se = entry$se,
        p_value = entry$p_value,
        intercept_p_value = entry$intercept_p_value,
        global_p_value = entry$global_p_value,
        stringsAsFactors = FALSE
      )
      row_id <- row_id + 1L
    }
  }
  components <- if (length(component_rows) == 0L) {
    data.frame()
  } else {
    do.call(rbind, component_rows)
  }
  pooled <- .mrwin_pool_component_methods(components)
  ivw_components <- components[components$method == "ivw" & is.finite(components$theta), , drop = FALSE]
  signs <- unique(sign(ivw_components$theta[ivw_components$theta != 0]))
  discordant <- length(signs) > 1L

  structure(
    list(
      exposure_summary = beta_x,
      components = components,
      pooled = pooled,
      discordant_components = discordant
    ),
    class = c("mrwin_component_benchmark", "list")
  )
}

print.mrwin_simulation_grid <- function(x, ...) {
  cat("mrwin simulation grid\n")
  cat("  scenarios:", length(x$scenarios), "\n")
  cat("  iterations per scenario:", x$n_iter, "\n")
  cat("  rows:", nrow(x$results), "\n")
  if (nrow(x$warnings) > 0L) {
    cat("  warnings:", nrow(x$warnings), "\n")
  }
  invisible(x)
}

summary.mrwin_simulation_grid <- function(object, ...) {
  mrwin_simulation_summary(object)
}

.mrwin_run_simulation_iteration <- function(
    scenario,
    iteration,
    seed,
    config,
    controls,
    run_sdpd,
    sdpd_scale,
    run_benchmark,
    benchmark_permutations
) {
  dat <- mrwin_simulate(config, seed = seed)
  fit <- tryCatch(
    mrwin(
      endpoint = mrwin_endpoint(dat$time, dat$status, colnames(dat$time)),
      genotype = dat$G,
      exposure = dat$X,
      gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
      controls = controls
    ),
    error = function(e) e
  )
  success <- !inherits(fit, "error")
  fit_warnings <- list()
  if (success) {
    fit_warnings <- .mrwin_simulation_fit_warnings(fit$warnings, scenario$name, iteration, seed)
  }

  sdpd_rejected <- NA
  sdpd_min_p <- NA_real_
  if (isTRUE(run_sdpd)) {
    sdpd <- tryCatch(
      mrwin_sdpd(dat$G, dat$X, dat$time[, 1L], dat$status[, 1L], scale = sdpd_scale, min_snps = 3L),
      error = function(e) NULL
    )
    if (!is.null(sdpd)) {
      sdpd_rejected <- isTRUE(sdpd$rejected)
      pvals <- vapply(sdpd$results, function(x) {
        if (is.null(x$intercept_p_value)) NA_real_ else x$intercept_p_value
      }, numeric(1))
      sdpd_min_p <- suppressWarnings(min(pvals, na.rm = TRUE))
      if (!is.finite(sdpd_min_p)) sdpd_min_p <- NA_real_
    }
  }

  benchmark <- NULL
  discordant <- identical(scenario$expected_direction, "mixed")
  pooled_ivw_theta <- pooled_ivw_p <- NA_real_
  if (isTRUE(run_benchmark)) {
    benchmark <- tryCatch(
      mrwin_per_component_benchmark(
        dat$G, dat$X, dat$time, dat$status,
        priorities = colnames(dat$time),
        n_perm = benchmark_permutations,
        seed = seed
      ),
      error = function(e) NULL
    )
    if (!is.null(benchmark)) {
      discordant <- isTRUE(discordant || benchmark$discordant_components)
      pooled <- benchmark$pooled
      pick <- if (nrow(pooled) > 0L && all(c("method", "pool") %in% names(pooled))) {
        pooled$method == "ivw" & pooled$pool == "inverse_variance"
      } else {
        logical()
      }
      if (length(pick) > 0L && any(pick)) {
        pooled_ivw_theta <- pooled$theta[pick][1L]
        pooled_ivw_p <- pooled$p_value[pick][1L]
      }
    }
  }
  extra_warnings <- list()
  if (discordant) {
    extra_warnings[[1L]] <- data.frame(
      scenario = scenario$name,
      iteration = iteration,
      seed = seed,
      code = "discordant_components",
      message = "Scenario or component benchmark indicates discordant endpoint directions.",
      stringsAsFactors = FALSE
    )
  }
  if (!success) {
    extra_warnings[[length(extra_warnings) + 1L]] <- data.frame(
      scenario = scenario$name,
      iteration = iteration,
      seed = seed,
      code = "fit_failed",
      message = conditionMessage(fit),
      stringsAsFactors = FALSE
    )
  }

  delta <- se <- dscwr <- ci_low <- ci_high <- q <- q_p <- p_delta <- NA_real_
  warning_count <- 0L
  if (success) {
    delta <- fit$point$delta_gls
    se <- fit$inference$se_delta_gls
    dscwr <- fit$point$dscwr
    ci_low <- fit$inference$ci95_dscwr[1L]
    ci_high <- fit$inference$ci95_dscwr[2L]
    q <- fit$heterogeneity$q
    q_p <- fit$heterogeneity$q_p_value
    p_delta <- 2 * stats::pnorm(-abs(delta / pmax(se, 1e-12)))
    warning_count <- length(fit$warnings)
  }

  row <- data.frame(
    scenario = scenario$name,
    package = scenario$package,
    iteration = iteration,
    seed = seed,
    success = success,
    expected_direction = scenario$expected_direction,
    true_effect = scenario$true_effect,
    delta_gls = delta,
    se_delta = se,
    dscwr = dscwr,
    ci_low = ci_low,
    ci_high = ci_high,
    p_delta = p_delta,
    q = q,
    q_p_value = q_p,
    sdpd_rejected = sdpd_rejected,
    sdpd_min_p_value = sdpd_min_p,
    component_ivw_theta = pooled_ivw_theta,
    component_ivw_p_value = pooled_ivw_p,
    discordant_components = discordant,
    warning_count = warning_count,
    stringsAsFactors = FALSE
  )
  list(row = row, warnings = c(fit_warnings, extra_warnings))
}

.mrwin_summary_ivw <- function(beta_x, beta_y, se_y) {
  keep <- .mrwin_benchmark_keep(beta_x, beta_y, se_y)
  if (sum(keep) < 2L) return(.mrwin_empty_benchmark("ivw"))
  bx <- beta_x[keep]
  by <- beta_y[keep]
  sy <- se_y[keep]
  w <- bx^2 / pmax(sy^2, 1e-20)
  theta <- sum(w * by / bx) / sum(w)
  se <- 1 / sqrt(sum(w))
  p <- 2 * stats::pnorm(-abs(theta / pmax(se, 1e-12)))
  list(method = "ivw", theta = theta, se = se, p_value = p, intercept_p_value = NA_real_, global_p_value = NA_real_)
}

.mrwin_summary_egger <- function(beta_x, beta_y, se_y) {
  fit <- tryCatch(mrwin_mr_egger(beta_x, beta_y, se_y), error = function(e) NULL)
  if (is.null(fit)) return(.mrwin_empty_benchmark("egger"))
  p <- 2 * stats::pnorm(-abs(fit$slope / pmax(fit$slope_se, 1e-12)))
  list(
    method = "egger",
    theta = fit$slope,
    se = fit$slope_se,
    p_value = p,
    intercept_p_value = fit$intercept_p_value,
    global_p_value = NA_real_
  )
}

.mrwin_summary_weighted_median <- function(beta_x, beta_y, se_y) {
  keep <- .mrwin_benchmark_keep(beta_x, beta_y, se_y)
  if (sum(keep) < 3L) return(.mrwin_empty_benchmark("weighted_median"))
  bx <- beta_x[keep]
  by <- beta_y[keep]
  sy <- se_y[keep]
  theta_snp <- by / bx
  w <- bx^2 / pmax(sy^2, 1e-20)
  theta <- .mrwin_weighted_quantile(theta_snp, w, 0.5)
  scale <- .mrwin_weighted_quantile(abs(theta_snp - theta), w, 0.5)
  se <- pmax(1.4826 * scale / sqrt(sum(keep)), 1e-8)
  p <- 2 * stats::pnorm(-abs(theta / se))
  list(method = "weighted_median", theta = theta, se = se, p_value = p, intercept_p_value = NA_real_, global_p_value = NA_real_)
}

.mrwin_summary_presso <- function(beta_x, beta_y, se_y, n_perm, seed) {
  keep <- .mrwin_benchmark_keep(beta_x, beta_y, se_y)
  if (sum(keep) < 3L) return(.mrwin_empty_benchmark("presso_global"))
  bx <- beta_x[keep]
  by <- beta_y[keep]
  sy <- se_y[keep]
  ivw <- .mrwin_summary_ivw(bx, by, sy)
  rss_obs <- sum(((by - ivw$theta * bx) / pmax(sy, 1e-20))^2)
  set.seed(seed)
  rss_perm <- numeric(n_perm)
  for (b in seq_len(n_perm)) {
    by_p <- sample(by, length(by), replace = FALSE)
    theta_p <- .mrwin_summary_ivw(bx, by_p, sy)$theta
    rss_perm[b] <- sum(((by_p - theta_p * bx) / pmax(sy, 1e-20))^2)
  }
  global_p <- mean(rss_perm >= rss_obs)
  list(method = "presso_global", theta = ivw$theta, se = ivw$se, p_value = ivw$p_value, intercept_p_value = NA_real_, global_p_value = global_p)
}

.mrwin_pool_component_methods <- function(components) {
  if (nrow(components) == 0L) {
    return(data.frame())
  }
  methods <- unique(components$method)
  rows <- list()
  row_id <- 1L
  for (method in methods) {
    dat <- components[components$method == method, , drop = FALSE]
    dat <- dat[is.finite(dat$theta) & is.finite(dat$se) & is.finite(dat$p_value), , drop = FALSE]
    if (nrow(dat) == 0L) next
    inv <- .mrwin_pool_inverse_variance(dat)
    bonf <- .mrwin_pool_bonferroni(dat)
    fish <- .mrwin_pool_fisher(dat)
    for (entry in list(inv, bonf, fish)) {
      rows[[row_id]] <- data.frame(
        method = method,
        pool = entry$pool,
        theta = entry$theta,
        se = entry$se,
        p_value = entry$p_value,
        stringsAsFactors = FALSE
      )
      row_id <- row_id + 1L
    }
  }
  if (length(rows) == 0L) data.frame() else do.call(rbind, rows)
}

.mrwin_pool_inverse_variance <- function(dat) {
  w <- 1 / pmax(dat$se^2, 1e-20)
  theta <- sum(w * dat$theta) / sum(w)
  se <- 1 / sqrt(sum(w))
  p <- 2 * stats::pnorm(-abs(theta / pmax(se, 1e-12)))
  list(pool = "inverse_variance", theta = theta, se = se, p_value = p)
}

.mrwin_pool_bonferroni <- function(dat) {
  idx <- which.min(dat$p_value)
  list(
    pool = "bonferroni",
    theta = dat$theta[idx],
    se = dat$se[idx],
    p_value = min(1, dat$p_value[idx] * nrow(dat))
  )
}

.mrwin_pool_fisher <- function(dat) {
  ps <- pmin(pmax(dat$p_value, 1e-300), 1 - 1e-12)
  stat <- -2 * sum(log(ps))
  p <- stats::pchisq(stat, df = 2 * length(ps), lower.tail = FALSE)
  list(pool = "fisher", theta = NA_real_, se = NA_real_, p_value = p)
}

.mrwin_benchmark_keep <- function(beta_x, beta_y, se_y) {
  is.finite(beta_x) & is.finite(beta_y) & is.finite(se_y) & se_y > 0 & abs(beta_x) > 1e-8
}

.mrwin_empty_benchmark <- function(method) {
  list(method = method, theta = NA_real_, se = NA_real_, p_value = NA_real_, intercept_p_value = NA_real_, global_p_value = NA_real_)
}

.mrwin_weighted_quantile <- function(x, w, prob) {
  ord <- order(x)
  x <- x[ord]
  w <- w[ord] / sum(w[ord])
  x[which(cumsum(w) >= prob)[1L]]
}

.mrwin_grid_integer <- function(x, name, lower) {
  x <- suppressWarnings(as.integer(x))
  if (length(x) != 1L || is.na(x) || x < lower) {
    stop("`", name, "` must be a single integer at least ", lower, ".", call. = FALSE)
  }
  x
}

.mrwin_simulation_fit_warnings <- function(warnings, scenario, iteration, seed) {
  if (length(warnings) == 0L) {
    return(list())
  }
  lapply(warnings, function(w) {
    data.frame(
      scenario = scenario,
      iteration = iteration,
      seed = seed,
      code = w$code,
      message = w$message,
      stringsAsFactors = FALSE
    )
  })
}

.mrwin_simulation_warnings_df <- function(warnings) {
  if (length(warnings) == 0L) {
    return(data.frame(
      scenario = character(),
      iteration = integer(),
      seed = integer(),
      code = character(),
      message = character(),
      stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, warnings)
  rownames(out) <- NULL
  out
}

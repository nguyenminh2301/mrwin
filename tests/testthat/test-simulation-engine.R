test_that("WP8 config validation and DGP produce stable package inputs", {
  expect_error(mrwin_config(alpha_x = c(0, 0)), "`alpha_x`")
  expect_error(mrwin_config(maf_low = 0.4, maf_high = 0.2), "`maf_low`")

  cfg <- mrwin_config(
    n_outcome = 36,
    m_snps = 5,
    baseline_haz = c(0.12, 0.18, 0.24),
    censoring_rate = 0.01,
    max_follow_up = 6,
    theta_f = 0.2,
    seed = 801
  )
  dat <- mrwin_simulate(cfg)

  expect_s3_class(cfg, "mrwin_config")
  expect_equal(dim(dat$G), c(36L, 5L))
  expect_equal(dim(dat$time), c(36L, 3L))
  expect_equal(dim(dat$status), c(36L, 3L))
  expect_true(all(is.finite(dat$G)))
  expect_true(all(vapply(seq_len(ncol(dat$G)), function(j) stats::var(dat$G[, j]) > 0, logical(1))))
  expect_true(all(dat$status %in% c(0L, 1L)))
  expect_true(all(dat$time[, 2L] <= dat$time[, 1L]))
  expect_true(all(dat$time[, 3L] <= dat$time[, 1L]))
})

test_that("WP8 scenarios define Package A-D metadata", {
  base <- mrwin_config(n_outcome = 30, m_snps = 4, seed = 802)
  scenarios <- mrwin_scenarios(base)

  expect_s3_class(scenarios, "mrwin_scenario_set")
  expect_named(scenarios, c("A_null", "B_valid_IV", "C_pleiotropy", "D_hierarchy_discordant"))
  expect_equal(scenarios$A_null$config$alpha_x, c(0, 0, 0))
  expect_equal(scenarios$C_pleiotropy$config$gamma_direct, c(0.05, 0, 0))
  expect_equal(scenarios$D_hierarchy_discordant$expected_direction, "mixed")
  expect_equal(scenarios$D_hierarchy_discordant$package, "D")
})

test_that("WP8 per-component benchmark returns component and pooled schemas", {
  cfg <- mrwin_config(
    n_outcome = 56,
    m_snps = 6,
    baseline_haz = c(0.15, 0.22, 0.30),
    censoring_rate = 0.01,
    max_follow_up = 6,
    seed = 803
  )
  dat <- mrwin_simulate(cfg)

  bench <- mrwin_per_component_benchmark(
    dat$G,
    dat$X,
    dat$time,
    dat$status,
    priorities = colnames(dat$time),
    n_perm = 10,
    seed = 804
  )

  expect_s3_class(bench, "mrwin_component_benchmark")
  expect_true(all(c("priority", "method", "theta", "se", "p_value") %in% names(bench$components)))
  expect_true(all(c("method", "pool", "theta", "p_value") %in% names(bench$pooled)))
  expect_true("ivw" %in% bench$components$method)
  expect_true("inverse_variance" %in% bench$pooled$pool)
  expect_type(bench$discordant_components, "logical")
})

test_that("WP8 simulation grid runs small scenarios and records caveats", {
  base <- mrwin_config(
    n_outcome = 42,
    m_snps = 5,
    baseline_haz = c(0.15, 0.22, 0.30),
    censoring_rate = 0.01,
    max_follow_up = 6,
    theta_f = 0.2,
    seed = 805
  )
  scenarios <- mrwin_scenarios(base, scenarios = c("A_null", "D_hierarchy_discordant"))

  grid <- mrwin_run_simulation_grid(
    scenarios = scenarios,
    n_iter = 1,
    controls = mrwin_controls(n_strata = 3, bootstrap = 6, seed = 806, run_sdpd = FALSE),
    seed = 807,
    run_sdpd = TRUE,
    sdpd_scale = "aalen",
    run_benchmark = TRUE,
    benchmark_permutations = 10
  )
  summary <- mrwin_simulation_summary(grid)

  expect_s3_class(grid, "mrwin_simulation_grid")
  expect_equal(nrow(grid$results), 2L)
  expect_true(all(grid$results$success))
  expect_true(all(c("delta_gls", "sdpd_rejected", "component_ivw_theta") %in% names(grid$results)))
  expect_equal(nrow(summary), 2L)
  expect_true("discordant_components" %in% grid$warnings$code)
  expect_equal(summary$expected_direction[summary$scenario == "D_hierarchy_discordant"], "mixed")
})

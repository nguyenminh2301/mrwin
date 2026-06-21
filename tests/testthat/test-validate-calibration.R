# Smoke test for the WP19 calibration harness that produces the external
# validation numbers recorded in inst/spec/validation-findings.md. Kept small
# (M = 3, tiny N) — this checks the harness contract, not the calibration rate
# (the rate itself is established by the recorded grid, not by CI).

test_that(".mrwin_validate_calibration returns a well-formed summary", {
  vc <- getFromNamespace(".mrwin_validate_calibration", "mrwin")
  res <- vc(
    M = 3L, N = 600L, m_snps = 20L, n_strata = 5L,
    scenario = "null", inference = "bootstrap", bootstrap = 40L, seed0 = 10L
  )
  expect_true(all(c(
    "scenario", "inference", "sigma_beta", "adjustment", "stratification",
    "M", "valid", "rate", "se", "weak_frac"
  ) %in% names(res)))
  expect_identical(res$scenario, "null")
  expect_identical(res$M, 3L)
  expect_gte(res$valid, 0L)
  expect_lte(res$valid, 3L)
  if (res$valid > 0L) {
    expect_gte(res$rate, 0)
    expect_lte(res$rate, 1)
  }
})

test_that(".mrwin_validate_calibration requires truth for coverage scenarios", {
  vc <- getFromNamespace(".mrwin_validate_calibration", "mrwin")
  expect_error(
    vc(M = 2L, N = 400L, scenario = "valid_iv", truth = NULL),
    "truth"
  )
})

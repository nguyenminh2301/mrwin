# WP6 Covariate Adjustment

Status: complete as of 2026-05-10.

## Scope

WP6 adds covariate adjustment to the R package without making Python a runtime dependency. The implemented production path is ordinal-IPTW over ordered PRS strata. GPS remains a reserved option and is blocked in `mrwin()` until a later work package defines the continuous-score model and fallback rules.

## Implemented Components

- `mrwin_propensity_weights()` fits a proportional-odds propensity model for stratum assignment using covariates.
- Stabilized weights use the weighted marginal stratum probability divided by the fitted stratum probability.
- Weights are truncated by probability limits, defaulting to the 1st and 99th percentiles.
- Per-stratum effective sample size is reported as `(sum w)^2 / sum(w^2)`.
- Strata with `ESS_d < ess_fraction * N_d` are dropped from the active contrast set.
- `mrwin_estimate()` accepts `active_strata` and builds adjacent or bridged contrasts across surviving labels.
- `mrwin_multiplier_bootstrap()` refits the propensity model inside each multiplier bootstrap iteration using the multiplier weights in the likelihood.
- High-level `mrwin()` now supports `adjustment = "ordinal_iptw"` and records ESS, active/dropped strata, weights, balance diagnostics, and positivity/bridging warnings.

## Acceptance Tests

The WP6 test fixture checks:

- ordinal-IPTW returns finite stabilized weights, balance diagnostics, and ESS-threshold failures;
- an internal dropped stratum creates a bridged contrast, for example `s3_vs_s1`;
- bootstrap refitting preserves the point-estimate bridged contrast schema and records per-iteration dropped strata;
- missing covariates, invalid truncation, zero-variance covariates, and deferred GPS mode fail explicitly.

## Remaining Review Items

- Add simulation stress tests for severe non-positivity and consecutive dropped strata.
- Review proportional-odds adequacy and GPS fallback design after WP8 simulation coverage is expanded.
- Evaluate whether balance diagnostics should include pairwise SMD tables or only max absolute SMD summaries in user reports.

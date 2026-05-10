# WP7 SDPD Diagnostics

Status: complete as of 2026-05-10.

## Scope

WP7 completes the package-level diagnostic layer for priority-1 direct pleiotropy. The implementation keeps R as the user interface and uses the Python prototype only as a reference design.

## Implemented Components

- `mrwin_sdpd()` is the public Step-Down Pleiotropy Diagnostic wrapper.
- `mrwin_aalen_per_snp()` extracts per-SNP additive-rate outcome summaries with HC0 sandwich standard errors.
- `mrwin_cox_per_snp()` extracts per-SNP Cox-PH log-HR summaries through a Breslow-style Newton solver.
- `mrwin_mr_egger()` performs random-effects MR-Egger and reports intercept, standard error, z statistic, p-value, slope, Q, phi, and valid SNP count.
- `mrwin_pleiotropy_bounded_ci()` validates inputs and reports sampling-only and pleiotropy-bounded intervals on log and DS-CWR scales.
- `mrwin()` now records SDPD underpower warnings when the valid SNP count is below the configured threshold.
- `mrwin_controls()` adds `sdpd_min_snps` and `pleiotropy_bias_radius`.

## Acceptance Tests

The WP7 tests check:

- deterministic MR-Egger null behavior and pleiotropy rejection behavior;
- Aalen per-SNP output against a hand sandwich calculation;
- Cox per-SNP output for finite log-HR estimates;
- public SDPD wrapper schema for both Aalen and Cox scales;
- validation failure modes for malformed inputs;
- high-level `mrwin()` SDPD underpower warnings and optional pleiotropy-bounded output.

## Deferred To WP8

Full Monte Carlo Type-I and power reproduction across v5 scenario grids belongs in WP8 because it requires a broader simulation engine, scenario grid runner, benchmark outputs, and runtime controls. WP7 locks the diagnostic functions and their contracts so WP8 can scale them.

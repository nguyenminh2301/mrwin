# WP8 Simulation Engine

Status: complete as of 2026-05-10.

## Scope

WP8 turns the existing DGP into a small, testable scenario-grid engine. It is intentionally not a heavy Monte Carlo runner by default. Slow full-size Type-I/power reproduction remains a release-validation task so ordinary package tests stay fast and deterministic.

## Implemented Components

- `mrwin_config()` now validates scalar/vector inputs with clear errors.
- `mrwin_simulate()` guards against zero-variance simulated SNP columns in small cohorts.
- `mrwin_scenarios()` defines Package A-D scenario metadata:
  - `A_null`
  - `B_valid_IV`
  - `C_pleiotropy`
  - `D_hierarchy_discordant`
- `mrwin_run_simulation_grid()` runs scenario/iteration combinations and records cCWR, SDPD, benchmark, and warning fields.
- `mrwin_simulation_summary()` aggregates scenario-level estimates, rejection rates, SDPD rejection rates, and discordance warning rates.
- `mrwin_per_component_benchmark()` creates a component-wise MR benchmark schema with IVW, MR-Egger, weighted median approximation, MR-PRESSO-style global test, inverse-variance pooling, Bonferroni pooling, and Fisher p-value pooling.

## Acceptance Tests

The WP8 tests check:

- DGP output dimensions, finite genotype matrix, binary statuses, and death-censoring of lower-priority times.
- scenario definitions and metadata for Package A-D;
- per-component benchmark component and pooled output schemas;
- small simulation-grid execution over null and hierarchy-discordant scenarios;
- explicit `discordant_components` warning capture.

## Remaining Slow Validation

These are not default unit tests yet:

- full Monte Carlo Type-I and power grids;
- large weak-instrument and pleiotropy grids;
- full Python/R parity fixtures for bootstrap, SDPD, DGP, and benchmark outputs;
- benchmark runtime and memory sweeps.

Those items should be handled through WP9 performance work and WP12 release QA gates.

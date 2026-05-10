# Project Checkpoint

Last updated: 2026-05-10.

This file tracks current implementation state. The roadmap remains canonical for work-package definitions: `inst/spec/work-package-roadmap.md`.

## Current State

| Area | State |
|---|---|
| Completed work packages | WP0-WP7 |
| Current work package | WP8 |
| Current branch | `main` |
| Remote | `origin/main` |
| Primary interface | R package |
| Python role | Reference/oracle harness only |
| Release readiness | Not release-ready; core path works, validation/reporting/performance/CI still pending |

## Completed Commits

| WP | Commit | Scope |
|---|---|---|
| Scaffold | `6d47b05` | Code-only package scaffold |
| Initial R implementation | `201e813` | Initial R package functions |
| WP0 | `919d3e3` | Algorithm specification |
| WP1 | `48f6875` | High-level user API |
| WP2 | `73fb895` | Data validation layer |
| WP3 | `bb0eff6` | Sparse kernel primitives |
| Roadmap | `bfec94b` | Canonical WP0-WP12 roadmap |
| WP4 | `519c753` | DS-CWR estimator |
| WP5 | `b5f54c2` | Bootstrap inference |
| WP6 | `9ee3dc7` | Covariate adjustment |
| WP7 | `792a7ce` | SDPD diagnostics |

## Latest Verification

Latest full verification was performed at the WP7 checkpoint:

- `testthat::test_local('tests/testthat')`: 156 passing tests.
- `pytest -q`: 39 passed, 22 skipped.
- `R CMD check --no-manual --no-build-vignettes` on source tarball: `Status: OK`.

Re-run these gates after each implementation work package and before any release candidate.

## Implemented Capability

- High-level `mrwin()` dense workflow.
- Input constructors: `mrwin_endpoint()`, `mrwin_gwas()`, `mrwin_controls()`.
- Data validation: endpoints, genotype, exposure, GWAS, covariates, strata, terminal-event consistency.
- Core kernel: dense pairwise hierarchical kernel and sparse stratum-pair primitives.
- Estimation: adjacent and bridged log-CWR, Delta-X, ISG, GLS, Q statistic.
- Inference: multiplier bootstrap, bivariate Delta covariance, Fieller sensitivity CI, bootstrap moment diagnostics.
- Adjustment: ordinal-IPTW, truncation, ESS diagnostics, positivity filtering, bridged active strata.
- Diagnostics: SDPD Aalen/Cox per-SNP summaries, MR-Egger intercept test, underpower/rejection warnings, pleiotropy-bounded CI helper.
- User methods: `print()`, `summary()`, `plot()` for `mrwin_fit`.

## Open Strategic Gaps

These are the main blockers before a release candidate:

- WP8: scenario-grid simulation engine and small-scale reproduction of Package A-D.
- WP8: Python/R oracle fixtures for DGP, bootstrap, SDPD, and benchmark schemas.
- WP8: simulation checks for Type-I, power, weak instruments, pleiotropy, and discordant components.
- WP9: performance backend and memory/time benchmarks.
- WP9: high-level sparse backend wiring; optional Rcpp/data.table path remains undecided.
- WP10: tidy outputs, richer plots, markdown report, clearer user caveats.
- WP11: README expansion, vignettes, examples, reference documentation.
- WP12: CI, coverage target, release checklist, install-from-GitHub verification.
- Statistical review: exact or documented Ledoit-Wolf shrinkage, AL-CWR secondary diagnostic, v5 Table 3 bias interpolation, GPS fallback decision.

## Next WP8 Work Plan

1. Freeze the simulation result schema: scenarios, seeds, metrics, and output table names.
2. Expand `mrwin_simulate()` or add a scenario runner without breaking existing toy tests.
3. Add small deterministic fixtures for null, valid-IV, pleiotropy, weak-instrument, and discordant-component scenarios.
4. Add Package A-D benchmark summaries at small scale.
5. Store oracle outputs under `inst/extdata/oracle/` when stable.
6. Keep slow Monte Carlo tests out of default `testthat` until CI strategy is defined.

## Documentation Rules

- README should stay user-facing and brief.
- `work-package-roadmap.md` defines the strategic sequence.
- `algorithm-spec.md` defines formulas, contracts, gaps, and review criteria.
- `project-checkpoint.md` records current progress and verification evidence.

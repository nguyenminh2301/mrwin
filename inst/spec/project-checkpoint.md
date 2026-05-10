# Project Checkpoint

Last updated: 2026-05-10.

This file tracks current implementation state. The roadmap remains canonical for work-package definitions: `inst/spec/work-package-roadmap.md`.

## Current State

| Area | State |
|---|---|
| Completed work packages | WP0-WP8 |
| Current work package | WP9 |
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
| WP8 | current checkpoint | Simulation engine and Package A-D small-grid schema |

## Latest Verification

Latest full verification was performed at the WP8 checkpoint:

- `testthat::test_local('tests/testthat')`: 186 passing tests.
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
- Simulation engine: Package A-D scenario definitions, small scenario-grid runner, scenario summary, per-component MR benchmark schema.
- User methods: `print()`, `summary()`, `plot()` for `mrwin_fit`.

## Open Strategic Gaps

These are the main blockers before a release candidate:

- WP9: performance backend and memory/time benchmarks.
- WP9: high-level sparse backend wiring; optional Rcpp/data.table path remains undecided.
- WP10: tidy outputs, richer plots, markdown report, clearer user caveats.
- WP11: README expansion, vignettes, examples, reference documentation.
- WP12: CI, coverage target, release checklist, install-from-GitHub verification.
- Slow validation: full-size Monte Carlo Type-I, power, weak-instrument, pleiotropy, and discordant-component grids.
- Statistical review: exact or documented Ledoit-Wolf shrinkage, AL-CWR secondary diagnostic, v5 Table 3 bias interpolation, GPS fallback decision.

## Next WP9 Work Plan

1. Decide the high-level sparse backend contract for `mrwin()`.
2. Benchmark dense kernel, sparse pair kernel, bootstrap, SDPD, and scenario-grid runtime at increasing N.
3. Add memory/time benchmark helpers with stable output schema.
4. Verify sparse and dense outputs match within tolerance on small fixtures.
5. Decide whether optional Rcpp/data.table is justified before release.
6. Keep large benchmarks out of default `testthat` until CI strategy is defined.

## Documentation Rules

- README should stay user-facing and brief.
- `work-package-roadmap.md` defines the strategic sequence.
- `algorithm-spec.md` defines formulas, contracts, gaps, and review criteria.
- `project-checkpoint.md` records current progress and verification evidence.

# Project Checkpoint

Last updated: 2026-05-10.

This file tracks current implementation state. The roadmap remains canonical for work-package definitions: `inst/spec/work-package-roadmap.md`.

## Current State

| Area | State |
|---|---|
| Completed work packages | WP0-WP12 (Phase I) |
| Current work package | Phase II (WP13-WP19): scalability, continuous estimator, theory |
| Current branch | `C` |
| Integration branch | `C` (Phase II); `main` holds Phase I |
| Primary interface | R package |
| Python role | Reference/oracle harness only |
| Release readiness | Phase I core path works; Phase II makes it biobank-scale and adds the continuous estimator + analytic variance before external release |

## Phase II Pointer

The second program of work is canonically defined in
`inst/spec/acceleration-roadmap.md`, with per-work-package specs in
`inst/spec/wp13-fast-kernel.md` … `inst/spec/wp19-scalability-validation.md`.
Phase II addresses the one structural gap Phase I left open: the win/loss pair
sweep is `Θ(N²/D)` per bootstrap iteration, which is infeasible at biobank
scale. Phase II replaces it with a subquadratic algorithm, removes the arbitrary
stratum-count `D` via a continuous estimator, and derives an analytic variance.

### Phase II progress

- **S1 (C1, single-endpoint fast kernel) — landed + R-verified 2026-06-14** on
  branch `C-wp13`. `R/kernel_fast.R` + `python/p1_engine_v5/kernel_fast.py`;
  parity bit-for-bit (Python) and full testthat suite under **R 4.3.3** (71
  groups, 0 failures) including `test-backend-fast.R` end-to-end; benchmark shows
  fast tail exponent 1.21 vs dense 2.30 (`inst/spec/benchmark-results.md`);
  `backend = "fast"` opt-in for K=1.
- **S2 (C2, hierarchical K=2 fast kernel) — landed + R-verified 2026-06-14** on
  `C-wp13`. Fenwick 2D dominance counter + four-regime tie-split; parity vs dense
  oracle (Python 28k random cohorts, 0 mismatches; R full suite 74 groups, 0
  failures); scaling exponent 1.18; `backend = "fast"` covers K∈{1,2}.
- **K=3 fast path (the v5 flagship) — landed + R-verified 2026-06-14** on
  `C-wp13`. CDQ 3D dominance counter + 16-regime level-3 tie-split; parity vs
  dense oracle; `backend = "fast"` covers K∈{1,2,3}; K≥4 falls back to dense.
- **Rcpp port — landed + R-verified 2026-06-18.** `src/fast_kernel.cpp`
  (`mrwin_fast_pair_cpp`), exact port of the validated kernel;
  `mrwin_fast_pair_win_loss` dispatches to it (pure-R fallback retained). The
  package now has compiled code (`LinkingTo: Rcpp`). Wall-clock: K=3 is 72×
  (n=2000) / 300× (n=5000) faster than dense, K=3 N=80k in 0.78 s, K=1 N=200k in
  0.083 s — biobank-scale. Differential testing (16k cohorts, 0 mismatches)
  caught and fixed a double-`eq` regime bug. Full suite 82 groups, 0 failures.

- **M3 (analytic variance) — sampling part landed + R-verified 2026-06-18.**
  `R/analytic_variance.R`: `mrwin_analytic_covariance()` (influence-function
  `cov_u = CᵀC`) + `mrwin_analytic_inference()` reusing the existing
  ISG/GLS/Fieller machinery. Validated vs the fixed-strata multiplier bootstrap:
  `se` ratio 0.99–1.00 (N∈{500,1500,4000}), no systematic bias.
  **Σ_gwas added 2026-06-18** via the exact GWAS-only resample
  (`mrwin_gwas_resample_covariance`; point estimate is piecewise-constant in β so
  has no pointwise gradient — the resample is the exact term, the dominant `xi`
  part stays analytic). Combined `[Σ_sampling + Σ_gwas]` matches the full
  bootstrap: `se` ratio 0.997–1.009 with σ_β∈{0.05,0.15,0.30} (GWAS share ~50%).
  **Wired end-to-end 2026-06-18:** `mrwin(inference="analytic")` (via
  `mrwin_analytic_bootstrap`, output-compatible with the bootstrap so the S3
  methods work). 7× faster than bootstrap at σ_β=0 (pure closed-form), ~1.5× at
  σ_β>0. **IPTW done 2026-06-18:** `adjustment="ordinal_iptw"` supported (weighted
  influence function + per-`β*` propensity-refit GWAS resample); it omits the
  estimated-weights correction (first-order), measured gap ≤~3% vs the bootstrap
  (either sign, shrinks with N; small because MR strata are ~independent of
  covariates). Bootstrap stays the exact default. **M3 complete.** Full suite 92
  groups, 0 failures.

### CRITICAL — R2 calibration finding (2026-06-18)

First external validation (type-I / coverage) found the package's **primary
bivariate-Delta CI is mis-calibrated** (type-I 0.000, coverage 0.995 — ~2× too
wide, no power); the **Fieller CI is correctly calibrated** (type-I 0.055,
coverage 0.955). Point estimator is consistent/unbiased with adequate
instruments. Bug is in the original WP4/WP5 ratio propagation, not Phase II;
`cov_u` is correct (Fieller uses it). See `inst/spec/validation-findings.md`.
**Fix (validated): make Fieller the primary interval — required before M1.**
Internal `analytic==bootstrap` parity never caught it (both reproduce the same
mis-calibrated interval).

### Direction (2026-06-18): the O(N²) goal is solved. The speed goal
is done; remaining value is methodological. Recommended path: **fix calibration
(Fieller) → M1 (continuous ISG) → WP19 (validation) → release**, with M2-wiring
in parallel; **WP14 and K≥4 are deferred** (WP14 is now constant-factor only and
largely superseded by M3 removing the bootstrap `B`-loop; v5 is K=3). ~23 steps
remain to a publishable + released package.
- **M2 (doubly-ranked stratification) core — landed + R-verified 2026-06-14.**
  `mrwin_doubly_ranked_strata()` in `R/strata.R`, exported, tested
  (`test-doubly-ranked.R`). Pipeline wiring (WP16 T2) still to do.
- Next options: K≥3 structural fast kernel; WP14 (build-once/reuse +
  incremental re-stratification); M2 T2 (wire doubly-ranked through the
  pipeline); or M3 (analytic variance).
- Verification environment: R 4.3.3 installed in the dev container; full
  testthat suite (79 groups, 0 failures) + Python differential tests green.

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
| WP9 | current checkpoint | Sparse backend, benchmark helpers, parity verification |
| WP10 | current checkpoint | User reporting: tidy output, multi-type plots, markdown report, formatted summary |
| WP11 | current checkpoint | Documentation: README expansion, 3 vignettes, reference manual |
| WP12 | current checkpoint | QA/Release: CI pipeline, coverage, release checklist, --as-cran verification |

## Latest Verification

Latest full verification was performed at the WP12 checkpoint:

- `testthat::test_dir('tests/testthat')`: 273 passing tests.
- `R CMD check --as-cran` on source tarball: `Status: 2 WARNINGs (expected vignette pre-built), 3 NOTEs (standard dev submission)`.
- CI pipeline: GitHub Actions configured for R CMD check (3 OS x 2 R versions), coverage via covr, and Python tests.

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
- Sparse backend: `mrwin_sparse_estimate()`, `mrwin_sparse_bootstrap()`, `mrwin_precompute_pair_kernels()`; wired into `mrwin()` via `backend = "sparse"`.
- Benchmark helpers: `mrwin_benchmark()` for runtime comparison across backends/sizes; `mrwin_verify_sparse_dense_parity()` for numerical equivalence checks.
- User methods: `print()`, `summary()`, `plot()` for `mrwin_fit`.
- User reporting: `tidy()` for broom-compatible data frames; `mrwin_report()` for text/markdown reports; `plot()` with three types: `isg`, `forest`, `bootstrap`; formatted `summary()` with Fieller, diagnostics, and caveats.

## Open Strategic Gaps

These are the remaining items after WP12 (deferred to post-release):

- Full-size Monte Carlo Type-I, power, weak-instrument, pleiotropy, and discordant-component validation grids.
- Statistical review: exact or documented Ledoit-Wolf shrinkage, AL-CWR secondary diagnostic, v5 Table 3 bias interpolation, GPS fallback decision.
- Rcpp decision: benchmarks show dense R backend is faster than pure R sparse up to N=2000 due to vectorization; Rcpp justified only for biobank-scale (N>10000) where memory becomes limiting.
- Coverage measurement: needs CI to run covr; target >= 80%.

## WP9 Completion Summary

1. Sparse backend contract decided: `backend = "sparse"` in `mrwin_controls()` routes to `mrwin_sparse_bootstrap()`.
2. Sparse backend wired into `mrwin()` workflow; removes the N×N kernel materialization.
3. Benchmark helpers added: `mrwin_benchmark()` and `mrwin_verify_sparse_dense_parity()`.
4. Parity verified: sparse and dense outputs match within `1e-10` tolerance.
5. Benchmark results (N=200..2000, B=50, M=20, D=5): dense is 2-3x faster than sparse due to R vectorization; sparse saves memory.
6. Rcpp decision: deferred to post-release; dense R is sufficient for typical cohort sizes (N<5000). Rcpp justified only for biobank-scale (N>10000).

## WP11 Completion Summary

1. README expanded with WP9/WP10 features: sparse backend, tidy output, multiple plot types, markdown reports.
2. Three vignettes created: simulation quickstart, real cohort template, interpreting diagnostics.
3. All exported functions documented via .Rd files (mrwin.Rd covers core API; dedicated .Rd files for WP9/WP10 functions).
4. DESCRIPTION updated with knitr/rmarkdown Suggests and VignetteBuilder.
5. All vignettes pass R CMD check (running R code and re-building outputs).

## WP12 Completion Summary

1. CI pipeline created: GitHub Actions with R CMD check on 3 OS x 2 R versions (release + devel).
2. Coverage workflow: covr integration with codecov upload.
3. Python tests: CI runs smoke, kernel, and replication tests on 3 OS x 2 Python versions.
4. Release checklist: `inst/spec/release-checklist.md` with 10 gate categories.
5. `R CMD check --as-cran` passes (2 expected WARNINGs for vignettes, 3 standard NOTEs).
6. covr added to Suggests for coverage tracking.

## Post-Release Validation Plan

1. Run full Monte Carlo Type-I error grid (null scenario, 1000 iterations).
2. Run power grid (valid-IV scenario, varying effect sizes).
3. Run pleiotropy detection grid (SDPD rejection rate vs gamma).
4. Run weak-instrument grid (Fieller bounded/unbounded rate).
5. Run discordant-component grid (warning emission rate).
6. Measure and report coverage via covr on CI.
7. Consider Rcpp sparse backend if biobank-scale demand emerges.

## Documentation Rules

- README should stay user-facing and brief.
- `work-package-roadmap.md` defines the strategic sequence.
- `algorithm-spec.md` defines formulas, contracts, gaps, and review criteria.
- `project-checkpoint.md` records current progress and verification evidence.

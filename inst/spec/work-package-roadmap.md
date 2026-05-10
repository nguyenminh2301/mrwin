# Work Package Roadmap

Status: canonical implementation roadmap.

Date recorded: 2026-05-09.

Progress: WP0, WP1, WP2, WP3, and WP4 are complete. The next implementation package is WP5.

Important sequencing rule: this roadmap defines WP0 plus 12 main work packages, WP1 through WP12. When the user asks for `WP4`, `WP5`, and later package names, use the package definitions in this file.

## Design Principles

1. R is the primary user interface. Python remains only as an oracle/reference during verification.
2. Every public function must have strong input validation, standardized output, and clear warnings.
3. The implementation is split into three layers:
   - `core`: kernel, CWR, ISG, GLS, bootstrap.
   - `diagnostics`: SDPD, MR-Egger, pleiotropy-bound, Q statistic, weak instrument checks.
   - `workflow`: end-user helper functions, plotting, summaries, and reporting.
4. Do not optimize prematurely with C++ before the algorithm is locked, but keep the design ready for a biobank-scale backend.
5. Every algorithm must have a test oracle: toy hand-calculated examples, Python parity, and simulation recovery.

## Work Packages

| WP | Content | Deliverable | Acceptance Criteria | Status |
|---|---|---|---|---|
| WP0 | Freeze specification from paper v5.1/v5.2 | `inst/spec/algorithm-spec.md` | All formulas and APIs map clearly to code | Complete |
| WP1 | User API design | `mrwin()`, `summary()`, `plot()`, S3 classes | User can run the pipeline through one main function | Complete |
| WP2 | Data validation layer | `mrwin_validate_data()` | Catches time/status, missing data, endpoint order, genotype, and PRS issues | Complete |
| WP3 | Core kernel | Dense kernel plus sparse pair kernel | Toy tests, antisymmetry, zero diagonal, Python parity | Complete |
| WP4 | DS-CWR estimator | Adjacent CWR, ISG, GLS, Q statistic | Matches Python on fixed seeds | Complete |
| WP5 | Bootstrap inference | Multiplier bootstrap, bivariate Delta, Fieller | Coverage simulation reaches the predefined threshold | Next |
| WP6 | Covariate adjustment | IPTW, ESS, truncation, bridging dropped strata | Balance diagnostics and positivity-failure tests pass | Pending |
| WP7 | SDPD diagnostics | Aalen/Cox per-SNP, MR-Egger, pleiotropy bounded CI | Type-I and power simulations reproduce the paper | Pending |
| WP8 | Simulation engine | v5 DGP, scenario grid, benchmarks | Reproduces Package A-D tables at small scale | Pending |
| WP9 | Performance backend | Sparse kernel, optional Rcpp/data.table | Memory/time benchmark passes with unchanged results | Pending |
| WP10 | User reporting | Tidy outputs, plots, markdown report | Results clearly show estimates, CI, Q, SDPD, and caveats | Pending |
| WP11 | Documentation | README, vignettes, examples, reference manual | New user can run the package in under 10 minutes | Pending |
| WP12 | QA/release | CI, R CMD check, coverage, review checklist | Release candidate passes all gates | Pending |

## Target API

Low-level API:

```r
mrwin_kernel()
mrwin_stratum_win_loss()
mrwin_estimate()
mrwin_multiplier_bootstrap()
mrwin_gls_pool()
mrwin_mr_egger()
mrwin_pleiotropy_bounded_ci()
```

End-user API:

```r
fit <- mrwin(
  data = dat,
  endpoint = mrwin_endpoint(
    time = c("t_death", "t_hf", "t_renal"),
    status = c("d_death", "d_hf", "d_renal"),
    priority = c("death", "hf", "renal")
  ),
  genotype = snp_cols,
  exposure = "X",
  beta_gwas = beta_hat,
  se_gwas = se_beta,
  covariates = c("age", "sex", "PC1", "PC2"),
  n_strata = 10,
  bootstrap = 500
)

summary(fit)
plot(fit)
mrwin_report(fit)
```

## Accuracy Tests

Unit tests:

- Kernel antisymmetry.
- Win/loss/tie hand-calculated examples.
- PRS strata deterministic behavior.
- GLS pooling with known covariance.
- MR-Egger with synthetic exact intercept.

Cross-language parity:

- Store five small seeds from the Python prototype.
- R output must match Python within tolerance.
- Oracle files live under `inst/extdata/oracle/`.

Simulation tests:

- Null: Type-I near 0.05.
- Valid IV: DS-CWR has the expected direction.
- Pleiotropy: SDPD detects signal when gamma is large enough.
- Weak instrument: Fieller becomes wide or unbounded as expected.
- Discordant components: output includes direction-interpretation warnings.

Data validation tests:

- Negative event time, missing values, and non-binary status.
- Death occurs before a lower-priority event while lower status is still 1.
- Empty stratum, small ESS, and extreme weights.
- Duplicate SNPs, zero-variance genotype, beta/se mismatch.

Performance tests:

- Dense backend at `N = 500`, `N = 1000`, and `N = 5000`.
- Sparse backend benchmark.
- Memory ceiling test.
- Bootstrap reproducibility with fixed seed.

## Review Plan

Code review:

- Do not merge public functions without tests.
- Every public function must define input contract, output contract, and failure modes.

Statistical review:

- Match every formula to the v5.2 methods.
- Recheck non-collapsibility caveat, SDPD caveat, and discordant-component caveat.
- Give special review to bivariate Delta, Fieller, Ledoit-Wolf, and IPTW bootstrap.

User review:

- Add a "minimal simulation" vignette.
- Add a "real cohort template" vignette.
- Add an "interpreting diagnostics" vignette.
- Users should understand outputs without reading the paper.

Release review:

- `R CMD check --as-cran`.
- `testthat`.
- Minimum 80% coverage target for R code.
- Snapshot output summaries.
- Clean package install from GitHub.

## Proposed Timeline

Phase 1, 1-2 weeks: lock spec, complete API, data validation, and docs skeleton.

Phase 2, 2-3 weeks: complete core estimator, bootstrap, and Python parity tests.

Phase 3, 2 weeks: complete SDPD, pleiotropy diagnostics, and simulation engine.

Phase 4, 2 weeks: complete IPTW, ESS, bridging, sparse kernel, and performance.

Phase 5, 1-2 weeks: complete vignettes, reports, CI, and release candidate.

## Release Gates

Do not release if any of these remain true:

- R package does not pass `R CMD check`.
- R results do not match Python oracle for toy/small simulation cases.
- Main function does not validate bad data.
- Bootstrap is not reproducible with a fixed seed.
- Output does not state caveats for weak instruments, pleiotropy risk, positivity failure, or discordant components.

The optimized sequence is: freeze specification, stabilize API, verify core, expand diagnostics, then optimize performance. This keeps the package convenient for users while reducing statistical error risk during expansion.

# Baseline Verification

Local verification baseline for the completion roadmap (stage G0.1,
`completion-roadmap-academic.vi.md` §6.1). Every gate runs on a local machine; no
cloud service is used. Reproduce with:

```bash
bash tools/baseline/run-baseline.sh            # all gates
bash tools/baseline/run-baseline.sh --offline  # skip network-dependent CRAN-incoming checks
```

Raw logs go to `out/baseline-<date>/` (gitignored). Record a new baseline here
after any change that should become the reference, and before every release
candidate.

## Baseline 2026-09-29 (G0)

Package version `0.1.1.9000`. The run used the working tree of commit `53e1089`
plus the G0 documentation and build-configuration changes, which are committed
together with this file. No R, C++ or test source changed in G0.

| Gate | Result | Verdict |
|---|---|---|
| testthat (source tree) | 15 files, 99 blocks, **887 expectations**; 0 failed, 0 errors, 0 skipped, 0 warnings | PASS |
| `R CMD build` | `mrwin_0.1.1.9000.tar.gz` | PASS |
| `R CMD check --as-cran --no-manual` | **0 ERROR, 0 WARNING, 3 NOTE** (classified below) | PASS |
| pytest (Python oracle) | 45 passed, 22 skipped | PASS |
| covr line coverage | **88.30%** (target >= 80%) | PASS |

### `R CMD check` NOTEs

| NOTE | Cause | Class | Action |
|---|---|---|---|
| CRAN incoming feasibility | "New submission"; "Version contains large components (0.1.1.9000)"; licence `MIT + file LICENSE` | Expected for an unreleased development version | Version becomes `1.0.0` at release (G3.5); see the LICENSE finding below |
| Future file timestamps: "unable to verify current time" | The check could not reach its time source | Environment only | None |
| Compilation flags: `-mno-omit-leaf-frame-pointer` | Default compiler flags of the distribution's R build | Environment only | None; CRAN's own flags do not include it |

The PDF reference manual was not built because the baseline machine has no
LaTeX installation (`--no-manual`). The Rd files themselves are checked
(`checking Rd files`, `checking Rd cross-references`, `checking for missing
documentation entries`: all OK). Build the manual on a machine with LaTeX
before submission (release checklist).

Fixed during G0 before this run (the first run had 1 WARNING and 4 NOTEs):

- A top-level-files NOTE listed the six README translations and the baseline
  output directory `out/`. Both were being bundled into the tarball. They are
  now excluded in `.Rbuildignore`.
- A WARNING came from the missing `en_US.UTF-8` locale on the baseline machine
  (environment only). The locale was installed.

### Python oracle skips

All 22 skips are in `tests/python/test_replication.py`. They check
`replication_outputs/`, which is generated only by `run_all.sh` (a long
simulation run) and is deliberately untracked. This is by design, not a failure.

### Coverage by file

| File | Line coverage |
|---|---:|
| `R/validate_calibration.R` | 0.00% |
| `R/validation.R` | 77.16% |
| `R/adjustment.R` | 79.33% |
| `R/api.R` | 81.36% |
| `R/kernel.R` | 82.07% |
| `R/benchmark.R` | 82.39% |
| `R/continuous_isg.R` | 88.73% |
| `R/methods.R` | 89.03% |
| `R/simulation_engine.R` | 89.18% |
| `R/sdpd.R` | 89.87% |
| `R/kernel_fast.R` | 90.12% |
| `R/estimate.R` | 91.05% |
| `R/config.R` | 92.59% |
| `R/analytic_variance.R` | 93.31% |
| `R/bootstrap.R` | 94.26% |
| `R/strata.R` | 95.56% |
| `R/backend_sparse.R` | 96.05% |
| `R/simulate.R` | 97.26% |
| `src/fast_kernel.cpp` | 99.58% |
| **Total** | **88.30%** |

`R/validate_calibration.R` is the internal external-calibration harness. It is
exercised by `tools/validation-scripts/`, not by unit tests. G2 will either add a
small smoke test or exclude it from coverage explicitly.

### Findings for later stages

| Finding | Stage |
|---|---|
| `print.mrwin_fit` prints `Bootstrap:100(100valid )` (missing spaces) | G3 |
| When the Fieller interval is unbounded, `tidy()` reports `ci_low`/`ci_high` as `exp(-50)`/`exp(50)` instead of `0`/`Inf`, and the row name is `low` | G3 |
| `LICENSE` contains the full MIT text. For `MIT + file LICENSE`, CRAN expects the two-line `YEAR:` / `COPYRIGHT HOLDER:` template | G3.5 |
| `LICENSE` names "Minh Nguyen and M. Aziz" as copyright holders, while `Authors@R` lists one author and no `cph` role; the maintainer should decide authorship and copyright roles | G3.5 |
| No vignettes exist (never committed); gate 4 of the release checklist cannot run yet | G3.4 |
| `checking examples ... NONE`: no Rd file has an `\examples` section. CRAN review usually asks for runnable examples for exported functions | G3.5 |

### Environment (`sessionInfo()` excerpt)

```
R version 4.3.3 (2024-02-29)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 24.04.4 LTS
BLAS/LAPACK: reference (libblas 3.12.0 / liblapack 3.12.0)
Rcpp 1.0.12, testthat 3.2.1, knitr 1.45, rmarkdown 2.25, covr 3.6.4, xml2 1.3.6
pandoc 3.1.3, qpdf 11.9.0, pdflatex: not installed
Python 3.11.15: numpy 2.4.6, scipy 1.17.1, statsmodels 0.15.0, pytest 9.1.1
```

This baseline was produced in the development container used for G0. The team
should re-run the same script on its own workstation. Matching gate results
(the same counts, coverage within rounding, and no new NOTE classes) confirm
that the baseline is portable.

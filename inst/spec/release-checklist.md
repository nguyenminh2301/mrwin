# Release Checklist

Last updated: 2026-09-29.

Use this checklist before any release candidate. All items must pass.

Gates 1-3 and 8 are automated by the local baseline script
`tools/baseline/run-baseline.sh` (no cloud services); its latest result is
recorded in `inst/spec/baseline-verification.md`.

## 1. R CMD check

```bash
R CMD build .
R CMD check mrwin_*.tar.gz --as-cran
```

- [ ] `Status: OK` (no ERRORs, no WARNINGs)
- [ ] All vignettes build successfully
- [ ] All tests pass

## 2. Local test gate

```r
devtools::load_all()
testthat::test_dir("tests/testthat")
```

- [ ] All tests pass (currently 887 expectations in 99 test blocks)
- [ ] No warnings from test execution

## 3. Coverage

```r
cov <- covr::package_coverage()
print(cov)
covr::report(cov)
```

- [ ] Line coverage >= 80%
- [ ] All exported functions have test coverage

## 4. Vignettes

```r
devtools::build_vignettes()
```

- [ ] `simulation-quickstart` builds without error
- [ ] `real-cohort-template` builds without error
- [ ] `interpreting-diagnostics` builds without error

## 5. Install from GitHub

```r
# In a fresh R session:
# install.packages("devtools")
devtools::install_github("nguyenminh2301/mrwin")
library(mrwin)
```

- [ ] Package installs without error
- [ ] `library(mrwin)` loads without error
- [ ] Quick start example runs end-to-end

## 6. Quick start verification

```r
library(mrwin)
cfg <- mrwin_config(n_outcome = 500, m_snps = 20, seed = 1)
dat <- mrwin_simulate(cfg)
fit <- mrwin(
  endpoint = mrwin_endpoint(dat$time, dat$status, c("death", "hf", "renal")),
  genotype = dat$G,
  exposure = dat$X,
  gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
  controls = mrwin_controls(n_strata = 5, bootstrap = 100, seed = 2)
)
print(fit)
summary(fit)
plot(fit)
tidy(fit)
mrwin_report(fit)
```

- [ ] `print(fit)` shows DS-CWR and CI
- [ ] `summary(fit)` shows formatted output
- [ ] `plot(fit)` produces ISG plot
- [ ] `tidy(fit)` returns one-row data frame
- [ ] `mrwin_report(fit)` produces report

## 7. Sparse backend verification

```r
fit_sparse <- mrwin(
  endpoint = mrwin_endpoint(dat$time, dat$status, c("death", "hf", "renal")),
  genotype = dat$G,
  exposure = dat$X,
  gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
  controls = mrwin_controls(n_strata = 5, bootstrap = 100, seed = 2, backend = "sparse")
)
```

- [ ] Sparse backend produces valid results
- [ ] Dense and sparse DS-CWR match within 1e-10

## 8. Python parity (optional for R release)

```bash
python -m pytest tests/python/ -v
```

- [ ] Python smoke tests pass
- [ ] Python kernel tests pass

## 9. Documentation

- [ ] README is up to date
- [ ] All exported functions have .Rd documentation
- [ ] No `R CMD check` warnings about undocumented objects

## 10. Statistical review

- [ ] DS-CWR formula matches paper v5.2
- [ ] Q statistic reference distribution is chisq(D-2)
- [ ] SDPD interpretation follows S8/S9 caveats
- [ ] Discordant-component warning follows v5.2 Section 6.4

## Current Status

As of 2026-09-29 (G0 local baseline, version `0.1.1.9000`):

| Gate | Status |
|---|---|
| R CMD check --as-cran | No ERROR; no package-caused WARNING; remaining NOTEs are expected for a development version or environment-only (see `baseline-verification.md`). PDF manual not built (no LaTeX on the baseline machine) |
| Tests (887 expectations, 99 blocks) | PASS (0 failures, 0 warnings, 0 skips) |
| Vignettes (3) | **Missing**: never committed (the `*.Rmd` gitignore rule excluded them; fixed in G0). To be written in G3.4 |
| Documentation | All 42 exports have Rd documentation; READMEs synchronised in G0 |
| Coverage | 88.30% line coverage (target >= 80%) |
| Python oracle tests | 45 passed, 22 skipped by design (replication outputs are generated only by `run_all.sh`) |
| Install from GitHub | Not yet tested in a clean environment |

## Remaining Work

- [ ] Write the three vignettes (G3.4) and restore gate 4
- [ ] Build the PDF manual on a machine with LaTeX before submission
- [ ] Verify install from GitHub in a clean environment
- [ ] Full Monte Carlo validation (G2, ADEMP protocol) before release
- [ ] Reporting polish found in G0: `print()` spacing, `tidy()` bounds and row name when the Fieller interval is unbounded (G3)

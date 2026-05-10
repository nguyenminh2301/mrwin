# Release Checklist

Last updated: 2026-05-10.

Use this checklist before any release candidate. All items must pass.

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

- [ ] All tests pass (currently 273)
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

As of 2026-05-10 (WP12 checkpoint):

| Gate | Status |
|---|---|
| R CMD check --as-cran | PASS (local) |
| Tests (273) | PASS |
| Vignettes (3) | PASS |
| Documentation | Complete |
| Coverage | Not yet measured (needs CI) |
| Install from GitHub | Not yet tested |

## Remaining Work

- [ ] Run `R CMD check --as-cran` on CI (multiple OS, R versions)
- [ ] Measure and report coverage via covr
- [ ] Verify install from GitHub in clean environment
- [ ] Full Monte Carlo validation (deferred to post-release)

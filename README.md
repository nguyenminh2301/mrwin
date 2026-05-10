# mrwin

`mrwin` is a work-in-progress R package for causal win-statistic methods on hierarchical composite endpoints.

The primary implementation is now under `R/`. The Python implementation is retained under `python/` as a prototype and validation harness while the R package matures.

## Current layout

```text
mrwin/
├── DESCRIPTION          R package metadata scaffold
├── NAMESPACE            R namespace scaffold
├── R/                   R package implementation
├── python/              Python prototype code used for validation
├── tests/
│   ├── python/          Python prototype tests
│   └── testthat/        R package tests
├── requirements.txt     Python prototype dependencies
├── pyproject.toml       Python prototype packaging/test config
└── run_all.sh           Optional Python prototype runner
```

## R API

Core exported functions:

- `mrwin()` runs the high-level individual-level workflow.
- `mrwin_endpoint()`, `mrwin_gwas()`, and `mrwin_controls()` define analysis inputs explicitly.
- `mrwin_validate_data()` checks endpoint, genotype, exposure, GWAS, covariates, strata, and terminal-event consistency before estimation.
- `mrwin_config()` and `mrwin_simulate()` create v5-style simulation inputs.
- `mrwin_kernel()` computes the hierarchical pairwise win/loss/tie kernel.
- `mrwin_pair_kernel()` and related pair functions compute sparse stratum-pair kernel blocks equivalent to dense kernel slices.
- `mrwin_estimate()` computes adjacent-stratum log-CWR and ISG point estimates.
- `mrwin_gls_pool()` pools adjacent ISG contrasts into DS-CWR and reports Q heterogeneity diagnostics.
- `mrwin_multiplier_bootstrap()` estimates DS-CWR, GLS standard error, bivariate-Delta CI, Fieller CI, Q heterogeneity diagnostics, and bootstrap moment diagnostics.
- `mrwin_propensity_weights()` computes ordinal-IPTW stabilized weights, ESS diagnostics, balance summaries, and positivity-filtered active strata.
- `mrwin_sdpd()` runs the Step-Down Pleiotropy Diagnostic on Aalen, Cox, or both scales.
- `mrwin_aalen_per_snp()`, `mrwin_cox_per_snp()`, and `mrwin_mr_egger()` expose the SDPD components.
- `mrwin_pleiotropy_bounded_ci()` widens the DS-CWR CI by an SDPD-implied bias band.

Minimal example:

```r
cfg <- mrwin_config(n_outcome = 500, m_snps = 20, seed = 1)
dat <- mrwin_simulate(cfg)

fit <- mrwin(
  endpoint = mrwin_endpoint(dat$time, dat$status, c("death", "hf", "renal")),
  genotype = dat$G,
  exposure = dat$X,
  gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
  controls = mrwin_controls(n_strata = 5, bootstrap = 100, seed = 2)
)

summary(fit)
plot(fit)
```

## Python Prototype

The Python code is not the intended long-term package interface. It is kept to preserve the current computational checks while the R implementation is written.

```bash
python -m pip install -r requirements.txt
python -m pytest
```

## Scope

The repository intentionally excludes manuscript drafts, replication outputs, and generated report files. Code and package tests belong here; paper text and publication artifacts should stay outside the repository or in a separate manuscript repository.

## Implementation Spec

The package implementation plan is frozen in [inst/spec/algorithm-spec.md](inst/spec/algorithm-spec.md). This document maps the paper v5.1/v5.2 methods to R package functions, test oracles, edge cases, and remaining implementation gaps.

The canonical work-package sequence is tracked in [inst/spec/work-package-roadmap.md](inst/spec/work-package-roadmap.md). It records WP0 plus the 12 main work packages WP1-WP12, with progress complete through WP7.

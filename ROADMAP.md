# mrwin roadmap

`mrwin` is built as a **programme**, not a single tool: one growing R package
that will host several methods papers and the functions each one introduces.
This file is the public place where planned work is parked so the structure is
ready before the code lands. Status legend: ✅ shipped · 🟡 in progress ·
⬜ planned · 💡 idea.

## Released (v0.1.0)

- ✅ Dose-standardized causal win ratio (DS-CWR): PRS strata → adjacent-stratum
  instrument-standardized gradients → GLS pooling.
- ✅ Subquadratic, compiled hierarchical win/loss kernel (biobank-scale).
- ✅ Analytic influence-function variance as a fast alternative to the bootstrap.
- ✅ Fieller confidence intervals for the ratio estimand.
- ✅ Ordinal inverse-probability-of-treatment weighting (covariate adjustment).
- ✅ Step-down pleiotropy diagnostic (SDPD) and per-SNP MR-Egger style checks.
- ✅ Simulation engine + scenario grids for calibration study.

## Next (planned)

These have a reserved home in the API and the documentation site
(`_pkgdown.yml`) so they slot in without reorganizing:

| status | feature | notes |
|--------|---------|-------|
| 🟡 | **Complete function reference** | roxygen2 docs for every export; `R CMD check` doc-clean |
| ⬜ | **Two-sample / summary-data mode** | run DS-CWR from GWAS summary statistics, no individual data |
| ⬜ | **Multivariable & mediation** | multiple exposures; decompose direct vs mediated win effects |
| ⬜ | **Plotting helpers** | forest plot of ISGs, dose–response curve, calibration plots |
| ⬜ | **Cluster-robust / family-structure variance** | related individuals in biobanks |
| 💡 | **K>3 priority levels on the fast path** | extend the compiled kernel beyond 3 endpoints |
| 💡 | **Shiny / web front-end** | point-and-click for non-programmers |

## Papers in the programme

Drafts live privately under `papers/<slug>/` (one folder per manuscript); only
public-facing summaries appear as documentation articles.

| status | slug | working title |
|--------|------|---------------|
| 🟡 | `01-methods-scalable-cwr` | Scalable & correctly-calibrated causal win statistics for hierarchical composite endpoints in MR |
| ⬜ | `02-…` | (reserved) |

## How to propose or claim an item

Open an issue describing the feature and the paper (if any) it supports, or see
`dev/` for the internal work-package process. New functionality follows one
rule from the project roadmap: **no documented claim before its test passes.**

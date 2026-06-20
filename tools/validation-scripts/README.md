# Validation / reproducibility scripts

Exploratory harnesses that produced the external-validation findings recorded in
`inst/spec/validation-findings.md`, `inst/spec/wp15-continuous-isg.md`, and the
WP18 ledger. They are **not** part of the package build (excluded via
`.Rbuildignore`); they are kept under version control so every number in the
specs and the manuscript is reproducible.

Each script expects the package installed (`R CMD INSTALL .`) and uses the
public API plus a few internals via `getFromNamespace(..., "mrwin")`. They are
single-core and seed-controlled.

| script | what it checks |
|---|---|
| `r2_type1.R`, `r2_coverage.R`, `r2_strong.R` | the original R2 calibration run (type-I / coverage) that exposed the mis-calibrated delta-method interval |
| `r2_cov_fieller.R`, `r2_fieller.R`, `r2_diag.R`, `r2_z.R` | Fieller construction, the ~2x delta-method over-coverage diagnosis, robust z-score |
| `m3_validate.R`, `gwas_validate.R` | analytic influence-function variance vs bootstrap (adjustment = none, sigma_beta > 0) |
| `iptw_gap.R`, `iptw_stress.R` | ordinal-IPTW analytic-vs-bootstrap agreement and stress cases |
| `m1_step2.R`, `m1_scale.R`, `m1_diag.R` | continuous-ISG (M1) exploration: boxcar reduction, variance vs the decile (negative finding) |
| `wp19_grid.R`, `wp19_cells.R` | the widened WP19 calibration grid; `wp19_cells.R` also drives the bootstrap-IPTW and doubly-ranked cells |
| `wp19_cells.results.txt` | recorded output of the WP19 additional cells (bootstrap-IPTW 0.017; doubly-ranked type-I 1.000) |

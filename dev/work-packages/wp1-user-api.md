# WP1 User API Specification

Status: implemented as the first high-level API layer.

## Goal

WP1 makes the package usable through one main function while preserving the lower-level numerical functions for testing and development.

Primary entry point:

```r
fit <- mrwin(
  data = cohort,
  endpoint = mrwin_endpoint(
    time = c("t_death", "t_hf", "t_renal"),
    status = c("d_death", "d_hf", "d_renal"),
    priority = c("death", "hf", "renal")
  ),
  genotype = snp_columns,
  exposure = "X",
  gwas = mrwin_gwas(beta = beta_hat, se = se_beta),
  controls = mrwin_controls(n_strata = 10, bootstrap = 500)
)
```

## Public Constructors

| Function | Purpose |
|---|---|
| `mrwin_endpoint()` | Declares endpoint time/status columns or matrices and clinical priority order. |
| `mrwin_gwas()` | Declares external GWAS weights, standard errors, SNP labels, and optional covariance. |
| `mrwin_controls()` | Declares strata, bootstrap count, seed, backend, adjustment mode, SDPD options, and tolerances. |

Constructors validate shape and intent before numerical estimation begins.

## Main Fit Object

`mrwin()` returns an object with class:

```r
c("mrwin_fit", "list")
```

Required fields:

- `call`
- `data_info`
- `endpoint_info`
- `instrument_info`
- `point`
- `inference`
- `heterogeneity`
- `diagnostics`
- `sdpd`
- `warnings`
- `controls`
- `bootstrap`
- `session_info`

The fit object intentionally stores the bootstrap object during early development so parity tests and debugging can inspect internals. A later performance release may add a `keep_bootstrap = FALSE` option.

## S3 Methods

| Method | Behavior |
|---|---|
| `print.mrwin_fit()` | Compact one-screen estimate, CI, and run metadata. |
| `summary.mrwin_fit()` | Returns estimate table, adjacent gradients, Q statistic, SDPD table, structured warnings. |
| `print.summary.mrwin_fit()` | Human-readable summary output. |
| `plot.mrwin_fit()` | Base R adjacent-gradient plot with GLS reference line. |

## WP1 Boundaries

WP1 deliberately blocks unsupported options:

- `backend = "sparse"` and `backend = "rcpp"` stop with a clear message.
- `adjustment = "ordinal_iptw"` and `adjustment = "gps"` stop with a clear message.
- Full GWAS covariance can be stored in `mrwin_gwas()` but WP1 bootstrap uses diagonal standard errors and adds a structured warning.

This is safer than silently running an incomplete method.

## Structured Warnings

`fit$warnings` is a list of records with:

- `code`
- `message`

Currently emitted codes:

| Code | Trigger |
|---|---|
| `covariates_ignored` | Covariates supplied while `adjustment = "none"`. |
| `full_covariance_not_used` | Full GWAS covariance supplied but diagonal bootstrap is used. |
| `weak_instrument` | Near-zero adjacent `Delta_X` or unbounded Fieller interval. |
| `sdpd_rejected` | SDPD MR-Egger intercept rejects on at least one scale. |

## Acceptance Criteria

WP1 is accepted when:

- users can run `mrwin()` with matrices or data-frame column names.
- fit object has stable `mrwin_fit` schema.
- `print()`, `summary()`, and `plot()` dispatch.
- unsupported advanced methods fail loudly.
- testthat coverage exists for constructors, one high-level fit, and unsupported options.

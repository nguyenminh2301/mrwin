# WP4 DS-CWR Estimator Specification

Status: implemented.

## Goal

WP4 makes the point-estimation layer a stable core API. It covers adjacent-stratum CWR, log-CWR, phenotypic shift, instrument-standardized gradients, GLS pooling, and Q heterogeneity diagnostics.

## Public Functions

| Function | WP4 role |
|---|---|
| `mrwin_estimate()` | Computes adjacent wins, losses, totals, CWR, log-CWR, Delta-X, and ISG. |
| `mrwin_gls_pool()` | Pools ISG contrasts into DS-CWR and computes standard error, confidence interval, GLS weights, and Q statistic. |

## Estimator Contract

`mrwin_estimate()` validates:

- matching `time` and `status` dimensions.
- matching rows across endpoint, genotype, and exposure.
- genotype columns matching `beta_hat`.
- finite `time`, `G`, `X`, and `beta_hat`.
- binary `status`.
- valid `n_strata`.
- optional `kernel` dimensions and values.
- optional non-negative finite weights.

The returned object has class `mrwin_estimate` and contains:

- `strata`, `score`, and `kernel`.
- adjacent `wins`, `losses`, `total`.
- adjacent `cwr`, `log_theta`, `delta_x`, and `delta_isg`.
- `weak_delta_x` flags for near-zero phenotypic shifts.

## GLS and Q Contract

`mrwin_gls_pool()` validates finite ISG values and a finite symmetric positive semi-definite covariance matrix. With `shrink = FALSE`, it performs exact GLS pooling:

```text
delta_GLS = (1^T Sigma^-1 1)^-1 1^T Sigma^-1 delta_ISG
Var(delta_GLS) = (1^T Sigma^-1 1)^-1
Q = (delta_ISG - 1 * delta_GLS)^T Sigma^-1 (delta_ISG - 1 * delta_GLS)
```

The Q reference degrees of freedom are `D - 2`, equivalently `length(delta_ISG) - 1`. With only one contrast, Q is reported as 0 and the p-value is `NA` because heterogeneity is not testable.

## Oracle

The fixed Python oracle is stored at `inst/extdata/oracle/wp4_estimator_small.R`. It was generated from deterministic arrays using the Python prototype kernel and the same GLS/Q formula. R tests verify:

- PRS score and strata labels.
- adjacent wins, losses, totals.
- CWR, log-CWR, Delta-X, and ISG.
- no-shrink GLS DS-CWR and Q statistic.

## Scope Boundary

WP4 does not certify multiplier-bootstrap inference, Fieller coverage, IPTW, ESS, bridging, SDPD, or performance backends. Those are covered by WP5-WP9.

# WP5 Bootstrap Inference Specification

Status: implemented.

## Goal

WP5 stabilizes multiplier-bootstrap inference for DS-CWR. It verifies that uncertainty is propagated through the joint bootstrap distribution of adjacent log-CWR numerators and Delta-X denominators, then extracted by bivariate Delta at the point estimates.

## Public Function

| Function | WP5 role |
|---|---|
| `mrwin_multiplier_bootstrap()` | Runs GWAS beta perturbation, PRS re-stratification, exponential multiplier weighting, bivariate-Delta covariance extraction, GLS pooling, Fieller sensitivity interval, Q statistic, and bootstrap moment diagnostics. |

## Bootstrap Contract

For each bootstrap iteration:

1. Draw or receive `beta_star`.
2. Recompute PRS and strata from `G %*% beta_star`.
3. Draw or receive multiplier weights `xi`.
4. Compute adjacent weighted log-CWR and Delta-X using pair weights `xi_i * xi_j`.

The returned object has class `mrwin_bootstrap` and contains:

- point estimates: `point_log_theta`, `point_cwr`, `point_delta_x`, `delta_isg`.
- pooled inference: `delta_gls`, `se_delta_gls`, `dscwr`, bivariate-Delta CI.
- Fieller sensitivity interval: `ci95_delta_fieller`, `ci95_dscwr_fieller`, `fieller_unbounded`.
- heterogeneity: `q`, `q_df`, `q_p_value`.
- covariance matrices: `cov_u`, `sigma_isg`, `sigma_lw`.
- bootstrap matrices: `bootstrap_log_theta`, `bootstrap_delta_x`, `valid_bootstrap`.
- moment diagnostics: kurtosis and skewness for log-CWR and Delta-X.

## Bivariate Delta Rule

WP5 explicitly does not compute covariance from per-iteration ratios. It forms:

```text
U = [LT | DX]
```

and computes `Cov(U)` across valid bootstrap iterations. The ISG covariance is then:

```text
Sigma_ISG[a,c] = g_a^T Cov((LT_a, DX_a), (LT_c, DX_c)) g_c
```

where:

```text
g_a = (1 / Delta_X_a, -log_theta_a / Delta_X_a^2)
```

and gradients are evaluated at point estimates.

## Fieller Rule

Fieller is reported as a sensitivity interval for the pooled gradient. If the quadratic has no bounded finite solution, the interval is `(-Inf, Inf)` and `fieller_unbounded = TRUE`.

## Test Hooks

`mrwin_multiplier_bootstrap()` accepts optional deterministic `beta_draws` and `multiplier_weights` matrices. These are for tests and oracle fixtures, not ordinary user workflows. They make the bootstrap path exactly reproducible without relying on cross-language RNG identity.

## Scope Boundary

WP5 does not implement IPTW, ESS, positivity bridging, production coverage simulation grids, or sparse/Rcpp bootstrap backends. Those remain WP6-WP9 work.

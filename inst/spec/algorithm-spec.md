# mrwin WP0 Algorithm Specification

Version: WP0 freeze, 2026-05-09
Repository target: `mrwin` R package
Primary manuscript source: Paper 1 v5.2, with v5 and v5.1 review deltas reconciled
Implementation policy: R is the user-facing package; Python remains a reference oracle until parity is proven.

## 1. Source Audit

The package specification is frozen from these local sources:

| Source | Role in package spec |
|---|---|
| `paper1_v5/CausalWinStatistic_v5.md` | Baseline main-text methods: kernel, estimand, GLS/ISG, multiplier bootstrap, IPTW, bridging, SDPD, DGP, application. |
| `paper1_v5/CausalWinStatistic_v5_Supplement.md` | Baseline supplement S1-S7, especially S7 algorithmic pseudocode. |
| `paper1_v5.1/docs/CausalWinStatistic_v5_2.md` | Current main-text source of truth. Adds benchmark caveats and v5.2 wording. |
| `paper1_v5.1/docs/CausalWinStatistic_v5_2_Supplement.md` | Current supplement source of truth, S1-S13. |
| `paper1_v5.1/docs/v5_1_to_v5_2_changelog.md` | Required v5.2 deltas: Package D, discordant-component caveat, bridging edge case. |
| `paper1_v5.1/docs/p1_v51_adversarial_review.md` | Review-driven failure modes that package must address. |
| `paper1_v5.1/docs/Package_D_Aziz_Benchmark.md` | Per-component MR benchmark design and interpretation requirements. |
| `mrwin/python/p1_engine_v5/*.py` | Current executable reference implementation for small-scale oracle tests. |

No manuscript drafts or generated outputs are package artifacts. This spec is the bridge from manuscript methods to production package behavior.

## 2. Package Scope

The package must support two workflows.

1. Individual-level cCWR workflow:
   - Inputs: genotype matrix `G`, exposure `X`, hierarchical composite endpoint times/statuses, optional covariates `Z`, GWAS weights `beta_hat`, GWAS uncertainty `sigma_beta` or covariance.
   - Outputs: DS-CWR, log-scale GLS estimate, bivariate-Delta CI, Fieller sensitivity CI, Q heterogeneity statistic, SDPD diagnostics, warnings, and reproducibility metadata.

2. Simulation and validation workflow:
   - Inputs: v5 DGP parameters and scenario grids.
   - Outputs: simulation cohorts, benchmark summaries, Python parity fixtures, coverage/type-I/power checks.

Summary-data-only estimation is out of scope for the first complete release except as a future extension hook.

## 3. Core Objects and Data Contracts

### 3.1 Endpoint Contract

Endpoint data are matrices:

- `time`: numeric `N x K`, observed times per priority, finite and non-negative.
- `status`: integer/logical `N x K`, event indicators in `{0, 1}`.
- Columns are ordered from highest clinical priority to lowest priority.

For cardiorenal v5 simulations:

1. priority 1: death, absorbing.
2. priority 2: heart failure hospitalization.
3. priority 3: renal decline.

Validation requirements:

- `dim(time) == dim(status)`.
- no missing `time` or `status` in the analysis set.
- all `time > 0` for simulated survival data; real data may permit `time == 0` only if explicitly allowed.
- `status` must be binary.
- priority order must be declared by the user-facing workflow; do not infer clinical order from column names.
- if a terminal event occurs before a lower-priority event, the lower-priority event must either be censored or explicitly handled by a competing-risk preprocessing rule.

### 3.2 Genotype and Instrument Contract

- `G`: numeric `N x M` matrix, SNPs in columns.
- `beta_hat`: numeric length `M`, external GWAS weights used to compute PRS.
- `sigma_beta`: numeric length `M` for diagonal GWAS covariance, or future full `M x M` covariance.
- SNPs used with diagonal `sigma_beta` must be LD-clumped; package must warn that residual LD can understate uncertainty.

Validation requirements:

- no zero-variance SNP columns.
- `ncol(G) == length(beta_hat)`.
- `length(sigma_beta) == length(beta_hat)` for diagonal mode.
- reject or warn on missing genotype values unless an imputation rule is explicitly supplied.

### 3.3 Covariate Contract

`Z` contains PCs, age, sex, or other covariates used for IPTW/GPS adjustment.

Validation requirements:

- `nrow(Z) == N`.
- no missing values unless an imputation/preprocessing object is supplied.
- categorical covariates must be encoded through a documented model matrix transformation.

## 4. Algorithm Specification

### 4.1 Hierarchical Comparison Kernel

For individuals `i` and `j`, evaluate endpoint priorities `k = 1, ..., K`.

At priority `k`:

- `i` wins if `j` experienced event `k` and `time_i,k > time_j,k`.
- `i` loses if `i` experienced event `k` and `time_j,k > time_i,k`.
- otherwise the pair is tied at priority `k`.

The comparison terminates at the first informative priority. If no priority is informative, the pair ties.

Required properties:

- output values in `{-1, 0, +1}`.
- strict antisymmetry: `H[i, j] == -H[j, i]`.
- zero diagonal.
- no use of stratum assignment inside the kernel itself.

Package API:

- low-level dense: `mrwin_kernel(time, status, block_size)`.
- low-level block/sparse backend: to be added, equivalent to Python `kernel_sparse.py`.
- weighted aggregation: `mrwin_stratum_win_loss(kernel, idx_high, idx_low, weights)`.
- log contrast: `mrwin_stratum_log_cwr(...)`.

Acceptance tests:

- hand-calculated 3-person endpoint examples.
- antisymmetry and diagonal tests.
- dense-vs-sparse equality on the same stratum pair.
- Python parity against `python/p1_engine_v5/kernel.py`.

### 4.2 PRS and Strata

Compute:

```text
S_i = G_i^T beta_hat
```

Assign individuals to `D` ordered strata by PRS rank/quantile. Default `D = 10`.

Required behavior:

- deterministic tie handling.
- strata labels are 1-indexed in R.
- each stratum should be non-empty under normal operation; empty strata trigger a validation error or bridging workflow depending on context.

Package API:

- `mrwin_prs_strata(G, beta, n_strata)`.

Acceptance tests:

- deterministic ties.
- exact stratum counts for toy vectors.
- Python parity for fixed PRS values.

### 4.3 Adjacent-Stratum CWR and ISG

For adjacent strata `(d, d-1)`:

```text
tau_1(d,d-1)  = weighted proportion of wins by stratum d over d-1
tau_-1(d,d-1) = weighted proportion of losses by stratum d over d-1
theta(d,d-1)  = tau_1 / tau_-1
log_theta     = log(theta)
Delta_X       = mean(X | d) - mean(X | d-1)
delta_ISG     = log_theta / Delta_X
```

Implementation notes:

- win/loss counts are enough because total pair weights cancel in `tau_1 / tau_-1`.
- protect logs with a tiny floor only for computational stability; expose when floors are used.
- if `Delta_X` is near zero, flag weak-instrument behavior.

Package API:

- point estimates: `mrwin_estimate(...)`.
- user-facing result class later: `mrwin_fit`.

Acceptance tests:

- hand-calculated CWR example.
- weighted vs unweighted consistency when all weights are 1.
- weak-instrument warning when `abs(Delta_X)` is below tolerance.

### 4.4 Multiplier Bootstrap

For bootstrap iteration `b = 1, ..., B`:

1. Draw GWAS weights:

```text
beta_star^(b) ~ Normal(beta_hat, Sigma_GWAS)
```

Initial implementation supports diagonal `Sigma_GWAS = diag(sigma_beta^2)`.

2. Recompute PRS and strata:

```text
S_star_i = G_i^T beta_star
d_star_i = quantile_stratum(S_star_i)
```

3. Draw multiplier weights:

```text
xi_i ~ Exp(1)
```

4. If covariate adjustment is enabled, refit the propensity model inside the bootstrap with `xi`-weighted likelihood and compute stabilized weights.

5. For each surviving adjacent/bridged contrast, compute:

```text
log_theta_star^(b,d)
Delta_X_star^(b,d)
```

using pair weights `(xi_i * omega_i) * (xi_j * omega_j)`.

Required output from bootstrap:

- matrix `LT`: `B x (D-1)` perturbed log-theta values.
- matrix `DX`: `B x (D-1)` perturbed phenotypic shifts.
- valid iteration count.
- moment diagnostics for `LT` and `DX`: kurtosis and skewness.

Package API:

- `mrwin_multiplier_bootstrap(...)`.

Acceptance tests:

- reproducible output for fixed seed.
- no per-iteration ratio used for covariance.
- valid iteration count is reported.
- Python parity on small fixed seeds.

### 4.5 Bivariate Extraction and ISG Covariance

Form the joint bootstrap matrix:

```text
U = [LT | DX]
```

Then compute the empirical covariance:

```text
Cov(U) = Sigma_U, dimension 2(D-1) x 2(D-1)
```

Use first-order multivariate Delta for:

```text
delta_ISG_d = LT_d / DX_d
```

For contrasts `a` and `c`, define gradients:

```text
g_a = (1 / Delta_X_a, -log_theta_a / Delta_X_a^2)
g_c = (1 / Delta_X_c, -log_theta_c / Delta_X_c^2)
```

Then:

```text
Sigma_ISG[a,c] = g_a^T Cov((LT_a, DX_a), (LT_c, DX_c)) g_c
```

This covariance is evaluated at point estimates, not per-bootstrap ratios.

Acceptance tests:

- compare to manual 2-contrast covariance calculation.
- reject/flag covariance with non-finite entries.
- confirm dimensions under ordinary and bridged strata.

### 4.6 Ledoit-Wolf Shrinkage

Use a shrinkage covariance:

```text
mu = trace(Sigma_ISG) / p
Sigma_LW = (1 - rho) Sigma_ISG + rho * mu * I
```

where `p = D - 1` or the number of surviving/bridged contrasts.

The current R implementation has a simple stabilized shrinkage approximation. The complete package must either:

- implement the Ledoit-Wolf estimator exactly enough for production, or
- document the approximation and validate coverage by simulation.

Acceptance tests:

- symmetric positive semi-definite output.
- improved condition number when `Sigma_ISG` is ill-conditioned.
- parity or documented tolerance against Python reference.

### 4.7 GLS Pooling and DS-CWR

Given `delta_ISG` and `Sigma_LW`:

```text
delta_GLS = (1^T Sigma_LW^-1 1)^-1 1^T Sigma_LW^-1 delta_ISG
Var(delta_GLS) = (1^T Sigma_LW^-1 1)^-1
DS-CWR = exp(delta_GLS)
```

Primary CI:

```text
delta_GLS +/- 1.96 * se_delta_GLS
exp(bounds)
```

Secondary diagnostic:

- AL-CWR, pooling raw log-CWR without ISG scaling, must be added.
- divergence between AL-CWR and DS-CWR should be reported as sensitivity to `D`.

Package API:

- `mrwin_gls_pool(delta_isg, sigma, shrink = TRUE)`.

Acceptance tests:

- manual GLS example with diagonal covariance.
- pseudoinverse behavior for near-singular covariance.
- Q statistic p-value against `chisq(df = D-2)`.

### 4.8 Heterogeneity Q Statistic

Compute:

```text
Q = (delta_ISG - 1 * delta_GLS)^T Sigma_LW^-1 (delta_ISG - 1 * delta_GLS)
```

Reference distribution:

```text
Q ~ chisq(D - 2)
```

Interpretation:

- The null is gradient homogeneity across adjacent local gradients.
- This does not imply transitivity of the win-ratio relation.
- Do not use label permutation as the reference distribution because it breaks adjacent-stratum covariance.

Acceptance tests:

- Q = 0 when all `delta_ISG` values equal.
- simulated Gaussian null check in a slow/integration test.

### 4.9 Fieller Sensitivity CI

The primary interval is bivariate-Delta. Fieller is a sensitivity diagnostic for weak instruments.

For pooled numerator `U1_bar` and denominator `U2_bar`, solve:

```text
(U1_bar - delta * U2_bar)^2 < z^2 Var(U1_bar - delta * U2_bar)
```

If the quadratic has no bounded finite solution, report an unbounded Fieller interval and flag weak instrument behavior.

Acceptance tests:

- bounded interval in non-weak toy case.
- unbounded interval when denominator variance makes the Fieller denominator non-positive.

### 4.10 IPTW and Positivity

Default adjustment model:

- ordinal logistic propensity of stratum on covariates.
- weights:

```text
omega_i = pi_d / P(d_i | Z_i)
```

Bootstrap requirement:

- refit propensity in every bootstrap iteration using `xi`-weighted likelihood.
- do not estimate propensity once and reuse it through bootstrap.

Diagnostics:

- stabilized weight distribution.
- truncation at 1st/99th percentiles by default.
- effective sample size:

```text
ESS_d = (sum_i omega_i)^2 / sum_i omega_i^2
```

Failure rule:

- if `ESS_d < 0.5 * N_d`, excise stratum `d` and apply bridging.

Current status:

- current R implementation does not yet implement production IPTW/GPS.
- next work package must implement `mrwin_propensity_weights()`, ESS diagnostics, and bootstrap refitting.

Acceptance tests:

- balance diagnostics before/after weighting.
- extreme weights trigger warnings.
- proportional-odds failure path supports GPS fallback later.

### 4.11 Bridging Under Positivity Failure

If an internal stratum `d` is excised, replace contrasts:

```text
(d+1, d) and (d, d-1) -> (d+1, d-1)
```

Scale ISG by:

```text
Delta_X(d+1, d-1)
```

If consecutive strata are excised, bridge across nearest surviving neighbors, e.g.:

```text
(d+1, d-2)
```

Covariance is not manually specified. It is recovered from the same bivariate bootstrap stack over surviving/bridged contrasts.

Acceptance tests:

- single dropped stratum gives one bridged contrast.
- consecutive dropped strata reduce dimension correctly.
- off-diagonal vanishes when no stratum overlap exists.
- bootstrap covariance dimensions match surviving contrast count.

### 4.12 SDPD Diagnostic

The Step-Down Pleiotropy Diagnostic evaluates priority-1 direct pleiotropy.

Algorithm:

1. Extract per-SNP exposure summaries `beta_X`.
2. Extract per-SNP priority-1 outcome summaries:
   - Aalen additive hazard scale for collapsibility.
   - Cox-PH variant as sensitivity and preferred scale when prevalent-sample ascertainment is plausible.
3. Fit random-effects MR-Egger:

```text
beta_Y_m = theta_0 + theta_1 beta_X_m + error_m
```

4. Test:

```text
H0: theta_0 = 0
```

If rejected, mark the cCWR analysis invalid due to hierarchy-contaminating pleiotropy unless a domain-specific override is explicitly supplied.

Important interpretation:

- SDPD is necessary but not sufficient.
- It has low power at moderate SNP counts.
- It is most meaningful for priority 1, the absorbing endpoint.
- Lower-priority marginal MR-Egger can be confounded by competing-risk structure.

Package API:

- `mrwin_aalen_per_snp(G, time, status)`.
- `mrwin_cox_per_snp(G, time, status)`.
- `mrwin_mr_egger(beta_x, beta_y, se_y)`.
- future high-level: `mrwin_sdpd(...)`.

Acceptance tests:

- MR-Egger exact synthetic intercept.
- Type-I simulation under `gamma_1 = 0`.
- power simulation under `gamma_1 > 0`.
- ascertained vs unascertained scenario reproduces v5.2 S8 direction.

### 4.13 Pleiotropy-Bounded CI

If SDPD does not reject but power is incomplete, widen the sampling CI by the bias band implied by the detectable `gamma_1` threshold.

Inputs:

- `delta_hat`.
- `se_delta`.
- `bias_radius`, or later a lookup/interpolation from v5 Table 3 bias curve.

Output:

- sampling CI on log scale and DS-CWR scale.
- pleiotropy-bounded CI on both scales.
- null-crossing indicators.

Package API:

- `mrwin_pleiotropy_bounded_ci(delta_hat, se_delta, bias_radius, level)`.

Acceptance tests:

- bounded CI strictly contains sampling CI.
- null-crossing flags correct.
- v5 Table 3 interpolation fixture to be added.

### 4.14 DGP Specification

The simulation DGP is an illness-death model with shared gamma frailty.

State structure:

- healthy -> HF.
- healthy -> renal decline.
- healthy -> death.
- HF -> death.
- renal decline -> death.

Operational simplified package DGP currently samples latent priority event times from Weibull transition-style intensities and censors lower-priority events at death.

Parameters:

- `G_im ~ Binomial(2, maf_m)` and standardized.
- `maf_m ~ Uniform(0.1, 0.4)`.
- `beta_m ~ Normal(0, sigma_beta^2)`.
- `S_i = G_i^T beta`.
- `U_i ~ Normal(0, 1)`.
- `W_i ~ Gamma(1/theta_F, theta_F)` with mean 1 and variance `theta_F`; if `theta_F = 0`, `W_i = 1`.
- `X_i = alpha_S S_i + alpha_U U_i + error_i`.
- event hazards follow Weibull-PH style linear predictors with `alpha_X`, `nu_U`, `gamma_direct`, and `log(W)`.
- administrative censoring is exponential capped at `max_follow_up`.

Package API:

- `mrwin_config(...)`.
- `mrwin_simulate(config, seed)`.

Acceptance tests:

- output dimensions.
- binary status.
- deterministic with fixed seed.
- death censors lower-priority events.
- parameter scenarios reproduce expected directional changes.

### 4.15 Per-Component MR Benchmark

v5.2 adds Package D to compare cCWR against ordinary component-wise MR.

Benchmark estimators:

- IVW.
- MR-Egger.
- weighted median.
- MR-PRESSO global/outlier workflow.

Pooling rules:

- Bonferroni: smallest component p-value adjusted by `K`.
- inverse-variance pooling across components.
- Fisher combined p-value.

Required scenarios:

| Scenario | `alpha_X` | `gamma_1` | Meaning |
|---|---:|---:|---|
| A_null | `(0, 0, 0)` | 0 | calibration |
| B_valid_IV | `(-0.4, -0.4, -0.4)` | 0 | valid causal effect |
| C_pleiotropy | `(-0.4, -0.4, -0.4)` | 0.05 | direct priority-1 pleiotropy |
| D_hierarchy_discordant | `(+0.4, 0, -0.4)` | 0 | opposing component effects |

Interpretation requirement:

- if components move in conflicting directions across priorities, cCWR direction must be qualified by component-specific MR.
- pooled per-component MR can give misleading directional signals under hierarchy discordance.

Current status:

- Python benchmark exists.
- R package does not yet implement benchmark estimators beyond MR-Egger core.

Acceptance tests:

- per-component Aalen summaries match Python for small seeds.
- benchmark result schema stable.
- discordant scenario emits interpretation warning.

## 5. Public API Target

### 5.1 Current Low-Level API

Already present or partially present:

| Function | Status | Notes |
|---|---|---|
| `mrwin()` | implemented | WP1 high-level dense/no-adjustment workflow. Sparse and IPTW/GPS are blocked until later work packages. |
| `mrwin_endpoint()` | implemented | Declares endpoint columns/matrices and priority order. |
| `mrwin_gwas()` | implemented | Declares GWAS beta/se/covariance; full covariance stored but diagonal bootstrap only in WP1. |
| `mrwin_controls()` | implemented | Declares strata/bootstrap/backend/SDPD options. |
| `mrwin_config()` | implemented | Needs validation tests for invalid parameter lengths/ranges. |
| `mrwin_simulate()` | implemented | Simplified v5 DGP; needs parity fixtures. |
| `mrwin_kernel()` | implemented | Dense R backend; needs sparse/Rcpp backend. |
| `mrwin_stratum_win_loss()` | implemented | Needs more weighted edge-case tests. |
| `mrwin_stratum_log_cwr()` | implemented | Needs floor warning metadata. |
| `mrwin_prs_strata()` | implemented | Needs deterministic tie tests. |
| `mrwin_estimate()` | implemented | Point adjacent contrasts; no AL-CWR yet. |
| `mrwin_gls_pool()` | implemented | Uses approximate shrinkage; needs production LW review. |
| `mrwin_multiplier_bootstrap()` | implemented | No IPTW/bridging yet. |
| `mrwin_aalen_per_snp()` | implemented | Needs validation against survival-package reference or Python. |
| `mrwin_cox_per_snp()` | implemented | Approximate Cox; needs parity/reference tests. |
| `mrwin_mr_egger()` | implemented | Needs robust input validation. |
| `mrwin_pleiotropy_bounded_ci()` | implemented | Needs v5 Table 3 interpolation helper. |

### 5.2 Required High-Level API

To be added:

```r
fit <- mrwin(
  data,
  endpoint,
  genotype,
  exposure,
  beta_gwas,
  se_gwas,
  covariates = NULL,
  n_strata = 10,
  bootstrap = 500,
  adjustment = c("none", "ordinal_iptw", "gps"),
  backend = c("dense", "sparse", "rcpp")
)
```

S3 methods:

- `print.mrwin_fit()`.
- `summary.mrwin_fit()`.
- `plot.mrwin_fit()`.
- `tidy.mrwin_fit()` if broom compatibility is desired.
- `autoplot.mrwin_fit()` only if ggplot2 is added as Suggests.

High-level helper constructors:

- `mrwin_endpoint(time, status, priority = NULL)`.
- `mrwin_gwas(beta, se, snp = NULL, covariance = NULL)`.
- `mrwin_controls(...)` for bootstrap/backend/tolerance settings.

## 6. Result Schema

Every high-level fit must contain:

```text
class: c("mrwin_fit", "list")
call
data_info
endpoint_info
instrument_info
point
bootstrap
sdpd
diagnostics
warnings
session_info
```

Minimum fields:

- `point$delta_isg`.
- `point$log_theta`.
- `point$delta_x`.
- `point$delta_gls`.
- `point$dscwr`.
- `inference$se_delta_gls`.
- `inference$ci95_delta`.
- `inference$ci95_dscwr`.
- `inference$ci95_delta_fieller`.
- `heterogeneity$q`, `q_df`, `q_p_value`.
- `sdpd$intercept`, `intercept_p_value`, `scale`.
- `diagnostics$ess`, `weights`, `weak_instrument`, `moment_kurtosis`, `moment_skewness`.
- `warnings` as structured machine-readable codes.

Warning codes:

| Code | Trigger |
|---|---|
| `weak_instrument` | near-zero `Delta_X`, low F-stat, or unbounded Fieller CI. |
| `sdpd_rejected` | SDPD MR-Egger intercept p-value below alpha. |
| `sdpd_underpowered` | SDPD non-rejection but detectable gamma threshold too high. |
| `positivity_failure` | ESS threshold failure. |
| `bridged_strata` | any stratum excised and bridged. |
| `moment_condition_warning` | high kurtosis/skewness in bootstrap components. |
| `discordant_components` | component-specific MR directions conflict. |
| `ld_covariance_warning` | diagonal GWAS covariance used without LD evidence. |

## 7. Testing and Oracle Plan

### 7.1 Unit Tests

Required before any release:

- kernel hand examples.
- PRS strata deterministic ties.
- weighted win/loss identity.
- adjacent contrast calculations.
- Delta covariance manual matrix example.
- GLS diagonal covariance example.
- Fieller bounded/unbounded examples.
- MR-Egger exact intercept example.
- pleiotropy-bounded CI widening.

### 7.2 Python Parity Fixtures

Create fixed-seed fixtures under `inst/extdata/oracle/`:

| Fixture | Python source | R target |
|---|---|---|
| `kernel_small.rds/json` | `kernel.py` | exact equality. |
| `dgp_small.rds/json` | `dgp.py` | dimension and summary parity, not exact RNG parity unless generated externally. |
| `bootstrap_small.json` | `multiplier_bootstrap.py` | tolerance parity for `delta_GLS`, `se`, covariance dimensions. |
| `sdpd_small.json` | `inference.py` | tolerance parity for Aalen/Cox/MR-Egger. |
| `benchmark_small.json` | `benchmark_mr.py` | schema and directional parity. |

Use Python as the oracle only until R tests and independent statistical checks are stable.

### 7.3 Simulation Validation

Slow/integration tests:

- null calibration: Type-I error near 0.05.
- valid IV scenario: DS-CWR direction matches `alpha_X`.
- priority-1 pleiotropy: SDPD rejection increases with `gamma_1`.
- weak IV: Fieller can become unbounded and bivariate-Delta widens.
- discordant components: emit interpretation warning; cCWR should not overclaim a global direction.
- Q asymptotic: Gaussian parametric check approximates `chisq(D-2)`.

### 7.4 Data Validation Tests

Add tests for:

- missing time/status.
- non-binary status.
- negative times.
- mismatched dimensions.
- genotype/beta length mismatch.
- zero-variance SNP.
- empty stratum.
- zero/near-zero `Delta_X`.
- extreme weights and ESS failure.
- repeated or unordered priority labels.

## 8. Performance Plan

### 8.1 Backends

| Backend | Target |
|---|---|
| `dense_R` | small data, correctness, easiest debugging. |
| `sparse_R` | adjacent-stratum blocks without full `N x N` matrix. |
| `dense_cpp` | moderate data speedup. |
| `sparse_cpp` | biobank-scale operational default. |

Backend outputs must be identical within exact integer kernel equality and numeric tolerance for weighted sums.

### 8.2 Benchmarks

Benchmark grid:

- `N = 500, 1000, 5000, 10000`.
- `D = 5, 10`.
- `B = 50, 200`.
- dense vs sparse.

Required metrics:

- runtime.
- peak memory.
- valid bootstrap iterations.
- numerical equality of estimates.

Biobank-scale target:

- sparse backend must avoid full `N x N` materialization.
- output memory must scale with adjacent blocks and bootstrap summaries, not dense full matrix.

## 9. Documentation Plan

Required docs before full release:

- README minimal example.
- vignette 1: simulation quickstart.
- vignette 2: individual-level cohort template.
- vignette 3: diagnostics and interpretation.
- vignette 4: benchmark against per-component MR.
- article: caveats for non-collapsibility, weak instruments, pleiotropy, and discordant components.

User-facing summaries must avoid overclaiming:

- DS-CWR is a population-marginal local-gradient summary.
- It is not a structural per-unit conditional effect.
- It does not define a global transitive pairwise ranking.
- Component-discordant biology requires component-specific MR context.

## 10. Implementation Gap Register

| Gap | Priority | Blocking for full package? |
|---|---:|---|
| High-level `mrwin()` workflow | P0 | Yes |
| Data validation layer | P0 | Yes |
| Production IPTW with bootstrap refit | P0 | Yes |
| ESS and bridging | P0 | Yes |
| Sparse backend | P1 | Required for large data, not for first small release |
| Rcpp backend | P1 | Required for performance release |
| Exact/documented Ledoit-Wolf implementation | P0 | Yes |
| AL-CWR secondary diagnostic | P1 | Needed for paper completeness |
| v5 Table 3 bias interpolation | P1 | Needed for full pleiotropy-bounded workflow |
| Per-component benchmark estimators | P1 | Needed for Package D parity |
| Python oracle fixture generation | P0 | Yes |
| R CMD check CI | P0 | Yes |
| Vignettes and report output | P1 | Needed for user adoption |

## 11. Review Checklist

Statistical review:

- formulas match v5.2 Sections 2-5 and S7.
- bivariate extraction uses numerator/denominator covariance, not ratios.
- Q reference distribution is `chisq(D-2)`, not label permutation.
- SDPD interpretation follows S8/S9 caveats.
- discordant-component warning follows v5.2 Section 6.4/S13.

Software review:

- every exported function validates inputs.
- every warning is structured and testable.
- every random function has seed control.
- no generated manuscript outputs are committed.
- public result schema is stable.

Testing review:

- unit tests pass.
- Python parity fixtures pass.
- simulation smoke tests pass.
- R CMD check passes.
- coverage target is defined and tracked.

Release review:

- README example runs end-to-end.
- vignettes build.
- no hidden dependency on Python for normal R users.
- performance caveats are explicit for large data.

## 12. WP0 Decision Log

1. The v5.2 manuscript is the active source of truth; v5 is used only to understand baseline methods.
2. Package D is part of the full package roadmap because v5.2 made it a main-text/supplement requirement.
3. The R package must implement individual-level estimation first; summary-data mode is deferred.
4. Python remains the oracle for current reference behavior but must not be a runtime dependency for users.
5. The first production blocker is not more simulation code; it is data validation plus a stable high-level API.
6. Complete statistical correctness requires IPTW-with-bootstrap-refit and bridging; the current R implementation is intentionally incomplete until those are implemented and tested.

# Validation Findings — R2 (frequentist calibration)

Created: 2026-06-18. This records the first **external** validation of the
package (frequentist type-I error and CI coverage), run after the engineering
and analytic-variance work. It surfaced a real, quantified calibration problem
in the package's primary confidence interval that none of the internal
consistency checks (kernel parity, analytic==bootstrap) could detect.

## What was tested

- **Type-I error (null DGP):** `alpha_x=(0,0,0)`, `gamma_direct=(0,0,0)` — no
  causal effect and no pleiotropy, but the DGP still has confounding
  (`alpha_u=0.8`). A valid method must give a rejection rate ≈ 0.05.
- **Coverage (valid-IV DGP):** `alpha_x=(-0.4,-0.4,-0.4)`. The estimand "truth"
  was taken as the large-N (`N=20000`, averaged over 8 reps) pooled `delta_GLS`.
  A valid 95% CI must cover it ≈ 95% of the time.
- Inference: `inference="analytic"`, oracle GWAS weights (`sigma_beta=0`),
  `n_strata=5`. Both CIs the package reports — bivariate-Delta (primary) and
  Fieller (currently labelled a sensitivity diagnostic) — were compared.

## Results

Point estimator (null): **consistent and ~unbiased when the instrument is
adequate.** N=8000, 100 SNPs: median `delta_GLS` ≈ −0.02, sd 0.10; N=5000, 60
SNPs: median ≈ −0.01, sd 0.27. With a **weak** instrument (N=2000, 15 SNPs) it is
biased and extremely noisy (sd 1.37, 49% weak-instrument flag) — the well-known
fragility of IV ratio estimators (`logθ/ΔX`), already flagged by the package's
`weak_instrument` warning.

Inference calibration (N=6000, 80 SNPs, M=200 cohorts):

| Interval | type-I (null) | coverage (valid-IV) | verdict |
|---|---:|---:|---|
| **bivariate-Delta** (package primary) | **0.000** | **0.995** | mis-calibrated — ~2× too wide, no power |
| **Fieller** (labelled "sensitivity") | **0.055** | **0.955** | **correctly calibrated** |

Robust z-score diagnostic (bivariate-Delta): `median |delta/se| = 0.28`,
`q95 = 1.03`, `max = 1.57` (calibrated targets: 0.67, 1.96, →). Typical reported
`se` is ≈ 2× the empirical spread of the point estimate. The bivariate-Delta CI
can therefore **never** reach the rejection threshold → structural type-I = 0.

## Diagnosis

- The shared covariance `cov_u` (of `(logθ, ΔX)`) is **correct**: the Fieller CI
  uses `cov_u` directly and is calibrated. The analytic `cov_u` (validated
  `analytic == bootstrap`) is therefore vindicated.
- The bug is in the **bivariate-Delta ratio propagation** —
  `.mrwin_isg_covariance` (delta-method on `logθ/ΔX` via the gradient
  `(1/ΔX, −logθ/ΔX²)`) plus the GLS/Ledoit-Wolf pooling — which over-estimates
  the variance by ≈ 2×. The delta method is known to behave poorly for ratios;
  Fieller is the standard, correct construction for a ratio CI and here it works.
- This is a defect in the **original** WP4/WP5 inference design, not in the
  Phase-II fast/analytic work. It was never simulation-validated (the roadmap
  listed coverage as "deferred post-release"). Both the fast and analytic paths
  faithfully reproduced it.

## Honest lesson

Every prior "success" was **internal consistency** (fast == dense, analytic ==
bootstrap). The first **external** check (does the CI cover at 95%?) immediately
found a real bug. "analytic ≈ bootstrap" guarantees the two compute the same
number — not that the number is right; both reproduced the same mis-calibrated
interval.

## Recommended fix (validated)

**Make Fieller the primary interval** (it is already computed as
`ci95_delta_fieller`), or equivalently correct the bivariate-Delta variance.
Fieller is empirically calibrated (type-I 0.055, coverage 0.955). This must be
done **before M1**, because M1's inference builds on this layer.

Follow-ups: (a) identify and fix the ≈2× factor in the bivariate-Delta path so
both agree; (b) widen the validation grid (more N, D, instrument strengths,
`sigma_beta>0`, IPTW); (c) confirm Fieller handles weak instruments via its
unbounded-interval branch (it returned `unbounded=0` here, i.e. all bounded).

## Reproduce

`/tmp` scripts used: `r2_type1.R`, `r2_strong.R`, `r2_z.R`, `r2_fieller.R`,
`r2_cov_fieller.R` (to be folded into a committed WP19 harness
`R/validate_calibration.R`).

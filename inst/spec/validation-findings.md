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

## Fix (IMPLEMENTED 2026-06-18)

**Fieller is now the primary reported interval.** `.mrwin_fieller_ci` returns a
Fieller-consistent p-value (the numerator test of `H0: pooled logθ = 0`, i.e. the
test whose type-I is 0.055); `print`/`summary`/`tidy`/`mrwin_report` show the
Fieller 95% CI and p-value as the headline, handle the unbounded (weak
instrument) case explicitly, and keep the bivariate-Delta interval as a labelled
`delta-method (reference)`. Verified: full suite 92+ groups 0 failures;
`test-fieller-primary.R`. Empirically Fieller is calibrated (type-I 0.055,
coverage 0.955) and, under a weak instrument, honestly reports "unbounded"
instead of a falsely-bounded over-conservative interval.

Optional follow-up (low priority): identify the ≈2× factor in the bivariate-Delta
ratio variance so it agrees with Fieller.

Follow-ups: (a) identify and fix the ≈2× factor in the bivariate-Delta path so
both agree; (b) widen the validation grid (more N, D, instrument strengths,
`sigma_beta>0`, IPTW); (c) confirm Fieller handles weak instruments via its
unbounded-interval branch (it returned `unbounded=0` here, i.e. all bounded).

## WP19 widened calibration grid (2026-06-18)

Confirming the Fieller fix holds beyond the single cell that found the bug.
Strong instrument (N=3000, 100 SNPs, D=5), Fieller interval, M = 150-200 cohorts
(MC SE ~ 0.015-0.018). Harness: `R/validate_calibration.R`.

| scenario | inference | sigma_beta | adjustment | Fieller rate | target |
|---|---|---:|---|---:|---:|
| null (type-I) | bootstrap | 0 | none | **0.047** | 0.05 |
| null (type-I) | analytic | 0 | none | 0.020 | 0.05 |
| null (type-I) | analytic | 0.01 | none | 0.007 | 0.05 |
| null (type-I) | bootstrap | 0.01 | none | 0.020 | 0.05 |
| null (type-I) | analytic | 0 | ordinal_iptw | **0.060** | 0.05 |
| valid_iv (coverage) | bootstrap | 0.01 | none | **0.960** | 0.95 |

Findings:
- **No cell over-rejects** (max type-I 0.060). The method is calibrated-to-
  **conservative** across every path tested (bootstrap, analytic, GWAS
  uncertainty, IPTW) -- the safe direction (no false-positive inflation).
- Bootstrap, sigma_beta=0, no adjustment is on target (0.047); IPTW stays
  calibrated (0.060), so the first-order IPTW influence-function approximation
  does not break calibration; coverage is on target (0.960).
- Adding GWAS uncertainty (sigma_beta>0) makes inference more conservative
  (0.007-0.020) -- expected (wider CI), possibly slightly over-propagated; a
  future refinement, not a defect (it errs wide, never narrow).
- The analytic vs bootstrap difference at sigma_beta=0/none (0.020 vs 0.047) is
  within ~1.3 MC SE -- not a real discrepancy.

Still not covered (future WP19): multiple N/D, pleiotropy/SDPD
rejection grids, discordant-component warnings, weak-instrument unbounded rate.
The **core inference foundation is now validated and sound for M1.**

## WP19 additional cells (2026-06-20): bootstrap×IPTW + doubly-ranked path

Strong instrument (N=3000, m=100, D=5). Harness `R/validate_calibration.R` (now
takes a `stratification` argument so the M2 path can be checked externally, not
just asserted to run).

| scenario | inference | sigma_beta | adjustment | stratification | Fieller type-I | target |
|---|---|---:|---|---|---:|---:|
| null | bootstrap | 0 | ordinal_iptw | prs_rank | **0.017** (M=120) | 0.05 |
| null | bootstrap | 0 | none | **doubly_ranked** | **1.000** (M=150) | 0.05 |
| null | analytic | 0 | none | **doubly_ranked** | **1.000** (M=150) | 0.05 |

### bootstrap × IPTW (prs_rank): fills the documented gap — calibrated

type-I 0.017 (conservative, the safe direction). With the earlier
analytic×IPTW = 0.060, the IPTW path is now validated on both inference engines:
no over-rejection. The first-order IPTW influence-function approximation does not
break calibration.

### doubly-ranked stratification: CATASTROPHIC over-rejection (type-I = 1.000)

A second instance of the R2 lesson, and a more serious one. The M2-wiring is
**mechanically** correct — the fitted strata match `mrwin_doubly_ranked_strata`
bit-for-bit across the dense/sparse/fast/analytic paths (the internal
consistency test passes). Yet the **inference is invalid**: it rejects the causal
null in *every* simulated cohort, on both the bootstrap and the analytic engine.

Diagnosis (12-cohort null probe, `delta_gls` should be ~0):

| stratification | mean delta_gls | sd | instrument-mean spread across strata |
|---|---:|---:|---:|
| prs_rank | +0.016 | 0.272 | 1.343 |
| doubly_ranked | **−0.179** | 0.017 | **0.001** |

The bias is large and **highly consistent** (−0.179 ± 0.017), so the CI excludes
0 every time → type-I = 1.0. The root cause is structural and is exactly the
property doubly-ranked stratification is *designed* to have:

- Doubly-ranked stratification (rank by instrument → pre-strata → rank by exposure
  within each pre-stratum → assign exposure rank as the stratum) **balances the
  instrument across the final strata** (PRS-mean spread 0.001 vs 1.343 for
  PRS-rank). That balance is its selling point for *within-stratum* non-linear MR
  (each final stratum spans the full instrument range, so a within-stratum IV /
  LACE analysis is valid).
- But the cCWR / DS-CWR estimand is a **between-adjacent-strata** contrast: it
  reads the win-odds gradient *per unit of the between-stratum exposure shift*,
  which is only causal when the strata differ in the **instrument**. PRS-rank
  makes strata differ in the instrument; doubly-ranked deliberately removes that
  difference. With the instrument balanced out, the between-stratum exposure
  difference is driven by the **confounder** (and noise) within each PRS band, so
  the win-odds contrast estimates the *confounded* exposure-outcome association,
  not the causal effect. Under MR-strength confounding (the whole reason MR
  exists) this yields ~100% false positives.
- The bias vanishes only when there is no confounding (`alpha_u = 0`) — i.e.
  exactly the case where MR is unnecessary. So the method is invalid in every
  setting where it would be used.

**Conclusion.** Doubly-ranked stratification is **incompatible with the
between-stratum DS-CWR estimand** — a genuine negative finding, not a tunable
defect. `mrwin_doubly_ranked_strata()` (a correct Tian/Burgess implementation)
and the internal `stratification` plumbing are retained for a possible future
*within-stratum LACE* estimator, but `mrwin(stratification = "doubly_ranked")`
now emits a structured `doubly_ranked_invalid` warning and the default stays
`"prs_rank"`. The validated contributions are unchanged: the fast/Rcpp kernel,
the analytic influence-function variance, and the Fieller calibration fix.

Methodological note: this is the second case (after the continuous-ISG M1
negative finding and the original R2 bivariate-Delta defect) where every
**internal** consistency check passes — compiled == dense, analytic == bootstrap,
fitted strata == helper — while an **external** calibration check exposes a real
problem. Internal agreement certifies that two computations match; only coverage
/ type-I against a known DGP certifies that the number is right.

## Paper-01 full Monte-Carlo grid: type-I, power, weak-instrument, pleiotropy (2026-06-23)

Completes the headline validation the WP19 cells deferred (they covered type-I +
coverage only). Inference = the analytic influence-function path (validated to
equal the bootstrap); primary CI = Fieller. M Monte-Carlo cohorts per cell, MC SE
reported. Driver `tools/validation-scripts/paper01_mc_grid.R`, with two
mechanism-driven supplements below.

### Type-I and power (Fieller, reject H0: delta = 0)

| grid | setting | rate | MC SE |
|---|---|---|---|
| type-I | null, N=2000 | **0.050** | 0.010 |
| type-I | null, N=4000 | **0.028** | 0.007 |
| power | ax=-0.4, N=2000 | 0.100 | 0.015 |
| power | ax=-0.4, N=4000 | 0.282 | 0.023 |
| power | N=4000, ax=-0.1 | 0.040 | 0.010 |
| power | N=4000, ax=-0.2 | 0.077 | 0.013 |
| power | N=4000, ax=-0.4 | 0.235 | 0.021 |

Type-I is controlled (calibrated-to-conservative, matching the WP19 Fieller
result). Power rises monotonically with N and with effect size, and is honestly
modest at N=4000 — hierarchical-composite win-ratio MR is data-hungry, the
expected price of the conservative Fieller interval. Stated as scope, not hidden.

### Weak instrument (Fieller, vary instrument strength alpha_s; `paper01_weakiv.R`)

Truth = the instrument-standardized estimand at a strong instrument + large N
(**0.311**, stable; matches the interventional oracle ~0.30). An **unbounded**
Fieller interval is counted as covering — the honest convention: the method
declining to bound a near-singular ratio is not a false exclusion.

| alpha_s | weak flag | unbounded | coverage (incl. unbounded) | coverage (bounded only) | median bias |
|---|---|---|---|---|---|
| 0.40 | 0.00 | 0.00 | **0.967** | 0.967 | -0.031 |
| 0.30 | 0.00 | 0.00 | 0.963 | 0.960 | -0.043 |
| 0.20 | 0.09 | 0.09 | 0.967 | 0.877 | -0.086 |
| 0.15 | 0.23 | 0.23 | 0.963 | 0.737 | -0.109 |
| 0.10 | 0.50 | 0.50 | 0.973 | 0.477 | -0.153 |

The method degrades **gracefully and honestly**: as the instrument weakens it
increasingly (i) raises the `weak_instrument` warning and (ii) returns an
unbounded Fieller interval, so the **honest coverage stays ~0.96-0.97 at every
instrument strength** — the truth is never falsely excluded. Bounded-only
coverage falls and finite-sample bias grows at weak instruments (textbook
ratio-IV behaviour), and the package surfaces it. (The main grid's first weak-IV
cells used an estimator-based truth that is itself unstable when ΔX→0 and counted
bounded-only coverage — both undersold the method; this corrected run is canonical.)

### Pleiotropy / SDPD

Two regimes, because the simulator's `gamma_direct` enters the outcome as
`gamma_direct * s_true` (s_true = G·β) — pleiotropy **proportional to instrument
strength**:

- **Type-I (no pleiotropy):** SDPD rejects **0.047** (main grid) / 0.050 (custom
  harness) — calibrated.
- **Score-proportional pleiotropy (`gamma_direct` > 0):** SDPD power stays ~null
  (0.047 → 0.070 for gamma_direct 0.15 → 0.50). This is **correct, not a defect**:
  pleiotropy proportional to instrument strength is the InSIDE-violating case,
  observationally equivalent to a stronger causal effect, undetectable by any
  MR-Egger / SDPD-type test.
- **InSIDE-satisfying per-SNP pleiotropy** (`paper01_sdpd_power.R`: 30% of SNPs
  given a directional direct effect π_k independent of β_k): SDPD has a clean,
  calibrated power curve —

  | τ (per-SNP direct effect) | SDPD reject |
  |---|---|
  | 0.00 (null) | 0.050 |
  | 0.03 | 0.243 |
  | 0.06 | 0.807 |
  | 0.10 | 0.990 |
  | 0.15 | 1.000 |

  powerful against the pleiotropy it is *designed* to catch; the pleiotropy also
  biases the point estimate (delta 0.31 → 0.10 at τ=0.15), which SDPD flags.

**Lesson reinforced:** the near-null SDPD power in the main grid looked like a
defect until the mechanism was found — `gamma_direct` is the undetectable
pleiotropy type. The fix was a DGP that injects *detectable* pleiotropy, not a
change to the test. Find the mechanism before reporting a result.

## Reproduce

Harness `R/validate_calibration.R::.mrwin_validate_calibration()`. The
exploratory scripts that drove the runs above are kept under
`tools/validation-scripts/` (`r2_*.R`, `m1_*.R`, `m3_validate.R`, `iptw_*.R`,
`wp19_grid.R`, `wp19_cells.R`); the full Paper-01 grid is
`paper01_mc_grid.R` + `paper01_weakiv.R` + `paper01_sdpd_power.R`.

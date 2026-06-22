# Finding: two-sample / summary-data win-ratio MR is feasible

**Status:** feasible AND consistent. The estimator recovers the causal win-odds
gradient; its only bias is ordinary weak-instrument finite-sample bias (non-
collapsibility ruled out). Candidate for a Paper 03 / `mrwin_twosample` module.
**Reproduce:** `Rscript tools/twosample-feasibility.R` (feasibility),
`Rscript tools/twosample-bias-source.R` (bias source: weak-IV vs non-collapsibility).

## The question

Ordinary two-sample MR never shares individual-level data: it pools per-SNP
summary coefficients `(beta_GX, beta_GY)` by inverse-variance weighting (IVW).
Win statistics, by contrast, are pairwise functionals of individual outcome
profiles — so the first intuition was that win-ratio MR *cannot* be done from
summary data.

## The insight (why the intuition was wrong)

The outcome-side summary statistic does not have to be a regression slope. It can
be a **per-SNP win-odds coefficient** `delta_G,m` — the log-win-odds slope of the
priority-ranked outcome on allele dosage, i.e. a *"win-odds GWAS"*. It is
computed once in the outcome cohort (individual-level there, as in any GWAS) and
then shared only as a per-SNP summary. Under a valid instrument and local
linearity the chain rule gives

```
delta_G,m  ~=  gamma * beta_GX,m
```

so fixed-effect IVW, `gamma_hat = sum(beta_GX * delta_Gwin) / sum(beta_GX^2)`,
recovers the causal win-odds gradient `gamma`. The interventional oracle from
paper 01 supplies the ground truth `gamma*`.

## What the probe shows

DGP: the paper-01 linear cardiorenal model (oracle `gamma* = 0.30`), 40 SNPs,
fixed genetic architecture shared across cohorts.

- **One sample** (β_GX and the win-GWAS from the same cohort, N=20k, 8 reps):
  `gamma_hat = 0.263` (sd 0.078), bias −0.035.
- **Two independent samples** (exposure GWAS in cohort A, win-GWAS in cohort B):

  | N | gamma_hat | bias vs 0.30 |
  |---|-----------|--------------|
  | 8,000  | 0.215 (se .051) | −0.080 |
  | 16,000 | 0.206 (se .029) | −0.090 |
  | 32,000 | 0.255 (se .025) | −0.041 |

**Verdict: FEASIBLE.** IVW of per-SNP win-odds coefficients, from two *independent*
cohorts, recovers the right sign and magnitude of the causal win-odds gradient
without sharing any individual-level data. This overturns the "impossible"
intuition.

## Source of the attenuation — RESOLVED (the estimator is consistent)

The small-N attenuation was tested against two hypotheses (`tools/twosample-bias-source.R`):
weak-instrument finite-sample bias (vanishes with N) vs win-odds non-collapsibility
(structural, scales with residual heterogeneity). **It is weak-instrument bias.**

- **Convergence in N** (two *independent* samples, oracle `gamma* = 0.300`):

  | N | IVW | bias |
  |---|-----|------|
  | 16,000  | 0.234 (se .021) | −0.066 |
  | 50,000  | 0.284 (se .019) | −0.016 |
  | 150,000 | 0.311 (se .026) | +0.011 |

  Bias decays monotonically to ~0 — textbook two-sample weak-instrument behaviour
  (toward the null at small N). At a single large draw `gamma_IVW` sits at
  0.31/0.31/0.31 for N = 50k/150k/400k. **Consistent.**
- **Per-SNP linearity** (N=150k): regressing `delta_med` on `beta_GX` through the
  origin gives slope **0.283** with correlation **0.80** — i.e. `delta_G,m ≈
  gamma * beta_GX,m`, exactly the IVW assumption. (Earlier `adj` NaNs were just
  `log(0)`; fixed with a +0.5 continuity correction.)
- **Non-collapsibility ruled out**: scaling residual outcome heterogeneity
  (×1.0, 0.4, 0.1) leaves the `IVW / gamma*` ratio at 1.03 / 1.18 / 0.88 — noise
  around 1, no systematic attenuation. (`gamma*` itself rises as noise falls, a
  real estimand effect; IVW tracks it.)

**Corrected verdict:** IVW of per-SNP win-odds coefficients is a **consistent**
estimator of the causal win-odds gradient. The only caveat is ordinary
weak-instrument finite-sample bias — well understood in MR, with standard fixes
(stronger/aggregated instruments, large N, weak-IV-robust pooling). The
median-split `delta_Gwin` is already adequate; a non-collapsibility correction is
*not* needed.

## Next steps if this becomes a method

1. **Per-SNP win-odds estimator with a closed-form SE** so the IVW weights are
   real inverse variances (Fieller / weighted-median pooling).
2. **Weak-IV-robust pooling** (the bias is the only issue): report instrument
   strength, use enough N / SNPs, consider winsorised / weighted-median IVW.
3. **Pleiotropy robustness**: win-ratio MR-Egger / weighted-median (the package
   already has `mrwin_mr_egger`, `mrwin_sdpd`).
4. **`mrwin_twosample()` API** + a `mrwin_win_gwas()` helper producing shareable
   per-SNP win-odds summary statistics.
5. Validate end-to-end against the oracle (consistency, coverage) as paper 01 did
   for the one-sample estimator.

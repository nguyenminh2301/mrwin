# Finding: two-sample / summary-data win-ratio MR is feasible

**Status:** feasibility established (positive, with a caveat). Open question
flagged for a possible Paper 03 / `mrwin_twosample` module.
**Reproduce:** `Rscript tools/twosample-feasibility.R`

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

## The caveat (the open question)

There is a **modest downward attenuation (~15–30%)**. Across N the trend is noisy
at 6 reps (the 32k estimate, 0.255, is within ~1.6 se of `gamma*`), so part of it
is plausibly weak-instrument / finite-sample bias. But a residual structural
attenuation is likely, for a specific reason: the **win-odds is non-collapsible**,
so a per-SNP win-odds computed on a crude high/low dosage split is not exactly
`gamma * beta_GX,m`. The naive median-split `delta_Gwin` used here is a placeholder.

## Next steps if this becomes a method

1. **Principled per-SNP win-odds estimator** with a closed-form SE (so the IVW
   weights are real inverse variances, and Fieller/weighted-median pooling apply).
2. **Characterize the attenuation**: separate weak-instrument bias (vanishes with
   N / instrument strength) from non-collapsibility attenuation (structural);
   derive a collapsibility correction if needed.
3. **Pleiotropy robustness**: win-ratio analogues of MR-Egger / weighted-median
   (the package already has `mrwin_mr_egger`, `mrwin_sdpd`).
4. **Define and document a `mrwin_twosample()` API** + a `mrwin_win_gwas()` helper
   that produces shareable per-SNP win-odds summary statistics.
5. Validate end-to-end against the oracle (consistency, coverage) exactly as
   paper 01 did for the one-sample estimator.

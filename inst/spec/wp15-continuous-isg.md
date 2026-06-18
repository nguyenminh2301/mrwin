# WP15 — Continuous Kernel-Smoothed ISG (M1)

Status block: `T1 [todo] T2 [todo] T3 [todo] T4 [todo] T5 [todo]`
Branch: `C-wp15` (from `C`).
Depends on: WP13 (for the accelerated form; theory/prototype may use the dense
backend).
Blocks: WP19 (validation), feeds WP18 (theory) as a headline result.

The decile estimator forces an arbitrary stratum count `D`. The package already
flags "AL-CWR vs DS-CWR sensitivity to `D`" as an open gap
(`algorithm-spec.md` §4.7, `project-checkpoint.md`). This WP replaces the
discrete choice with a **continuous, bandwidth-controlled** gradient estimator
for which deciles are the boxcar special case.

---

## 1. The estimand

Let `r_i ∈ (0,1)` be the PRS-rank (or doubly-ranked position, WP16) of
individual `i`. Define the **local win-log-odds** and **local exposure** as
smooth functions of rank, and the continuous instrument-standardised gradient:

```
g(r)  = d/dr  E[ logit P(win | rank ≈ r) ]
x'(r) = d/dr  E[ X | rank ≈ r ]
δ(r)  = g(r) / x'(r)
DS-CWR_continuous = exp( ∫ w(r) δ(r) dr )     (w a normalised efficiency weight)
```

### 1.1 Estimation

Treat each ordered pair `(i,j)` with `r_i > r_j` as a Bernoulli "win" with
covariate equal to the rank gap `u = r_i − r_j` (and exposure gap `ΔX_{ij}`).
Fit a **local-linear / Nadaraya–Watson smoother** of the win indicator on `u`
near `u → 0⁺`; its slope at the origin estimates the local gradient. The
estimator is a **local U-statistic (U-process)**; the decile estimator is the
special case where `w` is a boxcar of width `1/D`.

Computationally this folds into the WP13 sweep: as the ordered sweep
accumulates win counts, bin/weight them by rank gap with kernel weights instead
of hard decile membership. No extra asymptotic cost over WP13.

### 1.2 Bandwidth selection

- Default: cross-validated bandwidth `b` minimising a leave-one-out win-loss
  prediction loss.
- Report a **sensitivity curve** `δ̂(b)` so users see stability vs the decile
  jumpiness.
- Provide `n_strata = Inf` / `estimator = "continuous"` as the trigger in
  `mrwin_controls()`.

---

## 2. Tasks

- **T1** — Theory note (to WP18): definition, identification, consistency under
  rank-preservation, U-process CLT sketch, boxcar-limit equivalence to the
  decile estimator. Must precede code sign-off but may be drafted in parallel.
- **T2** — Prototype `mrwin_continuous_isg(...)` on the dense backend; verify it
  reproduces the decile DS-CWR when the kernel is a boxcar of width `1/D`.
- **T3** — Accelerated implementation folded into the WP13 sweep (kernel-weighted
  win accumulation by rank gap).
- **T4** — Bandwidth CV + sensitivity-curve output; plotting via the existing
  `plot.mrwin_fit` (`type = "gradient"`).
- **T5** — Coverage simulation: under the WP8 DGP, the continuous estimator's
  CI achieves nominal coverage and lower `D`-sensitivity than the decile
  estimator.

---

## 3. Acceptance criteria

1. **Boxcar reduction**: with a boxcar kernel of width `1/D`, the continuous
   estimator equals the decile DS-CWR within `1e-8`.
2. **Coverage**: ≥ nominal (within MC error) on the null and valid-IV scenarios.
3. **Stability**: `δ̂(b)` sensitivity curve is smooth; reported. Decile
   estimator remains the default until coverage is validated across scenarios.

## 4. Definition of Done

§3 satisfied; status block all `[done]`; WP18 "Continuous ISG" section complete
with proof obligations discharged or explicitly flagged; checkpoint updated.

## 5. References

- Doubly-ranked stratification (Tian, Burgess, et al., 2023), *PLOS Genetics*
  19(6): e1010823 — rank-based strata without linearity/homogeneity (see WP16).
- Local U-statistics / U-processes: standard empirical-process theory; confirm
  exact citations in WP18.

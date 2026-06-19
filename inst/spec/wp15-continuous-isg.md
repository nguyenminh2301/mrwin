# WP15 — Continuous Kernel-Smoothed ISG (M1)

Status block: `EXPLORED -> NOT PURSUED (negative finding; see below)`. Reference
+ exact boxcar reduction done; consistency/sensitivity validation shows it does
not beat the decile. Decile + Fieller is the foundation.
Branch: developed on `C-wp13`.

## Design (2026-06-18) — concrete, validated estimand

`r_i = rank(PRS)/N`. At a centre `c` (quantile in (0,1)) with bandwidth `b` and
kernel `K`, define smooth half-open upper/lower membership weights (computed in
**integer rank space** `R = rank(PRS)` to avoid float boundary errors):

```
a_i(c) = K((R_i - cN)/(bN)) * 1{cN <  R_i <= cN + bN}      (upper)
b_i(c) = K((R_i - cN)/(bN)) * 1{cN-bN < R_i <= cN}          (lower)
W(c) = sum_ij a_i b_j 1{win};  L(c) = sum_ij a_i b_j 1{loss}
log_theta(c) = log(W/L);  Delta_X(c) = wmean(X|a) - wmean(X|b)
delta_isg(c) = log_theta(c)/Delta_X(c)
```

Pooled via the **validated calibrated inference**: the influence-function
covariance generalises to smooth weights —
`C[t,c] = a_t(P+_t/W0 - P-_t/L0) + b_t(Q+_t/W0 - Q-_t/L0)`,
`C[t,m+c] = a_t(X_t - Xbar_a)/sum(a) - b_t(X_t - Xbar_b)/sum(b)`,
`cov_u = C^T C` (cross-centre correlation via shared subjects) → `.mrwin_isg_covariance`
→ GLS + **Fieller** (the R2-validated, calibrated CI).

**Decile is the boxcar special case:** `K = boxcar`, `b = 1/D`, centres `c = k/D`
(k=1..D-1) reproduce the D-1 adjacent decile contrasts EXACTLY. Continuous = a
smooth kernel (Epanechnikov) + a centre grid + a continuous bandwidth `b`,
removing the arbitrary `D`.

Progress: `R/continuous_isg.R::mrwin_continuous_isg()` (internal reference, dense
kernel). **Exact boxcar reduction validated** vs the decile estimator
(`log_theta` diff 0, `cov_u` diff 0, `delta_gls` diff ~1e-14, Fieller match) over
N/D grid — `tests/testthat/test-continuous-isg.R`. The float-boundary bug found
in the first attempt (a single boundary subject mis-assigned) was fixed by
working in integer rank space — caught by the external decile reference, not
assumed (R2 lesson).

## NEGATIVE FINDING (2026-06-18) — M1 as designed does not deliver

External validation (the R2 lesson) shows the continuous estimator **does not**
achieve its headline goal and is **worse** than the decile:

- **Tuning sensitivity not reduced.** delta_gls varies as much across bandwidth
  `b` as the decile varies across `D` (CV ~1.9 vs ~2.1; at small `b` it is more
  unstable — ratio blow-up).
- **Higher variance.** At matched scale (decile D=5 ~ continuous b=0.10, N=8000,
  8 cohorts): decile sd 0.075; continuous GLS-pool sd 0.189; continuous simple
  Wald sd 0.278 — i.e. 2.5–3.7× noisier than the decile, regardless of pooling
  (the simple Wald is the noisiest, so pooling is NOT the cause; `sigma_isg`
  condition number ~2e3 is fine).
- **Root cause (the math).** The ISG is a RATIO `logθ/ΔX`. Smoothing into finer
  windows shrinks the denominator `ΔX` (upper- vs lower-half of a window are
  closer than adjacent disjoint decile bins), so the ratio gets NOISIER, not
  smoother. Coarse disjoint bins give larger, more stable `ΔX`. "Smoothing to
  remove D" trades the arbitrary-D problem for a worse variance problem.
- **The estimand is fine.** At large N the decile delta_gls is ~stable (~0.19,
  mild D-variation 0.14–0.24) — the cCWR is not fundamentally scale-broken; the
  problem is the continuous *estimator*, not the target.

**Decision:** the validated **decile estimator + Fieller inference** is the solid
foundation; M1 (continuous smoothing) is not worth pursuing as the headline. A
local-linear / single-regression-slope reformulation could be attempted, but the
root cause (ratio denominator shrinking at fine scale) would affect it too, so
expected payoff is low. The `D`-sensitivity concern is better handled by a
sensitivity analysis over `D` (report the range) than by a noisier continuous
estimator. `R/continuous_isg.R` is kept as a correct, boxcar-reduction-validated
reference of this exploration.

The project's genuine, validated contributions stand: the subquadratic + Rcpp
fast kernel (biobank-scale), the analytic influence-function variance, and the
R2 calibration finding + Fieller fix.

Original WP15 task list (T1..T5) follows below.
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

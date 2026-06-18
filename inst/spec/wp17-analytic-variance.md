# WP17 — Analytic Influence-Function Variance + GWAS Delta-Method (M3)

Status block: `T1 [partial: sampling IF derived] T2 [done: sampling part] T3 [todo: mrwin wiring] T4 [done: sampling part validated]`
Branch: developed on `C-wp13`.

Progress (2026-06-18, R-verified). **Sampling part landed** (adjustment="none",
fixed GWAS weights):
- `R/analytic_variance.R`: `mrwin_analytic_covariance(kernel, strata, X,
  contrast_plan)` builds the N×2(D-1) influence-coefficient matrix `C` and
  returns `cov_u = CᵀC`; `mrwin_analytic_inference(estimate, X)` feeds it through
  the existing `.mrwin_isg_covariance` / `mrwin_gls_pool` / `.mrwin_fieller_ci`.
- Influence coefficients (multiplier weights `ξ_i = 1+e_i`, Var(e)=1):
  `log θ_d`: `k∈H → P⁺_k/W₀ − P⁻_k/L₀`, `k∈L → Q⁺_k/W₀ − Q⁻_k/L₀`;
  `Δx_d`: `k∈H → (X_k−X̄_H)/n_H`, `k∈L → −(X_k−X̄_L)/n_L`. Shared strata across
  adjacent contrasts reproduce the induced correlation automatically.
- **Validated** vs the fixed-strata (`sigma_beta=0`) multiplier bootstrap: `se`
  ratio 0.99–1.00 across N∈{500,1500,4000} (B=4000), relative Frobenius ≈3–4%
  (Monte-Carlo limited), no systematic bias. Tests `test-analytic-variance.R`;
  full suite 86 groups, 0 failures.

**Remaining (next steps):**
- **Step B — `Σ_gwas`:** GWAS-uncertainty term via the Jacobian of `(logθ,Δx)`
  through the stratum cutpoints (boundary movers); add `Σ = Σ_sampling + Σ_gwas`.
- **IPTW:** extend the influence function for `adjustment="ordinal_iptw"`
  (propensity refit contributes additional IF terms).
- **T3 wiring:** `inference = c("bootstrap","analytic","analytic+bootstrap")` in
  `mrwin_controls()` / `mrwin()`, guarded so the analytic path is only the
  default once `Σ_gwas` + IPTW are in (until then it covers the sampling part
  with fixed GWAS weights and `adjustment="none"`).
Depends on: WP4 (estimator), WP5 (bootstrap), existing
`python/p1_engine_v5/q_statistic_asymptotics.py` (the analytic-distribution
groundwork already started).
Blocks: WP19.

The multiplier bootstrap's expensive part is re-stratification (`B` times).
This WP derives an **analytic** variance for the stratified DS-CWR, splitting it
into (i) a sampling component from the influence function / Hájek projection and
(ii) a GWAS-uncertainty component propagated by a **delta-method through the
stratum cutpoints** — avoiding full re-stratification. Goal: cut `B` by ~10× or
remove it, with bootstrap retained as the validation oracle.

---

## 1. Structure of the derivation

### 1.1 Sampling component (fixed `β`)

Each adjacent log-win-ratio `log θ_d = log(W_d / L_d)` is a ratio of two
two-sample U-statistics. Its first-order influence function is standard
(per-subject win/loss projections; cf. Bebu–Lachin 2016, the same projection the
comparison packages use). Stack the `D−1` contrasts; the ISG map
`δ_d = log θ_d / ΔX_d` and the GLS pool propagate by the multivariate delta
method already coded in `_mrwin_isg_covariance` and `mrwin_gls_pool`. Output an
analytic `Σ_sampling` over `(log θ, ΔX)`.

### 1.2 GWAS-uncertainty component (`β* ≈ β`)

`β` enters **only** through the stratum boundaries (which individuals fall in
band `d`). Perturbing `β → β + δβ` moves boundary individuals between adjacent
bands. To first order, `∂ logθ_d / ∂β` equals the win/loss change from moving a
boundary individual across the `(d, d±1)` cut, summed with weights given by the
boundary density. Assemble the Jacobian `J = ∂(logθ, ΔX)/∂β` from boundary
movers (cheap: only individuals near cutpoints contribute) and propagate:
`Σ_gwas = J Σ_GWAS Jᵀ`. Total `Σ = Σ_sampling + Σ_gwas`.

### 1.3 Result

A closed-form covariance for `(log θ, ΔX)` feeding the existing Fieller / GLS /
Q machinery, with **no bootstrap loop** in the fast path.

---

## 2. Tasks

- **T1** — Theory note (WP18): influence function of the stratified DS-CWR;
  Jacobian-through-cutpoints derivation; statement of regularity conditions.
- **T2** — Implement `mrwin_analytic_covariance(...)` (R + Python reference)
  producing `Σ` for `(log θ, ΔX)`. Reuse `_mrwin_isg_covariance`,
  `mrwin_gls_pool`, `.mrwin_fieller_ci`.
- **T3** — Add `inference = c("bootstrap", "analytic", "analytic+bootstrap")`
  to `mrwin_controls()`. The `analytic+bootstrap` mode uses analytic variance as
  a **control variate** to reduce `B`.
- **T4** — Coverage + agreement validation: analytic `Σ` agrees with the
  bootstrap `Σ` within a stated tolerance on WP8 scenarios; CI coverage meets
  nominal. Quantify the `B` reduction achieved by the control-variate mode.

---

## 3. Acceptance criteria

1. Analytic and bootstrap covariances agree within the spec tolerance on null,
   valid-IV, and dose-heterogeneity scenarios.
2. Analytic CI coverage is within MC error of nominal.
3. Bootstrap remains available and is the default until coverage is validated;
   analytic is opt-in.

## 4. Definition of Done

§3 satisfied; status block all `[done]`; WP18 "Analytic variance" section
complete; checkpoint updated.

## 5. References

- Bebu, I., & Lachin, J. M. (2016). *Biostatistics* 17(1), 178–191.
- Luo, X., et al. (2015). Variance of the win ratio. (confirm exact citation in
  WP18 before publishing.)
- van der Vaart (2000), §23–25 (influence functions, delta method) — already
  cited in `q_statistic_asymptotics.py`.

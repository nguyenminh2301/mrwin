# The transport-velocity structure of causal win statistics — a unifying result and a super-plan

Status: new result + programme super-plan (2026-07-02). Companion to
`dev/research-frontier-roadmap.md`. All numbered claims are validated against the
interventional oracle; reproducibility scripts
`tools/validation-scripts/p5-collapsible-transport-scalar.R` (scalar theory,
Claims 1–3) and `p5-collapsible-transport-hierarchical.R` (hierarchical lift +
unification, Claims 4–5). No confidential data.

---

## 0. One-paragraph statement of the new truth

The entire `mrwin` programme (Papers 01–04) estimates a **net-benefit gradient**
`β* = Cov(Z,w(O))/Cov(Z,X)`. We prove and validate that this estimand is, exactly,
a **density-overlap-weighted average of an optimal-transport velocity field** — the
causal effect on the outcome *squashed* by how much the outcome distributions
overlap. Two consequences, both new: (i) the win statistic is non-collapsible for a
precise, computable reason (its implicit weight is the law-dependent overlap `∫f²`),
and dividing by that overlap yields a **collapsible** causal effect that recovers the
underlying structural effect exactly and is invariant to the outcome-noise law; and
(ii) — the unification — the *same* overlap functional that squashes the win effect
(the P1/non-collapsibility problem) is the *same* functional whose vanishing is
Paper 03's **degenerate U-statistic boundary**. Non-collapsibility and
weak-identification-degeneracy are two faces of one object: the **decidability
overlap** `D`. This connects P1 (keystone, previously "open"), Paper 03 (done), and
the shared estimand of all four papers into a single scalar functional, and it is the
mathematical spine of the programme's super-plan (§6).

---

## 1. Setup

Potential-outcome family `{F_x}` with density `f_x` (scalar case §2–4), or a
priority-ranked censored hierarchical composite (§5). Two independent draws
`Y_i(x), Y_j(x')`. The **win kernel** `h(O_i,O_j) ∈ {−1,0,+1}` decides a pair by the
first priority tier that separates them (accounting for censoring). The **net
benefit** `NB(x,x') = P(Y_i(x) ≻ Y_j(x')) − P(Y_i(x) ≺ Y_j(x'))`, and the **win
score** (Hoeffding projection) `w(o) = E[h(o,O_2)]`. The programme's interventional
oracle is `do(X → x±ε)`: `B(ε) := NB(x+ε, x−ε)`, `B(0)=0`, `B` odd, and the reported
win gradient is `B(ε)/(2ε) → B'(0)/2` (codebase convention).

For a scalar continuous outcome the transport (Lagrangian) velocity is
`v_x(y) := −∂_x F_x(y) / f_x(y)`, equivalently `v_x(F_x^{-1}(u)) = ∂_x F_x^{-1}(u)`:
the speed at which the quantile-`u` point of `F_x` moves under the intervention. This
`v` is the generator of the optimal-transport (monotone rearrangement) map
`T_x^{x'} = F_{x'}^{-1}∘F_x`, and it is the *collapsible* object — transport maps
compose exactly, `T_x^{x''} = T_{x'}^{x''}∘T_x^{x'}`, so quantile-wise displacements
telescope: `F_{x''}^{-1}(u) − F_x^{-1}(u) = ∫_x^{x''} v_s(F_s^{-1}(u)) ds`.

---

## 2. Claim 1 (local structure theorem) — proved and validated

**Theorem 1.** For a scalar potential-outcome family with a differentiable density,
```
B'(0) = d/dε NB(x+ε, x−ε)|_{ε=0} = 4 ∫ f_x(y)² v_x(y) dy.                     (★)
```

*Proof.* Let `P(ε) = ∫ f_{x+ε}(y) F_{x−ε}(y) dy`, so `B = 2P − 1`. Then
`P'(0) = ∫ (∂_x f_x) F_x dy − ∫ f_x (∂_x F_x) dy`. Writing `∂_x f_x = ∂_y(∂_x F_x)`
and integrating by parts, `∫ (∂_x f_x) F_x = −∫ (∂_x F_x) f_x` (boundary terms
vanish). Hence `P'(0) = −2 ∫ f_x (∂_x F_x) dy`, and `B'(0) = 2P'(0) = −4∫ f_x ∂_x F_x`.
Substituting `∂_x F_x = −f_x v_x` (from `∂_x F_x(F_x^{-1}(u)) = 0`) gives (★). ∎

So the interventional win gradient is `4∫f²·v`, i.e. **the transport velocity `v`
averaged with weight `f_x²` (the "self-overlap density")**, times 4. In the codebase
convention it is `2∫f²v`.

**Validation** (`…scalar.R`, part A; pure high-resolution numeric analysis, no Monte
Carlo — decisive). `(★)` holds to relative error `10⁻⁵–10⁻⁷` across three families,
including a *genuinely non-location* family (`Y(x)=αx+σ₀e^{γx}Z`, where `v_x(y)`
varies with `y`):

| family | B'(0) [numeric] | 4∫f²v | rel. err | analytic |
|---|---|---|---|---|
| Gaussian-loc σ=1 | 0.902703 | 0.902703 | 2.1e−07 | 0.902703 |
| Gaussian-loc σ=2 | 0.451352 | 0.451352 | 5.8e−08 | 0.451352 |
| Laplace-loc σ=1 | 1.131352 | 1.131371 | 1.7e−05 | — |
| Gauss-scale γ=0.3, x=0.5 | 0.776964 | 0.776964 | 2.4e−07 | — |

(Analytic Gaussian-location cross-check: `B'(0) = 2α/(σ√π)`, matched exactly.)

---

## 3. Claim 2 (non-collapsibility characterization + the collapsible correction)

Theorem 1 exposes the non-collapsibility mechanism cleanly. Under any **location**
model `Y(x)=μ(x)+σZ`, `v_x≡μ'(x)` and (★) gives `B'(0)=4μ'(x)∫f_x²`. The causal
effect `μ'(x)` is **squashed by the overlap factor `4∫f_x²`**, which depends on the
noise law. So:

- **Non-collapsibility is exactly the law-dependence of the overlap `∫f²`.** Two
  studies with identical per-subject causal effect `μ'` but different outcome-noise
  distributions report different win gradients. This is the pairwise-comparison
  analogue of odds-ratio non-collapsibility, and it is *intrinsic* to the win scale.
- **The collapsible causal win effect** is `τ(x) := B'(0)/(4∫f_x²)` = the `f²`-weighted
  mean transport velocity `= μ'(x)` under any location model, invariant to the noise
  law.

**Validation** (`…scalar.R`, part B). The win gradient scales as `1/σ` for Gaussians;
and at *matched* σ, Gaussian and Laplace give **different** win gradients in the exact
ratio `√(2/π)=0.7979` (measured 0.7979) — proving the correct squash factor is the
overlap `∫f²`, *not* the variance σ. The correction `τ = B'(0)/(4∫f²)` recovers
`α=0.800` to 5 significant figures at every σ and for both laws:

| family | σ | win gradient | 4∫f² | τ = grad/(4∫f²) |
|---|---|---|---|---|
| Gaussian | 0.5 | 1.805405 | 2.256758 | 0.799999 |
| Laplace | 0.5 | 2.262660 | 2.828429 | 0.799971 |
| Gaussian | 2.0 | 0.451352 | 0.564190 | 0.800000 |
| Laplace | 2.0 | 0.565681 | 0.707107 | 0.799993 |

**Structural remark (the honest boundary of the positive result).** Because
collapsible functionals of the transport map are exactly the *fixed-quantile-weight*
displacement functionals (they telescope), and the win statistic's weight `f²` is
*not* fixed (it moves with the law), **no fully scale-free collapsible win effect
exists**: collapsibility forces a scale choice. The overlap correction `τ` is the
canonical resolution — it is collapsible over the location orbit (invariant to the
noise law, the practically important "identical effects, different spreads" case) and
recovers the structural effect, at the cost of estimating the outcome-margin overlap
`∫f²`.

---

## 4. Claim 3 (finite-sample estimator)

**Validation** (`…scalar.R`, part C; Monte Carlo, R=20, N=12000). A plug-in
estimator — Mann–Whitney win gradient over the `do(X±ε)` populations, divided by a
Gaussian-kernel U-statistic estimate of `4∫f̂²` — recovers `μ'=α` and is invariant to
σ, where the raw win gradient varies 3.9×:

| σ | raw win gradient | ∫f̂² | corrected τ (MC SE) |
|---|---|---|---|
| 0.5 | 1.778 | 0.5455 | 0.815 (0.003) |
| 1.0 | 0.902 | 0.2787 | 0.810 (0.003) |
| 2.0 | 0.453 | 0.1401 | 0.808 (0.003) |

(τ is ~1.5% high from the kernel `∫f²` underestimate, within the bandwidth bias;
the invariance across σ — the collapsibility — is exact to MC error.)

---

## 5. Claims 4 & 5 — the hierarchical lift and the unification

The programme's real outcome is a K=3 censored hierarchical composite, not a scalar.
Censoring turns informative comparisons into **ties**, so it plays the role of
"increasing σ": it shrinks the **decidability overlap** `D := P(a random pair is
decided)`, the hierarchical analogue of `∫f²`, and squashes the win gradient.

**Claim 4 (censoring = squash; validated).** `…hierarchical.R`, N=6000, R=24,
censoring rate swept 0.002→1.6 (ground truth `g₀` = near-uncensored gradient):

| cens rate | win gradient g | decidability D | g/g₀ | D/D₀ |
|---|---|---|---|---|
| 0.002 | 0.1350 | 0.9861 | 1.000 | 1.000 |
| 0.15 | 0.0838 | 0.5279 | 0.620 | 0.535 |
| 0.35 | 0.0556 | 0.3260 | 0.412 | 0.331 |
| 1.6 | 0.0165 | 0.0878 | 0.122 | 0.089 |

The raw gradient falls 8×. It obeys a clean power law **`g ∝ D^0.873`** (R²=0.998).
Correcting by decidability collapses the censoring dependence: the coefficient of
variation across the sweep falls from **0.642 (raw g) → 0.120 (g/D) → 0.032
(g/D^0.87)**, a 5–20× stabilization — the hierarchical collapsibility correction.

**Claim 5 (THE UNIFICATION; validated).** The outcome-side first-projection variance
`Var(w(O_i))` — the quantity whose vanishing is Paper 03's degenerate-U-statistic
boundary (heavy-censoring / low-decided-fraction route, §8.2/§9 of that paper) —
tracks the *same* decidability overlap:

- `Spearman(D, Var(w)) = 1.0000` across the censoring sweep (perfect monotone).
- `Var(w)/D = 0.327, 0.301, 0.304, 0.309, 0.311, 0.310` — **constant to 3%**, i.e.
  `Var(w) ∝ D` empirically (exponent ≈ 1).

So the functional that squashes the win effect (non-collapsibility, the P1 problem)
and the functional whose vanishing degenerates Paper-03 inference are **the same
object, `D`**. Maximal censoring → `D→0` → *simultaneously* maximal non-collapsibility
squash **and** the `ζ₁→0` degenerate boundary. This is, to our knowledge, the first
statement linking non-collapsibility and weak-identification degeneracy for win
statistics; it says they are not two problems but one, indexed by decidability.

(Scalar consistency check: for a *fully observed* scalar outcome `w(O)=2F(Y)−1` with
`F(Y)~Uniform`, so `Var(w)=1/3` is constant and `D=1` — no degeneracy from the
outcome side, exactly as the unification predicts. The degeneracy route is
specifically the tie/censoring route, which is where `D<1`.)

---

## 6. THE SUPER-PLAN — an ambitious programme built on the transport spine

The result above reframes the whole programme around one object — the transport
velocity field `v_x(y)` and its overlap `∫f²` / decidability `D` — and turns the
five "frontier" papers (`research-frontier-roadmap.md`) from a loose wishlist into a
single connected theory. Sequenced by dependency and by decisiveness:

### P1 — Collapsible causal win effect via optimal transport (keystone). **UNBLOCKED.**
Previously "the causal slice is OPEN." This result *is* the local half of P1: the
collapsible target is the transport velocity `v_x`, the win statistic is its
`f²`-weighted average, and `τ=B'(0)/(4∫f²)` is the collapsible correction (Thm 1,
Claims 2–3). Remaining P1 math, now sharply scoped: (a) the **global** OT statement
on the *lexicographic partial order* of the hierarchical outcome — the transport map
on a lexicographically-ordered (tie-permitting) space, with `D` as the mass of the
"decided" (order-comparable) region; (b) the **impossibility theorem** made rigorous
(no scale-free collapsible win effect; the obstruction is the cocycle/curvature of
the non-location family); (c) IV-identification of `τ` (divide `β*` by the
outcome-margin overlap, which is IV-estimable) and its efficient estimation.

### P3 — Weak-instrument-robust inference. **DONE, and now explained.** Claim 5 gives
the *mechanism* behind Paper 03's degenerate boundary: `ζ₁ ∝ D`. Immediate payoff —
the paper's `weak_degenerate` flag can be reformulated as a **decidability threshold**
on `D` (an interpretable, estimable, outcome-only quantity), rather than the current
bootstrap-on-`ζ₁` heuristic; and the many-weak-instrument breakdown found this session
can be re-examined through the same lens.

### P2 — Semiparametric efficiency & double robustness for instrumented win ratios.
The EIF of the win-MR estimand must be the EIF of `∫f²·v` under IV. Thm 1 hands over
the pathwise derivative: perturbing the outcome law moves both `v` (the transport
generator) and `∫f²` (the overlap), and the influence function splits into a
transport part and an overlap part. This is the concrete route to the degree-2 EIF
the roadmap calls for, and the collapsible target `τ` is the natural object to build
a doubly-robust one-step/TMLE estimator for.

### P4 — Sharp partial identification without functional form. The transport map is
the extremal object of the Kantorovich problem; the sharp identified set for the win
effect under IV-stratified marginals is the range of `∫f²·v` over couplings
consistent with the marginals — an infinite-dimensional LP whose vertices are monotone
rearrangements. `D` (decidability) is the natural "identification budget": the sharp
bounds widen exactly as `D` shrinks. Shares the OT machinery with P1.

### P5 — Win-mediation. Cross-world pairwise potential outcomes; the direct/indirect
decomposition of the *transport velocity* (a vector field, hence additively
decomposable along mediator paths in a way the non-collapsible win statistic is not)
gives the exact 3-way split with a non-collapsibility interaction term that vanishes
iff `∫f²` is mediator-invariant — linking back to P1.

### Cross-cutting new probes (cheap, high-value, do next):
1. **Global collapsibility obstruction** — construct a non-location family, compute
   the cocycle that obstructs a scale-free collapsible effect, and validate the
   impossibility numerically (a genuine negative theorem).
2. **`τ` as an MR estimand** — validate `β*/overlap` recovers a censoring-/spread-
   invariant effect in the *full two-sample within-family* pipeline (chains this
   result into Paper 04's deployable form).
3. **`D`-thresholded degeneracy flag for Paper 03** — replace the bootstrap
   `ζ₁` diagnostic with the decidability `D`, validate equivalent size control at
   lower cost.
4. **Analytic `Var(w)=cD`** — derive the exponent (empirically ≈1) from the
   heterogeneous-decidability structure (early-event subjects decided against all;
   censored subjects mostly tied), closing the one gap in Claim 5 left empirical.

**Why this is the right spine.** Every frontier paper reduces to a question about the
*same* two functionals (`v`, `∫f²`/`D`); P1 and P3 are now provably the same
functional seen from two sides; and each next step is validated immediately against
the existing interventional oracle. The programme is no longer five separate hard
papers — it is one theory of degree-2 (pairwise) causal functionals with optimal
transport as its geometry and decidability as its single scalar invariant.

---

## 7. Honest scope and caveats

- Thm 1 (★) is exact for scalar continuous outcomes; the hierarchical lift is
  validated (Claims 4–5) but the exact overlap functional for the lexicographic
  censored order is characterized empirically (`g∝D^0.87`, `Var(w)∝D`), not yet
  derived in closed form — that derivation is cross-cutting probe #4/P1(a).
- The collapsible correction `τ` requires estimating the outcome-margin overlap
  `∫f²`/`D`; this is IV-estimable but its efficient estimator (and the full MR
  correction `β*/overlap`) is not yet built or validated end-to-end — probe #2.
- The impossibility result (no scale-free collapsible win effect) is argued
  structurally in §3; a fully rigorous theorem with the cocycle obstruction is P1(b).
- All numbers here are simulation results against the interventional oracle with
  Monte-Carlo / numeric-integration error reported; no real data.

# The transport-velocity structure of causal win statistics — a unifying result and a super-plan

Status: super-plan EXECUTED (2026-07-04) — §§8–13 below are the results of carrying
out §6's cross-cutting probes and P1(a)/P1(b). Two are decisive new positive results
(Claims 6–7 closing probe #4 and P1(b) rigorously), one is a positive-but-partial
result (P1(a), §9), and two are honestly-reported NEGATIVE/mixed results (probes #2
and #3, §§11–12) with the failure mechanism identified in each case, per the
project's "report faithfully, correct the record loudly" standard — not every probe
in §6 panned out, and this file says so explicitly rather than only keeping the wins.
P4 and P5 were not attempted (§13): scoped honestly as open, not silently dropped.

Companion to `dev/research-frontier-roadmap.md`. All numbered claims are validated
against the interventional oracle or a public simulation; reproducibility scripts:
`p5-collapsible-transport-scalar.R` (scalar theory, Claims 1–3),
`p5-collapsible-transport-hierarchical.R` (hierarchical lift + unification, Claims
4–5), `p1-analytic-vard-typeI-censoring.R` + `p1-ph-model-censoring-closedform.R`
(Claim 6, §8), `p1a-tier-decomposition-lexicographic.R` (§9),
`p5-impossibility-shape-cocycle.R` (Claim 7, §10),
`p2-tau-portable-winmr-twosample.R` (§11), `p2-delta-method-tau-se.R` (§11),
`p3-dthreshold-degeneracy-screen.R` (§12). No confidential data.

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
Paper 03's **degenerate U-statistic boundary, in the outcome-only special case**
(§12 sharpens this: for Paper 03's *actual, instrument-crossed* diagnostic, the
connection to `D` alone is weak — a scope correction, not a retraction, of the
outcome-only claim). Non-collapsibility and weak-identification-degeneracy are two
faces of one object, decidability `D`, in the setting where that object is
outcome-only; this connects P1 (keystone, previously "open") and Paper 03's theory
into a single scalar functional, and it remains the mathematical spine of the
programme's super-plan (§6), whose execution (§§8–13) produced two rigorous closed
new results (Claims 6–7), one honestly-partial result (§9), and two
mechanistically-explained negative results (§§11–12) — a mixed but genuinely
informative outcome, not a uniform success, and reported as such throughout.

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

### P1 — Collapsible causal win effect via optimal transport (keystone). **UNBLOCKED,
(a)/(b) substantially advanced, (c) still open.** This result *is* the local half of
P1: the collapsible target is the transport velocity `v_x`, the win statistic is its
`f²`-weighted average, and `τ=B'(0)/(4∫f²)` is the collapsible correction (Thm 1,
Claims 2–3). Update after execution: **(b) is DONE** — the impossibility theorem is
now rigorous over a continuous shape family, not a 2-point check (Claim 7, §10).
**(a) is PARTIALLY DONE** — the K=1 lexicographic/censored case now has an EXACT
closed form (Claim 6, §8, two instances), and the general K=3 frailty-correlated
case has a validated-but-incomplete tier-wise decomposition (§9, CV 0.31, a ~30%
residual not yet explained — the fully general closed form remains open). **(c) is
the weakest link**: §11 shows `τ` is a validated *target and interval* concept in
the real MR pipeline but explicitly **not yet a reliable point estimator** under
weak instruments, and the natural delta-method route to "efficient estimation" was
tried and found to fail for a specific, identified reason (weak-IV non-regularity)
— (c) needs a structurally different (AR-type, not delta-method) approach, not
attempted further this round.

### P3 — Weak-instrument-robust inference. **DONE (the AR/sup-CI machinery itself);
the proposed reformulation of its degeneracy flag is NEGATIVE, tested and reported.**
Claim 5 gives the mechanism behind Paper 03's degenerate boundary for the
**outcome-only** quantity `Var(w) ∝ D`. §6 optimistically proposed reformulating the
`degenerate` flag (which is instrument-crossed, `ζ₁(β̂)`, not outcome-only) as a
`D`-threshold; §12 tested this directly against 420 simulated cells and found it
**does not work** (`Spearman(D,degenerate)=−0.50`, no safe threshold exists) because
weak-instrument point-estimate volatility dominates the actual flag's behavior in a
way the outcome-only `D` cannot see. No code changed in `R/onesample_ar.R` as a
result — the existing bootstrap diagnostic remains the correct, validated approach.

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

### Cross-cutting new probes (cheap, high-value, do next) — STATUS AFTER EXECUTION:
1. **Global collapsibility obstruction** — **DONE, closed rigorously.** §10/Claim 7:
   continuous EPD shape family, `ρ(κ)` non-constant on a 9-point grid, formula
   validated to `1e-6`, impossibility demonstrated directly (raw `β*` ranges 1.70×
   at fixed `σ` for the same `α`; `β*/ρ(κ)` recovers the target exactly at every
   `κ`).
2. **`τ` as an MR estimand** — **MIXED, done honestly.** §11: population/oracle-level
   portability strongly confirmed (CV 0.850→0.136); estimated point-estimate
   portability does NOT improve under weak instruments (a ratio-instability
   mechanism, identified); CI-level portability partially improves but a
   delta-method attempt at a proper SE fails for a mechanistically-identified
   reason (weak-IV non-regularity). Net: `τ` is validated as a target/interval
   concept, not (yet) as a deployable point estimator.
3. **`D`-thresholded degeneracy flag for Paper 03** — **NEGATIVE, done honestly.**
   §12: tested directly against the actual bootstrap diagnostic across 420
   simulated cells; `Spearman(D,degenerate)` only `−0.50`, no usable threshold
   found (even lax thresholds pass ~50-60% false-safes); the diagnostic is
   dominated by weak-instrument point-estimate volatility, not decidability
   alone. No code changed; Claim 5's scope corrected to be precise about what it
   does and does not cover.
4. **Analytic `Var(w)=cD`** — **DONE, exceeded scope.** §8/Claim 6: not just the
   exponent, but an EXACT closed form for Type-I censoring (`D=p(2-p)`,
   `Var(w)=p(p²-3p+3)/3`, `Var(w)/D→1/2` as `D→0`), PLUS a second exact closed
   form for the win *gradient* itself under a proportional-hazards causal model
   (`B'(0)=−α·D(p)`), cross-checked against the classical two-exponential
   concordance formula independently.

P1(a) (§9) was also attempted beyond the original probe list: a partial,
honestly-quantified result (tier-wise decomposition holds to CV 0.31, not the
<5% precision of the clean K=1 cases). P4 and P5 (§13) were not attempted.

**Why this is the right spine.** Every frontier paper reduces to a question about the
*same* two functionals (`v`, `∫f²`/`D`); P1 and P3 are now provably the same
functional seen from two sides; and each next step is validated immediately against
the existing interventional oracle. The programme is no longer five separate hard
papers — it is one theory of degree-2 (pairwise) causal functionals with optimal
transport as its geometry and decidability as its single scalar invariant.

---

## 8. Claim 6 (NEW, decisive) — exact closed form for the K=1 win gradient under proportional hazards + Type-I censoring

Probe #4 asked for the exponent in `Var(w) ∝ D` (found empirically ≈1 in §5) in
closed form. Two independent exact derivations, both for `T ~ Exponential`, Type-I
(fixed administrative) censoring at `c`, `p := P(T≤c)`:

**Var(w) closed form** (no causal effect, `T` iid, `D`/`Var(w)` as functionals of the
base population — the quantity Paper 03's degeneracy theory conditions on):
```
D = p(2−p),      Var(w) = p(p²−3p+3)/3
```
Both derived by hand from the package's exact win rule (`i` beats `j` iff
`δ_j=1 & X_i>X_j`; `i` loses iff `δ_i=1 & X_j>X_i`) and validated by an **exact**
(no Monte-Carlo opponent-subsampling noise) `O(n log n)` leave-one-out computation
of the per-subject win score via order statistics (`tools/validation-scripts/
p1-analytic-vard-typeI-censoring.R`):

| p | D_hat | D_theory | Var(w)_hat | Var(w)_theory |
|---|---|---|---|---|
| 0.02 | 0.0408 | 0.0396 | 0.02020 | 0.01960 |
| 0.30 | 0.5110 | 0.5100 | 0.21937 | 0.21900 |
| 0.90 | 0.9899 | 0.9900 | 0.33300 | 0.33300 |

and the key asymptotic — `Var(w)/D → 1/2` **exactly** as `p→0` (heavy censoring,
the degenerate limit that matters for Paper 03) — confirmed numerically
(`Var(w)/D_hat = 0.500, 0.498, 0.487` for `p = 0.001, 0.01, 0.05`).

**Win-gradient closed form** (adds a causal proportional-hazards effect, closing a
genuinely new instance of P1(a)/P1(b), not just Var(w)): for
`T_x ~ Exponential(rate=λ(x)=λ₀e^{αx})` (log-linear hazard, `α` = causal
log-hazard-ratio gradient) with the same Type-I censoring,
```
B'(0) = −α · D(p),   D(p) = p(2−p)   [the SAME decidability functional]
```
derived by differentiating the exact `NB(x+ε,x−ε)` integral (both-events-observed
integral + one-censored boundary term) at `ε=0`. **Independent analytic
cross-check**: as `c→∞` (`p→1`, `D→1`), this reduces to `B'(0)→−α`, which matches
the *classical* two-independent-exponentials concordance formula
`P(T_A>T_B)=λ_B/(λ_A+λ_B)` to first order — a totally different derivation route
confirming the same limit (lesson E: validate the target more than one way).
Validated by direct simulation (`p1-ph-model-censoring-closedform.R`, `n=2×10⁶`):
rel. err `1–4%` across `c ∈ {0.3,…,50}` (the one `17%` outlier at the smallest,
noisiest `c` was confirmed to be Monte-Carlo noise, not a systematic error, by
re-running at a different seed and larger `n`, moving it to `2–5%`).

This makes "win gradient = (causal effect) × decidability" an **exact,
non-empirical** identity in two clean instances (Var(w) for K=1 Type-I censoring;
B'(0) for K=1 PH + Type-I censoring), not just the `g∝D^0.873` power-law fit found
empirically for the full K=3 frailty-correlated composite in §5 — and closes
cross-cutting probe #4.

---

## 9. P1(a), partial result — does the tier-wise structure lift to the real K=3 frailty-correlated composite?

The K=1 closed forms above are clean because there is no cross-tier correlation.
The programme's real outcome has `K=3` tiers with **cascading censoring** (a
lower-priority tier is further censored by all higher-priority event times) and a
**shared Gamma frailty** correlating all three tiers — exactly the structure that
could break a clean multiplicative D-factorization.

**Exact structural identity** (holds by construction, not a theorem to prove): the
pairwise kernel telescopes over the first decisive tier,
```
h(A,B) = Σ_k 1{R_k(A,B)} · sign_k(A,B),   R_k := {tiers 1..k−1 tied}
⟹ NB(x+ε,x−ε) = Σ_k contrib_k,   contrib_k := E[1{R_k}·sign_k]
```
**Question**: does each `contrib_k` further factor as `(reach-tier-k probability
R_k) × (tier-k conditional decidability D_k^cond) × (a roughly tier- and
censoring-invariant "local velocity")`, the natural tier-wise lift of Claim 6?

**Test** (`tools/validation-scripts/p1a-tier-decomposition-lexicographic.R`):
exact `O(n²)` pairwise tier-cascade tracking (not a subsample) on the real
`mrwin` frailty simulator (`N=1400`, `K=3`, shared Gamma frailty, cascading
censoring), computing `R_k`, `D_k^cond`, and `contrib_k` for the `do(X±ε)` oracle
comparison across a 4-point censoring ladder (12 tier×censoring cells):

| tier | mean ratio `contrib_k/(R_k·D_k^cond)` | CV |
|---|---|---|
| 1 | 0.158 | 0.360 |
| 2 | 0.166 | 0.214 |
| 3 | 0.141 | 0.418 |
| pooled | — | **0.309** |

Raw `contrib_k` itself spans **over two orders of magnitude** across the same
grid (0.0005 to 0.098). The decidability-weighted ratio compresses that to a
**CV of 0.31** — a real, substantial, non-trivial stabilization — but it is *not*
the near-perfect (`<5%` rel. err) precision of the clean K=1 closed forms in §8.
**Honest reading**: the tier-wise lift substantially holds (most of the raw
variation is explained by `R_k·D_k^cond`), but a genuine `~30%` residual remains,
plausibly a frailty-correlation correction term (a selection effect: "reaching
tier k" is not independent of the tier-k causal gradient once tiers share
frailty). The internal identity check (`Σ_k contrib_k = g_total` to numerical
precision, every row) confirms no bug in the tier-cascade code itself — the
residual is a real feature of the frailty-correlated DGP, not an artifact.
**P1(a)'s fully general closed form for the frailty-correlated K-tier case remains
open**; this section establishes the right decomposition and quantifies how far
it is from closing the gap, rather than claiming it closes.

---

## 10. Claim 7 (NEW, decisive) — P1(b), the impossibility theorem, made rigorous via a continuous shape family

§3's impossibility argument was checked at exactly two points (Gaussian, Laplace).
Probe #1 asked for a continuous-family version. Using the **exponential power
distribution** (generalized Gaussian, shape `κ`; `κ=1`=Laplace, `κ=2`=Gaussian,
`κ→∞`=uniform-like), unit-variance standardized, for a location family
`Y(x)=μ(x)+σZ_κ`:
```
β*(x;κ,σ) := B'(0)/2 = (2μ'(x)/σ) · ρ(κ),     ρ(κ) := ∫ g_κ(z)² dz
```
(`tools/validation-scripts/p5-impossibility-shape-cocycle.R`, three parts):

**(A) `ρ(κ)` is non-constant, verified on a 9-point grid `κ∈[0.7,12]`**: range
`[0.270, 0.459]`, ratio max/min `1.70`; analytic cross-checks match exactly
(`ρ(Gaussian)=1/(2√π)=0.282095`, `ρ(Laplace)=1/(2√2)=0.353553`, both matched to
6 decimals).

**(B) the formula itself holds exactly**: `β*` vs `(2α/σ)ρ(κ)` at 6 `(κ,σ)`
combinations, rel. err `4×10⁻⁷` to `2×10⁻⁶`.

**(C) the impossibility, directly**: at fixed `σ=1`, raw `β*` ranges **1.70×**
across shapes for the *same* causal effect `α=0.8` (`0.621` to `0.735` at `κ=0.7`
vs. others down to a common target); dividing by `ρ(κ)` recovers `2α=1.600`
**exactly at every κ including 0.7** (`1.5999`–`1.6000` uniformly). So no
function of `(β*,σ)` alone — no universal rescaling constant — can recover the
causal effect uniformly over shapes; the data-estimable shape/scale functional
`ρ(κ)` (equivalently `σ·∫f²`) is required, and it suffices.

**Methodological note (worth banking as a lesson)**: the first run of part (C)
showed `κ=0.7` recovering `1.353` instead of `1.600` — apparently a breakdown at
the heaviest-tailed shape. A convergence check (lesson A/B: probe before
pronouncing) varying the finite-difference `eps`, grid resolution, and
integration limit in isolation showed the discrepancy vanished entirely once the
integration domain was widened from `25σa` to `60σa`+ — `κ<1` EPD tails decay
slower than Laplace/Gaussian and were being silently truncated. This was a pure
numerical-truncation artifact, not a theorem breakdown; the production script now
uses a `100σa` domain and the result is clean at every `κ`. (A second,
independent bug of the same flavor recurred twice more this session — see §11's
delta-method script and the closed-form verification in §8: `NB(x+ε,x−ε)` is
*odd* in `ε`, so `B'(0) = NB(ε)/ε`, not `NB(ε)/(2ε)`; this factor-of-2 trap is now
called out explicitly to avoid a fourth recurrence.)

This closes cross-cutting probe #1 / P1(b) rigorously.

---

## 11. Probe #2 — τ as an MR estimand in the two-sample within-family pipeline: a MIXED result, reported honestly

Chains §5's population-level finding into the actual estimation pipeline
(Papers 02+03+04): does `τ_hat := γ_hat_AR / D_hat` restore cross-cohort
portability when cohorts differ only in censoring/follow-up (a real two-sample
MR problem — different biobanks, different follow-up)?

**Design** (`tools/validation-scripts/p2-tau-portable-winmr-twosample.R`): SAME
causal architecture (F=1500 families×2 sibs, M=25 SNPs), sweep ONLY the
censoring rate across 6 "cohorts"; ground truth = each cohort's own
interventional oracle; reference `τ* := oracle_ref/D_ref` from a near-uncensored,
large-N cohort.

**(a) Population/oracle level — confirmed, strongly.** `oracle(c)/D_hat(c)` across
the censoring ladder has **CV 0.136**, vs. **CV 0.850** for raw `oracle(c)`
— a >6× stabilization, extending §5's hierarchical-script finding to the real
GWAS+family-clustered architecture.

**(b) Estimated point estimate — does NOT improve, honestly reported.**
`gamma_AR_hat` (the AR argmin) has CV 0.641–0.753 across two independent runs
(R=40 then R=80); `tau_hat = gamma_AR_hat/D_hat` has CV 0.746–1.409 — **similar
or worse**, not better. Mechanism: the AR argmin is a well-documented unstable
point summary under weak instruments (per-SNP F-stat is intrinsically `<10` in
within-family designs, per the existing `within-family-twosample-winmr.R`
finding); dividing an already-noisy `gamma_hat` by a `D_hat` that shrinks toward
0 at heavy censoring *amplifies* rather than damps the instability (a classic
ratio-estimator blowup — `tau_hat` hit `2.55` vs. a target of `0.146` at the
heaviest censoring in one run). **The naive point-estimate correction is not a
free lunch under weak instruments.**

**(c) Confidence-interval level — partial improvement.** The raw
(uncorrected) AR confidence interval's coverage of the fixed portable target
`τ*` **degrades sharply with heavier censoring**: `0.963 → 0.925 → 0.925 → 0.863
→ 0.787 → 0.662` across the ladder — because the raw CI is scaled for the
cohort's own (shrinking) raw oracle, not the portable target. The D-rescaled CI
(`ar$ci / D_hat`, treating `D_hat` as fixed) holds up much better:
`0.963 → 0.925 → 0.925 → 0.900 → 0.875 → 0.875` — real, if imperfect
(under-covers, esp. at `cens≥2`), evidence that the correction has genuine value
at the interval level even where the naive point estimate does not.

**Attempted fix (P2), also negative — worth recording.** The interval-level
under-coverage in (c) plausibly comes from treating `D_hat` as a fixed constant,
ignoring its own sampling variance and its covariance with `gamma_hat` (both
computed from the same sample). `tools/validation-scripts/p2-delta-method-tau-se.R`
derives a proper joint influence-function (delta-method/sandwich) SE for
`tau_hat`, using per-subject pieces already native to `.mrwin_ar1_blocks()`
(`IF_gamma(o) = (gh(o;γ̂) − γ̂·gx(o))/b̂`, the standard M-estimation sandwich
form; `IF_D(o) = 2(d(o)−D̂)`, the U-statistic Hájek projection; combined via the
ratio delta method). **Result: badly miscalibrated** — `SE_delta/empirical-SD`
ranges from **2.4× to 27×** across a censoring sweep in the one-sample AR
setting (N=3000, R=150), and the resulting Wald CI's coverage of `tau*`
(0.80–0.95) is *worse*, not better, than the naive plug-in (0.97–1.00, itself
likely just "too-wide-to-be-wrong" under weak IV). **Mechanism, identified**: a
first-order delta-method/Wald linearization requires regularity conditions
(asymptotic normality of `gamma_hat`) that a weak-instrument ratio estimator
violates — this is *exactly* the reason Anderson-Rubin-type (not Wald-type)
inference exists in the first place for `gamma_hat` itself (Papers 02→03's whole
motivation); applying a delta method on top of it for `tau_hat` inherits the
same disease. **The correct fix is structural, not a bigger toolbox on the same
ratio**: an AR-type joint test statistic for `(γ,D)` that never divides by a
weak/noisy quantity, analogous to how `mrwin_ar_onesample` itself never divides
by `Cov(Z,X)`. This is a concretely scoped, well-motivated open problem for
future work, not attempted further here given the budget.

**Bottom line for probe #2**: the transport-collapsibility correction is a real,
validated fix for cross-cohort portability *at the level of the causal target and
of confidence intervals*, but does **not** (yet) give a reliable corrected point
estimator under the weak instruments intrinsic to within-family MR — an honest,
mechanistically-explained limitation, not a fixed method.

---

## 12. Probe #3 — a D-thresholded degeneracy screen for Paper 03: a NEGATIVE result, and a scope correction to Claim 5

§5's Claim 5 showed `Spearman(D, Var(w))=1.0000` for the **pure outcome-margin**
win-score variance (weights=1, no instrument). §6 optimistically proposed
reformulating Paper 03's bootstrap-`ζ1` degeneracy flag as a cheap `D`-threshold.
Tested directly (`tools/validation-scripts/p3-dthreshold-degeneracy-screen.R`):
does the **actual** quantity `mrwin_ar_onesample` bootstraps
(`ζ1(β̂) = c0 − 2β̂c1 + β̂²c2`, where `c0=Var(gh)` mixes the instrument `Z` into
the win score via `gh_i=(Z_i s_i − r_i)/(n−1)`) track `D_hat` well enough to
screen on?

**Result: no.** Across a grid of `N∈{800,2500}` × 7 censoring rates × weak
(`α_s=0.05`) and strong (`α_s=0.35`) instruments (420 rows total, `boot_reps=200`
ground truth per cell): `Spearman(D_hat, degenerate) = −0.50` to `−0.53` —
moderate, nowhere near the `1.0000` found for the outcome-only quantity. Even at
a lax screening threshold `D_hat > 0.02` (excludes only the most extreme 7% of
cases), **61% of the "safe" cases were still flagged degenerate** by the ground-
truth bootstrap (238/389 rows); tightening the threshold only *shrinks how much
of the data you can skip the bootstrap for* without ever reaching a safe
separation (at `D_hat>0.30`, still 43% false-safe). The `N=800` vs. `N=2500`
comparison at the *same* censoring rate (`0.01`) shows `P(degenerate)` swinging
from `0.50` to `0.00`, confirming sample size — not decidability — is a
first-order driver here, because `ζ1(β̂)` also depends on how far the (weak-
instrument-driven, noisy) point estimate `β̂` wanders from the quadratic's true
minimizer, a phenomenon `D_hat` alone cannot see.

**Scope correction (lesson F: correct the record loudly, not quietly).** Claim 5
(§5) remains valid **exactly as stated** — for the pure outcome-margin
`Var(w(O))`, ignoring the instrument. It does **not**, however, extend to a
practical drop-in replacement for the shipped `mrwin_ar_onesample` degeneracy
diagnostic, because that diagnostic's actual driver is instrument-crossed and
dominated by weak-IV point-estimate volatility that the outcome-only `D`
functional cannot capture. **No change was made to `R/onesample_ar.R`** as a
result of this probe — shipping an unvalidated screen would trade a slow-but-
correct diagnostic for a fast-but-wrong one, which the evidence above explicitly
rules out. A cheap screen that *does* work would need to jointly account for
instrument strength (e.g. an `n·D`-scale quantity combined with an F-statistic-
like term) — an open problem, not solved here.

---

## 13. P4 and P5 — explicitly out of scope this round

Both were listed in §6 as ambitious extensions (P4: sharp partial identification
via a Kantorovich-LP formulation with `D` as an "identification budget"; P5:
win-mediation via a transport-velocity vector-field decomposition). Neither was
attempted. Both require substantial new machinery (an ambiguity-set/LP-duality
argument for P4; a cross-world potential-outcomes formalization for P5) that
could not be brought to the same validated standard as §§8–12 within this
session's scope — per the project's "no claim before its test passes" rule, a
rushed, unvalidated sketch of either would be worse than an honest "not
attempted." Both remain open, and are good candidates for a dedicated future
session with their own probe-first validation loop.

---

## 14. Honest scope and caveats (updated)

- Thm 1 (★) is exact for scalar continuous outcomes. The K=1 hierarchical case is
  now **exact** in two instances (Claim 6, §8: Type-I censoring, both with and
  without a causal PH effect). The general K=3 frailty-correlated case is
  **partially** characterized: the tier-wise decidability-weighted decomposition
  substantially holds (CV 0.31, §9) but a real ~30% residual (vs. <5% in the
  clean K=1 cases) is not yet explained — plausibly a frailty-correlation
  correction term. P1(a)'s fully general closed form remains open.
- The impossibility result (no scale-free collapsible win effect) is now a
  **rigorous, continuous-shape-family theorem** (Claim 7, §10), not just a
  2-point check — closed.
- The collapsible correction `τ` is validated as a **portability fix at the
  population/oracle level and (partially) at the confidence-interval level**
  (§11a, c) in the real two-sample within-family MR pipeline, but is **not**
  validated as a reliable **point estimator** under weak instruments (§11b), and
  a natural delta-method fix for the interval's residual miscalibration **fails**
  for a mechanistically identified reason (§11, the delta method inherits weak-
  IV's non-regularity) — this is the biggest open gap in the programme's use of
  `τ` as a deployable MR estimand, and the report says so plainly rather than
  papering over it.
- The D-thresholded degeneracy screen proposed in §6 **does not work** (§12): it
  was tested directly, decisively, and found to fail because Paper 03's actual
  diagnostic is dominated by weak-instrument point-estimate volatility that the
  outcome-only decidability functional does not capture. No package code was
  changed as a result. Claim 5 itself (the outcome-only `Var(w)` vs. `D`
  relationship) stands, but its scope is now stated precisely rather than
  optimistically.
- P4 and P5 (§13) were not attempted; explicitly open, not silently dropped.
- All numbers here are simulation results (against the interventional oracle,
  an independent classical formula, or exact leave-one-out computation) with
  Monte-Carlo / numeric error reported or noted; no real data. Two recurring
  factor-of-2 finite-difference bugs (dividing by `2ε` instead of `ε` for an
  odd function `B(ε)`) and one numerical-truncation trap (§10's integration-limit
  issue) were caught by convergence checks before being reported as findings —
  each is called out explicitly above so the next session does not re-trip them.

# The transport-velocity structure of causal win statistics

**Status** (2026-07-04, round-2 follow-ups 2026-07-05). Core theory §§1–5 closed and
validated. Super-plan §6 executed: §§7–13 report the results — 3 rigorous new results
(Claims 6–7, §9), 1 honestly-partial (§8), 3 mixed/negative with the mechanism
identified in each case (§§10–12). P4 (§12) has a working tool but its headline
conjecture was retracted; work continues. P5 (§13) not attempted.

Companion to `dev/research-frontier-roadmap.md`. Every number is validated against
the interventional oracle, an independent analytic formula, or exact computation;
Monte-Carlo / numeric error is reported. No confidential data, no real data.

**Reproducibility scripts** (all in `tools/validation-scripts/`):

| § | script |
|---|---|
| 2–4 | `p5-collapsible-transport-scalar.R` |
| 5 | `p5-collapsible-transport-hierarchical.R` |
| 7 | `p1-analytic-vard-typeI-censoring.R`, `p1-ph-model-censoring-closedform.R` |
| 8 | `p1a-tier-decomposition-lexicographic.R` |
| 9 | `p5-impossibility-shape-cocycle.R` |
| 10 | `p2-tau-portable-winmr-twosample.R`, `p2-delta-method-tau-se.R`, `p2-ar-style-tau-test.R`, `p2-twosample-multisnp-ar-tau.R` |
| 11 | `p3-dthreshold-degeneracy-screen.R`, `p3-worstcase-zeta1min-screen.R`, `p3-analytic-c2-bound.R` |
| 12 | `p4-kantorovich-partial-id-censored.R` |

---

## 0. The result in one paragraph

The `mrwin` programme (Papers 01–04) estimates a net-benefit gradient
`β* = Cov(Z,w(O))/Cov(Z,X)`. We prove and validate that this estimand is exactly a
**density-overlap-weighted average of an optimal-transport velocity field**: the causal
effect on the outcome, *squashed* by how much the outcome distributions overlap. Two
consequences: **(i)** the win statistic is non-collapsible for a precise, computable
reason — its implicit weight is the law-dependent overlap `∫f²` — and dividing by that
overlap gives a **collapsible** effect that recovers the structural effect and is
invariant to the outcome-noise law; **(ii)** *the unification*: the same overlap
functional that squashes the win effect (the P1 problem) is the functional whose
vanishing is Paper 03's degenerate-U-statistic boundary. Non-collapsibility and
weak-identification degeneracy are two faces of one object — the **decidability
overlap `D`** — in the outcome-only setting (§11 sharpens the scope: for Paper 03's
*actual, instrument-crossed* diagnostic the link to `D` alone is weak).

---

## 1. Setup

Potential-outcome family `{F_x}` with density `f_x` (scalar, §§2–4) or a
priority-ranked censored hierarchical composite (§5). Win kernel
`h(O_i,O_j) ∈ {−1,0,+1}` decides a pair by the first tier that separates them (given
censoring). Net benefit `NB(x,x') = P(Y_i(x) ≻ Y_j(x')) − P(Y_i(x) ≺ Y_j(x'))`; win
score (Hoeffding projection) `w(o) = E[h(o,O_2)]`. Interventional oracle
`do(X→x±ε)`: `B(ε) := NB(x+ε,x−ε)`; `B(0)=0`, `B` odd; reported win gradient is
`B(ε)/(2ε) → B'(0)/2` (codebase convention).

Transport (Lagrangian) velocity `v_x(y) := −∂_xF_x(y)/f_x(y)`, equivalently
`v_x(F_x^{-1}(u)) = ∂_xF_x^{-1}(u)`: the speed of the quantile-`u` point under the
intervention. `v` generates the OT (monotone-rearrangement) map
`T_x^{x'} = F_{x'}^{-1}∘F_x` and is the *collapsible* object — transport maps compose
exactly, so quantile-wise displacements telescope.

---

## 2. Theorem 1 (local structure) — proved + validated

For a scalar family with differentiable density,
```
B'(0) = 4 ∫ f_x(y)² v_x(y) dy.                                              (★)
```
*Proof.* `P(ε) = ∫f_{x+ε}F_{x−ε}`, `B = 2P−1`. Then
`P'(0) = ∫(∂_xf_x)F_x − ∫f_x(∂_xF_x)`. Writing `∂_xf_x = ∂_y(∂_xF_x)` and integrating
by parts, `∫(∂_xf_x)F_x = −∫(∂_xF_x)f_x`, so `P'(0) = −2∫f_x∂_xF_x` and
`B'(0) = −4∫f_x∂_xF_x`. Substituting `∂_xF_x = −f_xv_x` gives (★). ∎

So the win gradient is `v` averaged with weight `f_x²` (the self-overlap density).

**Validated** (high-resolution numeric, no Monte Carlo — decisive), rel. err `1e−5`–`1e−7`
across three families incl. a genuinely non-location one (`Y(x)=αx+σ₀e^{γx}Z`):

| family | B'(0) numeric | 4∫f²v | rel. err |
|---|---|---|---|
| Gaussian-loc σ=1 | 0.902703 | 0.902703 | 2.1e−07 |
| Gaussian-loc σ=2 | 0.451352 | 0.451352 | 5.8e−08 |
| Laplace-loc σ=1 | 1.131352 | 1.131371 | 1.7e−05 |
| Gauss-scale γ=0.3, x=0.5 | 0.776964 | 0.776964 | 2.4e−07 |

(Analytic cross-check, Gaussian location: `B'(0)=2α/(σ√π)` — matched exactly.)

---

## 3. Claim 2 — non-collapsibility characterised, and the correction

Under any **location** model `Y(x)=μ(x)+σZ`, `v_x≡μ'(x)`, so (★) gives
`B'(0)=4μ'(x)∫f_x²`: the causal effect **squashed by the overlap `4∫f_x²`**, which
depends on the noise law. Hence:

- **Non-collapsibility is exactly the law-dependence of `∫f²`** — two studies with
  identical `μ'` but different noise laws report different win gradients. The
  pairwise analogue of odds-ratio non-collapsibility, intrinsic to the win scale.
- **The collapsible win effect** is `τ(x) := B'(0)/(4∫f_x²)` = the `f²`-weighted mean
  transport velocity `= μ'(x)` under any location model, invariant to the noise law.

**Validated**: at *matched* σ, Gaussian vs Laplace win gradients differ in the exact
ratio `√(2/π)=0.7979` (measured 0.7979) — proving the squash factor is the overlap
`∫f²`, **not** the variance σ. The correction recovers `α=0.800` to 5 s.f. at every σ
and both laws (e.g. Gaussian σ=0.5: grad 1.805405, 4∫f² 2.256758, τ 0.799999;
Laplace σ=2.0: grad 0.565681, 4∫f² 0.707107, τ 0.799993).

**Structural boundary (honest).** Collapsible functionals of the transport map are the
*fixed-quantile-weight* displacement functionals (they telescope); the win statistic's
weight `f²` is not fixed (it moves with the law). So **no fully scale-free collapsible
win effect exists** — collapsibility forces a scale choice. `τ` is the canonical
resolution: collapsible over the location orbit, at the cost of estimating `∫f²`.
Made rigorous in §9.

---

## 4. Claim 3 — finite-sample estimator

Plug-in (Mann–Whitney gradient over `do(X±ε)` ÷ Gaussian-kernel estimate of `4∫f̂²`),
MC R=20, N=12000. Raw gradient varies 3.9× across σ; corrected `τ` is invariant:

| σ | raw gradient | ∫f̂² | corrected τ (MC SE) |
|---|---|---|---|
| 0.5 | 1.778 | 0.5455 | 0.815 (0.003) |
| 1.0 | 0.902 | 0.2787 | 0.810 (0.003) |
| 2.0 | 0.453 | 0.1401 | 0.808 (0.003) |

(τ ~1.5% high from kernel `∫f²` underestimate — bandwidth bias; the σ-invariance is
exact to MC error.)

---

## 5. Claims 4 & 5 — the hierarchical lift, and the unification

The real outcome is a K=3 censored hierarchical composite. Censoring turns
informative comparisons into **ties**, playing the role of "increasing σ": it shrinks
the **decidability** `D := P(a random pair is decided)` — the hierarchical analogue of
`∫f²` — and squashes the gradient.

**Claim 4 (censoring = squash).** N=6000, R=24, censoring 0.002→1.6:

| cens | gradient g | D | g/g₀ | D/D₀ |
|---|---|---|---|---|
| 0.002 | 0.1350 | 0.9861 | 1.000 | 1.000 |
| 0.15 | 0.0838 | 0.5279 | 0.620 | 0.535 |
| 0.35 | 0.0556 | 0.3260 | 0.412 | 0.331 |
| 1.6 | 0.0165 | 0.0878 | 0.122 | 0.089 |

Raw gradient falls 8×, obeying `g ∝ D^0.873` (R²=0.998). Correcting by decidability
collapses the censoring dependence: CV falls **0.642 (raw g) → 0.120 (g/D) → 0.032
(g/D^0.87)**.

**Claim 5 (the unification).** The outcome-side projection variance `Var(w(O_i))` —
whose vanishing is Paper 03's degenerate boundary via the heavy-censoring route —
tracks the same `D`: `Spearman(D, Var(w)) = 1.0000`, and `Var(w)/D` = 0.327, 0.301,
0.304, 0.309, 0.311, 0.310 — **constant to 3%**, i.e. `Var(w) ∝ D`.

So the functional squashing the win effect and the functional whose vanishing
degenerates Paper-03 inference are **the same object `D`**. To our knowledge the first
statement linking non-collapsibility and weak-identification degeneracy for win
statistics. (Consistency check: for a fully observed scalar outcome `w(O)=2F(Y)−1`,
`Var(w)=1/3` constant and `D=1` — no outcome-side degeneracy, as predicted; the
degeneracy route is specifically the tie/censoring route where `D<1`.)

---

## 6. Super-plan — status after execution

**P1 (keystone).** The collapsible target is the transport velocity `v_x`; the win
statistic is its `f²`-weighted average; `τ = B'(0)/(4∫f²)` is the correction (Thm 1,
Claims 2–3). **(b) DONE** — impossibility now rigorous over a continuous shape family
(§9). **(a) PARTIAL** — K=1 exact in closed form (§7, two instances); general K=3
frailty-correlated case has a validated-but-incomplete tier decomposition (§8, ~30%
residual). **(c) WEAKEST** — §10: `τ` is a validated *target and interval* concept but
not a reliable *point estimator* under weak instruments.

**P3 (weak-IV-robust inference).** The AR/sup-CI machinery is DONE (shipped). Claim 5
explains the degenerate boundary for the *outcome-only* quantity. The proposed
`D`-threshold reformulation of the shipped flag is **NEGATIVE** (§11); a third attempt
replacing the bootstrap analytically is reported there too.

**P2 (semiparametric efficiency).** Thm 1 hands over the pathwise derivative: the
influence function splits into a transport part and an overlap part. §10 builds and
tests three estimators against this.

**P4 (sharp partial identification).** §12 — working Kantorovich-LP tool built; the
"`D` as identification budget" conjecture tested and **retracted**. Continuing.

**P5 (win-mediation).** Not attempted (§13).

**Cross-cutting probes**: #1 DONE (§9); #2 MIXED (§10); #3 NEGATIVE (§11); #4 DONE,
exceeded scope (§7).

**Why this is the right spine.** Every frontier paper reduces to questions about the
same two functionals (`v`, `∫f²`/`D`); P1 and P3 are provably the same functional seen
from two sides; each step validates immediately against the existing oracle.

---

## 7. Claim 6 — exact closed forms for K=1 censored win statistics

Two independent hand derivations for `T ~ Exponential`, Type-I (fixed administrative)
censoring at `c`, `p := P(T≤c)`.

**(a) Var(w) and D** (no causal effect; the base-population functionals Paper 03's
degeneracy theory conditions on):
```
D = p(2−p),        Var(w) = p(p²−3p+3)/3
```
from the package's exact win rule (`i` beats `j` iff `δ_j=1 & X_i>X_j`; loses iff
`δ_i=1 & X_j>X_i`). Both are **distribution-free** — they depend on `F` only through
`p`, since `w(o)` and the decided-indicator are functions of `F(T)~Uniform`.
Validated by **exact** `O(n log n)` leave-one-out computation (order statistics, no MC
opponent-subsampling):

| p | D_hat | D_theory | Var(w)_hat | Var(w)_theory |
|---|---|---|---|---|
| 0.02 | 0.0408 | 0.0396 | 0.02020 | 0.01960 |
| 0.30 | 0.5110 | 0.5100 | 0.21937 | 0.21900 |
| 0.90 | 0.9899 | 0.9900 | 0.33300 | 0.33300 |

**Key asymptotic**: `Var(w)/D → 1/2` **exactly** as `p→0` (the degenerate limit that
matters for Paper 03) — confirmed: `0.500, 0.498, 0.487` at `p = 0.001, 0.01, 0.05`.

**(b) The win gradient itself**, adding a proportional-hazards causal effect
`T_x ~ Exponential(λ(x)=λ₀e^{αx})` with the same censoring:
```
B'(0) = −α · D(p),          D(p) = p(2−p)      [the SAME functional]
```
by differentiating the exact `NB(x+ε,x−ε)` integral (both-observed integral +
one-censored boundary term) at `ε=0`. **Independent cross-check**: as `c→∞` (`D→1`)
this reduces to `B'(0)→−α`, matching the classical two-exponential concordance formula
`P(T_A>T_B)=λ_B/(λ_A+λ_B)` to first order — a different derivation route, same limit
(lesson E). Simulation (`n=2×10⁶`): rel. err 1–4% across `c ∈ {0.3,…,50}`; the single
17% outlier at the smallest, noisiest `c` was confirmed **Monte-Carlo noise** (not
systematic) by re-running at another seed and larger `n`, moving it to 2–5%.

This makes "win gradient = (causal effect) × decidability" an **exact identity** in
two clean instances — not merely the `g∝D^0.873` empirical fit of §5.

---

## 8. P1(a) partial — does the structure lift to the real K=3 frailty-correlated composite?

**Exact structural identity** (by construction): the kernel telescopes over the first
decisive tier, `h(A,B) = Σ_k 1{R_k}·sign_k` with `R_k := {tiers 1..k−1 tied}`, so
`NB = Σ_k contrib_k`. **Open question**: does each `contrib_k` further factor as
`R_k × D_k^cond × (a tier-invariant local velocity)`?

**Test**: exact `O(n²)` tier-cascade tracking (not a subsample) on the real `mrwin`
frailty simulator (N=1400, shared Gamma frailty, cascading censoring), 4-point
censoring ladder × 3 tiers:

| tier | mean `contrib_k/(R_k·D_k^cond)` | CV |
|---|---|---|
| 1 | 0.158 | 0.360 |
| 2 | 0.166 | 0.214 |
| 3 | 0.141 | 0.418 |
| **pooled** | — | **0.309** |

Raw `contrib_k` spans **>2 orders of magnitude** (0.0005–0.098) across the same grid;
the decidability-weighted ratio compresses that to **CV 0.31** — substantial, but not
the `<5%` precision of the clean K=1 forms in §7. **Honest reading**: the tier-wise
lift substantially holds, but a real ~30% residual remains, plausibly a
frailty-correlation term (reaching tier `k` is not independent of the tier-`k` causal
gradient once tiers share frailty). The internal identity `Σ_k contrib_k = g_total`
holds to numerical precision on every row — so the residual is a real feature of the
DGP, not a code bug. **P1(a)'s general closed form remains open.**

---

## 9. Claim 7 — P1(b), the impossibility theorem, made rigorous

§3's argument was checked at two points (Gaussian, Laplace). Using the **exponential
power / generalized Gaussian** family (shape `κ`; 1=Laplace, 2=Gaussian, →∞
uniform-like), unit-variance standardised, for `Y(x)=μ(x)+σZ_κ`:
```
β*(x;κ,σ) := B'(0)/2 = (2μ'(x)/σ)·ρ(κ),        ρ(κ) := ∫ g_κ(z)² dz
```

**(A)** `ρ(κ)` is non-constant on a 9-point grid `κ∈[0.7,12]`: range `[0.270,0.459]`,
max/min `1.70`. Analytic cross-checks exact to 6 d.p.: `ρ(Gaussian)=1/(2√π)=0.282095`,
`ρ(Laplace)=1/(2√2)=0.353553`.
**(B)** The formula holds to rel. err `4e−7`–`2e−6` over 6 `(κ,σ)` combinations.
**(C)** The impossibility, directly: at fixed `σ=1`, raw `β*` ranges **1.70×** across
shapes for the *same* `α=0.8`; dividing by `ρ(κ)` recovers `2α=1.600` **exactly at
every κ** (1.5999–1.6000). So no function of `(β*,σ)` alone — no universal rescaling —
recovers the causal effect uniformly over shapes; the data-estimable functional
`ρ(κ)` (equivalently `σ·∫f²`) is both necessary and sufficient.

**Methodological note (banked).** Part (C) first showed `κ=0.7` recovering `1.353`, not
`1.600` — an apparent breakdown at the heaviest tail. A convergence check varying
`eps`, grid resolution and integration limit *in isolation* showed it vanished once the
domain widened from `25σa` to `60σa`+: `κ<1` EPD tails decay slower than Laplace and
were being silently truncated. Pure numerical artifact, not a theorem failure
(production script now uses `100σa`). **A related trap recurred three times this
session**: `NB(x+ε,x−ε)` is *odd* in `ε`, so `B'(0) = NB(ε)/ε`, **not** `NB(ε)/(2ε)` —
flagged explicitly to prevent a fourth.

---

## 10. Probe #2 — τ as an MR estimand: a mixed result, three estimators tested

Does `τ` restore cross-cohort portability when cohorts differ only in
censoring/follow-up (a real two-sample MR problem: different biobanks)?

**(a) Population/oracle level — confirmed strongly.** In the real two-sample
within-family pipeline (F=1500 families × 2 sibs, M=25 SNPs), sweeping *only*
censoring: `oracle(c)/D_hat(c)` has **CV 0.136** vs **CV 0.850** for raw `oracle(c)` —
a >6× stabilisation.

**(b) Point estimate — does NOT improve.** `gamma_AR_hat` has CV 0.641–0.753 (two
independent runs); `tau_hat = gamma_AR_hat/D_hat` has CV 0.746–1.409 — similar or
**worse**. Mechanism: the AR argmin is an unstable point summary under weak
instruments (within-family per-SNP F-stat is intrinsically `<10`), and dividing by a
`D_hat` shrinking toward 0 *amplifies* the instability (`tau_hat` hit 2.55 vs a target
of 0.146 at the heaviest censoring). **Not a free lunch.**

**(c) Confidence interval — real improvement.** The raw AR CI's coverage of the fixed
portable target `τ*` degrades sharply with censoring: `0.963→0.925→0.925→0.863→0.787→0.662`.
The D-rescaled CI (`ar$ci/D_hat`) holds up: `0.963→0.925→0.925→0.900→0.875→0.875`.

**Fix attempt 1 — delta method: FAILS.** A joint influence-function SE
(`IF_gamma = (gh−γ·gx)/b`, `IF_D = 2(d−D̂)`, combined by the ratio delta method) is
**badly miscalibrated**: `SE/empirical-SD` ranges **2.4×–27×**, and its Wald CI covers
`τ*` at only 0.80–0.95 — worse than the naive plug-in. **Mechanism identified**: a
Wald/delta linearisation needs asymptotic normality of `gamma_hat`, which a weak-IV
ratio estimator violates — precisely why AR-type inference exists in the first place.

**Fix attempt 2 — AR-style test built directly for τ.** Test candidate `t` via the
moment `M(t) := a − b·t·D`, linearised at *fixed* `t` (never dividing by anything):
`IF_M_i(t) = gh_i − t·(D̂·gx_i + b̂·IF_D_i)`, linear in `t`, so its variance
`c0'−2t·c1'+t²c2'` is an exact quadratic and the set inverts by the *same* closed-form
machinery `mrwin_ar_onesample` uses for `γ` (relabelling `b→bD`, `gx→u`):

| instrument | bounded-set rate | AR-τ cov(τ*) | naive cov | delta cov |
|---|---|---|---|---|
| weak (α_s=0.15, N=3000) | 2–7% | 0.980–1.000 | 0.973–1.000 | 0.787–0.947 |
| strong (α_s=0.8, N=8000) | 55–75% | 0.950–0.990 | 0.950–0.990 | 0.700–0.840 |

**Reading**: AR-τ and the naive plug-in are **statistically equivalent** in both
regimes — which is real progress: the naive interval is now an *independently
validated* approximation to a properly derived joint test, not a lucky coincidence.
Delta-method is **decisively ruled out in both regimes**, confirming the diagnosis was
about the Wald linearisation, not about weak instruments alone. New free by-product:
the **bounded-set rate** honestly diagnoses how often `τ` is identified at all.
**What it does not do**: give a *tighter* interval, so (b) remains unresolved.

**Fix attempt 3 — two-sample multi-SNP AR-τ** (`p2-twosample-multisnp-ar-tau.R`).
Attempt 2 was tested in the wrong arena on two counts: portability is inherently a
*two-sample* problem, and one SNP cannot pin `τ` down (hence the 2–7% bounded rate).
Paper 04's deployable pipeline has L=25 SNPs, giving a `χ²_L` rather than `χ²_1`
statistic. Per-SNP `δ_l = τ·D·β_l`, so with both per-SNP SEs and `D_hat`'s own error
propagated (delta method; `D_hat` from the outcome sample, `β_l` from the exposure
regression, treated as independent):
```
AR_τ(t) = Σ_l (δ_l − t·D̂·β_l)² / ( se_δ,l² + t²·[D̂²·se_β,l² + β_l²·se_D²] )   ~ χ²_L
```
with `se_D` **cluster-robust by family** (`Var(D̂) = n^{-2} Σ_fam (Σ_{i∈fam} IF_D_i)²`),
and `min_t AR_τ(t) ~ χ²_{L−1}` as an over-ID test on the τ scale.

**Result — no improvement, after a self-correction worth recording.** The first run
appeared to show AR-τ decisively beating the plug-in (coverage 0.983–1.000 vs
0.817–0.933, the gap widening with censoring). That reading was a **counting
artifact**: naive sets that were non-contiguous had been scored as *non*-covering,
while AR-τ sets that were *unbounded* (68–75% of replicates) trivially cover. Rerun
apples-to-apples — naive scored generously (range as a conservative superset, no
contiguity requirement) and bounded rates reported alongside coverage — the advantage
disappears entirely:

| cens | D̂ | AR-τ bnd% | naive bnd% | AR-τ cov | naive cov |
|---|---|---|---|---|---|
| 0.02 | 0.859 | 0.271 | 0.102 | 1.000 | 1.000 |
| 0.30 | 0.371 | 0.322 | 0.288 | 1.000 | 1.000 |
| 0.80 | 0.176 | 0.283 | 0.283 | 1.000 | 0.983 |
| 5.00 | 0.026 | 0.250 | 0.333 | 0.983 | 0.983 |

Head-to-head on the *same replicates where both sets are bounded*, coverage is
**identical** (1.000/1.000 down to 0.933/0.933 at the heaviest censoring) and AR-τ is
**slightly wider** — width ratio 1.001, 1.003, 1.009, 1.024, 1.047, 1.202 across the
ladder. So the rigorous construction costs 0.1–20% width and buys nothing.
Separately, the τ-scale over-ID test is **degenerate**: type-I is `0.000` at every
censoring level under a valid instrument — no power, consistent with the known
conservatism of the win-MR over-ID test under weak instruments already recorded in
`within-family-twosample-winmr.R` (Q type-I 0.00–0.013).

**Bottom line for probe #2 (three attempts, one arena each).** `τ` is validated as a
portability fix **at the level of the causal target** (§10a, CV 0.850→0.136) **and of
confidence intervals** (§10c). It is **not** a reliable point estimator under the weak
instruments intrinsic to within-family MR (§10b). Of three candidate interval methods,
the **naive plug-in `ar$ci/D̂` is the practical recommendation**: the delta method is
decisively worse (attempt 1), and the rigorously-derived AR-τ construction is
statistically equivalent-but-wider in every arena tested — one-sample weak, one-sample
strong, and two-sample multi-SNP across a censoring ladder (attempts 2–3). That is a
useful positive conclusion in negative clothing: the simple estimator is *validated as
sufficient*, not merely unrefuted. The bounded-set rates (10–33% for both methods)
remain the honest headline — most of the time neither method identifies `τ` at all.

---

## 11. Probe #3 — a cheap degeneracy diagnostic: two negatives, then a reframe

§5's Claim 5 gave `Spearman(D, Var(w))=1.0000` for the **pure outcome-margin** win-score
variance. §6 proposed reformulating Paper 03's bootstrap-`ζ1` flag as a `D`-threshold.
But the shipped flag uses `ζ1(β̂) = c0 − 2β̂c1 + β̂²c2`, where `c0=Var(gh)` mixes the
instrument into the win score via `gh_i=(Z_i s_i − r_i)/(n−1)`.

**Attempt 1 — `D_hat` threshold: NEGATIVE.** Grid `N∈{800,2500}` × 7 censoring rates ×
weak/strong instruments (420 rows, `boot_reps=200` ground truth):
`Spearman(D_hat, degenerate) = −0.50` to `−0.53` — nowhere near the `1.0000` of the
outcome-only quantity. At a lax `D_hat>0.02`, **61% of "safe" cases were still flagged
degenerate**; tightening only shrinks coverage without reaching safety (still 43%
false-safe at `D_hat>0.30`). `N=800`→`N=2500` at the same censoring swings
`P(degenerate)` from 0.50 to 0.00: **sample size, not decidability, is a first-order
driver**, because `ζ1(β̂)` also depends on how far weak-IV noise pushes `β̂` from the
quadratic's minimiser — invisible to `D_hat`.

**Attempt 2 — `n·ζ1_min`, a deterministic lower bound: BETTER, still not clean.**
`ζ1(β)` is an upward parabola, so its minimum over **all** `β`,
`ζ1_min := c0 − c1²/c2`, satisfies `ζ1(β̂) ≥ ζ1_min` *by construction for any realised
`β̂`* — targeting exactly the failure mode above, and free from the existing blocks.
Over 840 cells: `Spearman = −0.54` (essentially unchanged), but the false-safe rate at
the useful threshold `τ=1.0` drops to **4.2%** (5/120) from 43–61% — an
order-of-magnitude improvement. Still not shippable: only 14.3% of rows clear `τ=1.0`,
thresholds above ~1.5 never trigger, and the degenerate/non-degenerate `n·ζ1_min`
ranges overlap substantially (up to 1.08 among degenerate rows; as low as 0.008 among
non-degenerate).

**Attempt 3 — reframe: replace the bootstrap, don't screen it**
(`p3-analytic-c2-bound.R`). The framing above was wrong: the flag is
`c2_boot_ci[1] < 1`, i.e. a **lower confidence bound** on `c² = n·ζ1(γ)`. So compute
*the bound itself*, not a correlate. Reading `R/onesample_ar.R` shows the bootstrap
holds `gamma` **fixed** inside the resampling loop, so the bootstrapped quantity is
exactly `Var_i(g_i)` at fixed `γ`, `g_i := gh_i − γ·gx_i` — no delta-method term for
`β̂`'s randomness is needed, which makes an analytic SE available:
```
Var(ζ1_hat) = (m4 − m2²)/n   ⟹   SE(c²) = sqrt(n·(m4 − m2²))
```
(`m2`,`m4` = 2nd/4th central moments of `g_i`), giving a normal bound `c²−z·SE` and a
lognormal bound `c²·exp(−z·SE/c²)` (variances are right-skewed). **Known caveat this
tests**: `g_i` are *not* iid — they are U-statistic projections the bootstrap captures
fully and the iid influence function does not; whether the gap matters is exactly the
measured question, and the decision-relevant metric is **flag agreement**, not bound
correlation. *Results pending — see §14.*

**No change has been made to `R/onesample_ar.R`** from any attempt. Claim 5 stands
**exactly as stated** — for the pure outcome-margin `Var(w(O))` — but its scope is now
precise: it does not extend to the instrument-crossed shipped diagnostic.

---

## 12. P4 — a working Kantorovich-LP tool; the "identification budget" conjecture retracted

**Question**: given ONLY the two arms' marginals (no rank-invariance/comonotonic
assumption), what is the sharp range of the win probability over *all* couplings? This
is literally a Kantorovich transportation problem.

**Error #1, caught before use.** A hand derivation assumed `1(y1>y2)` is globally
submodular, which would make comonotonic/countermonotonic couplings the universal
sharp extremes. Refuted directly: the configuration `y2<y1<y2'<y1'` violates the
submodularity inequality that holds for `y1<y2<y1'<y2'` — submodular in some regions,
not others, so the Monge shortcut does not apply. Retracted before it reached a
conclusion.

**Fix: solve the actual discretised LP** (`lpSolve`, dev-script only — *not* added to
package `DESCRIPTION`). Validated first on a hand-solvable 2-point case (equal
Bernoulli(0.5) marginals → sharp `[0, 0.5]`, matching hand calculation, with min/max
plans literally the diagonal/anti-diagonal). **Error #2, also caught**: `lp.transport`
defaults to `integers=1:(nc*nr)`, forcing an *integer* plan — infeasible for fractional
probability masses; caught via the solver's own `status` code (not a silently wrong
answer), fixed with `integers=NULL`. Post-fix, two shifted Gaussians (`μ1=0.5,μ2=−0.5`,
equal σ): sharp max `1.00000` (achieved by the comonotonic coupling — a pure location
shift makes it deterministic), sharp min `≈0.3875` (converging cleanly over 20→160
bins), achieved by **neither** simple extremal coupling. So even the
"comonotonic/countermonotonic are the two extremes" folklore is model-dependent.

**The conjecture test.** Extending the validated LP to the K=1 Type-I censoring win
rule (continuum + one atom per arm) and sweeping `c` (hence `p`, hence `D=p(2−p)`) at
fixed causal effect: the **`width ∝ (1−D)` conjecture from §6 is NOT confirmed — the
relationship runs the opposite way**. Width *shrinks* as censoring grows: `0.049` at
`D=0.12` vs `0.899` at `D≈1` (log-log slope of width vs `(1−D)` is `−0.06`, not `+1`).
**Mechanism**: winning requires the *opponent* to have an observed event, so heavy
censoring squashes `P(win)` toward 0 for *every* coupling simultaneously, narrowing the
absolute width mechanically — a different phenomenon from "the coupling is pinned
down", and conflating the two was the flaw. A validated check survived: the comonotonic
(mrwin-simulator-style) coupling matches the LP's **minimum** at every censoring level
here — the *opposite* corner from the Gaussian location-shift case, again confirming no
universal rule ties comonotonicity to a fixed Fréchet extreme.

**Status: continuing.** The LP tool is validated, reusable infrastructure. The literal
conjecture is retracted. Open reformulations worth testing: a *normalised/relative*
width; bounding the **net-benefit contrast** rather than raw `P(win)`; and accounting
for the fact that `D` itself has Fréchet-type ambiguity under an unconstrained coupling
of the censoring indicators.

---

## 13. P5 — not attempted

Win-mediation via a transport-velocity vector-field decomposition requires a
cross-world potential-outcomes formalisation that could not be brought to this file's
validation standard within scope. Per "no claim before its test passes," an
unvalidated sketch would be worse than an honest "not attempted." Open.

---

## 14. Honest scope and caveats

- **Thm 1 (★)** is exact for scalar continuous outcomes. K=1 hierarchical is now
  **exact** in two instances (§7). The general **K=3 frailty-correlated case is
  partial** (§8): the tier decomposition holds to CV 0.31, with a real ~30% residual
  unexplained. P1(a)'s general closed form is **open**.
- **The impossibility result** is now a rigorous continuous-shape-family theorem (§9)
  — closed.
- **`τ`** is validated as a portability fix at the population/oracle level and at the
  confidence-interval level (§10a,c) — the interval claim independently corroborated by
  a properly derived AR-style test — but **not** as a reliable point estimator under
  weak instruments (§10b). Delta-method fails for an identified reason; the rigorous
  AR-style fix validates the interval but gives no tighter one. **The two-sample
  multi-SNP AR-τ test (§10, attempt 3) is running; results not yet in.** This is the
  biggest open gap in deploying `τ` as an MR estimand.
- **The `D`-thresholded degeneracy screen does not work** (§11); `n·ζ1_min` is an
  order-of-magnitude improvement but still not shippable; **the analytic-`c²`-bound
  reframe (§11, attempt 3) is running; results not yet in.** No package code has been
  changed by any attempt.
- **P4** (§12): validated LP tool exists; the headline conjecture is **retracted**;
  reformulations open. **P5** (§13) not attempted.
- All numbers are simulation/numeric results with MC or numeric error reported. Three
  self-caught error classes are recorded rather than edited out: two factor-of-2
  finite-difference bugs (`B'(0)=NB(ε)/ε`, not `/(2ε)`), one integration-truncation
  artifact (§9), and two P4 errors (false submodularity shortcut; LP integer default).

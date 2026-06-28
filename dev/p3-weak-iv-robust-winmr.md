# P3 — Identification-robust (weak-instrument) inference for win-ratio MR

Deep-dive spec (2026-06-23). Theory note for Paper 03. Extends Paper 01 (Fieller)
and Paper 02 (two-sample win-MR). Status: derivation + plan; **no claim is final
until its probe passes** (CLAUDE.md §1). All numbered results are targets to prove
or to falsify by simulation against the interventional oracle.

---

## 0. The problem and the one-line idea

Win-ratio MR is a **ratio** estimand: `γ = (instrument → win) / (instrument →
exposure)`. Like every IV ratio it is fragile when the instrument is weak
(denominator `→ 0`): the Wald/delta SE `∝ 1/|Cov(Z,X)|` explodes — exactly the
defect Paper 01 documented (`cor(log se, log min|ΔX|) = −0.91`) and patched with
Fieller. But Paper 01's Fieller is a *linear/Gaussian* construction on a
plug-in covariance; it is not the right object for a **pairwise U-statistic**
moment, and it has no theory for the over-identified (many-SNP) or
near-irrelevant-instrument regimes.

**Idea.** Do not estimate the ratio. Test `H0: γ = γ0` through a pairwise
**Anderson–Rubin (AR) moment** that equals 0 under the null *regardless of
instrument strength*, and invert the test. The moment is a degree-2 U-statistic;
the genuinely new mathematics is making the test **uniformly valid across the
strong → weak → degenerate identification continuum**, where the U-statistic can
lose its first-order (Hájek) signal and acquire a non-Gaussian limit — a regime
that has no analog in the linear-IV AR literature (whose moment is a sample mean,
never degenerate).

---

## 1. Setup, kernel, and the moment identity

I.i.d. `W_i = (O_i, X_i, Z_i)`: hierarchical outcome `O_i` (K priority levels,
times + status), exposure `X_i`, instrument `Z_i` (scalar PRS in §1–4; vector of
`L` SNP dosages in §5). The **win kernel** `h(O_i,O_j) ∈ {−1,0,+1}` is
antisymmetric (`h_ji = −h_ij`, `h(o,o)=0`). Define the **win score** (Hájek
projection of the kernel):

```
w(o) := E[ h(o, O_2) ] = P(o ≻ O_2) − P(o ≺ O_2),     E[w(O)] = 0.
```

**Scale choice — net benefit, not log-odds.** We target the **net-benefit
gradient** `β` (the per-unit-exposure shift in net-win-probability), not the
log-win-odds. This is deliberate: it makes the identifying moment *exactly linear*
in the kernel `h`, so the U-statistic theory is exact and the ratio-of-logs
nonlinearity that made Paper 01's delta method fragile never appears. The
log-win-odds `γ` of Paper 01 is recovered post hoc by the smooth monotone link
`γ = link(β)` (delta method for the point; the *test* lives on the β scale where it
is exact). Net benefit is Buyse's parameterization — position vs Buyse, generalized
pairwise comparisons.

**Identifying moment.** Under the MR assumptions (relevance, independence,
exclusion) plus the local structural model (M) of Paper 01 and instrument–confounder
independence,

```
ψ(β) := E[ {h(O_i,O_j) − β (X_i − X_j)} (Z_i − Z_j) ] = 0   at  β = β*.     (1)
```

Expanding, `ψ(β) = ψ_h − β ψ_x` with

```
ψ_h = E[h(O_1,O_2)(Z_1−Z_2)] = 2 Cov(Z, w(O)),   ψ_x = 2 Cov(Z, X),
```

so the point estimand is the **win-statistic Wald ratio**

```
β* = ψ_h / ψ_x = Cov(Z, w(O)) / Cov(Z, X).                                  (2)
```

(2) is the single-moment limit that unifies Paper 01's per-stratum ISG and
Paper 02's per-SNP IVW: instrument's covariance with the *win score* over its
covariance with the *exposure*. The win score `w(O_i)` is the per-subject
net-win-vs-everyone — precisely what `mrwin_subject_win_loss_cpp` (built this
session) computes.

---

## 2. The AR statistic and the CI by inversion

Symmetric U-statistic kernel `φ_β(W_i,W_j) := {h_ij − β(X_i−X_j)}(Z_i−Z_j)`
(symmetric: both `h` and `(Z_i−Z_j)` flip under swap, product invariant). Estimator

```
U_n(β) = binom(n,2)^{-1} Σ_{i<j} φ_β(W_i,W_j) = U_n^h − β U_n^x,             (3)
```

linear in β, with `U_n^h` the U-stat of `h_ij(Z_i−Z_j)` and `U_n^x` of
`(X_i−X_j)(Z_i−Z_j)` (note `U_n^x = 2·\widehat{Cov}(Z,X)`).

**Hoeffding CLT.** With first projection `g_β(w) := E[φ_β(w,W_2)]` and
`ζ_1(β) := Var(g_β(W_1))`, if `ζ_1(β0) > 0`,

```
√n ( U_n(β0) − ψ(β0) ) →d N(0, 4 ζ_1(β0)).
```

Under `H0: β = β0` (structural model true), `ψ(β0)=0` **exactly**, for any
instrument strength. Hence the **AR test**

```
AR_n(β0) := n U_n(β0)^2 / (4 ζ̂_1(β0))  →d  χ²_1   under H0,                (4)
```

and the **identification-robust CI**

```
CI_AR = { β0 : AR_n(β0) ≤ q_{χ²_1, 0.95} }.                                 (5)
```

Because `U_n(β0)` is linear and `ζ̂_1(β0)` quadratic in β0, (5) is the solution of a
quadratic inequality — a bounded interval, a two-ray complement, or the whole line:
**the same geometry as Fieller, but with the U-statistic-correct variance.** P3 is
therefore the rigorous nonparametric generalization of Paper 01's Fieller interval.

**Why it fixes the weak instrument (the point).** The Wald CI inverts
`(β̂ − β0)/se(β̂)` with `se(β̂) ∝ 1/|U_n^x|` → ∞ as `Cov(Z,X) → 0` (denominator
blow-up). The AR CI (5) never divides by `U_n^x`; `ζ̂_1(β0)` stays bounded away from
0 generically (its `Z_i w(O_i)` term does not vanish with `Cov(Z,X)`), so the CI
keeps **correct coverage**, returning a wide/unbounded set *honestly* under a weak
instrument. This is the rigorous justification for the empirical fact validated
this session — "unbounded covers; coverage ≈ 0.96 at every instrument strength" —
which Fieller produced only heuristically.

---

## 3. Variance estimation (reuses the shipped C++ primitive)

First projection in closed form (`μ_X,μ_Z` means, `ρ(o):=E_j[h(o,O_j)Z_j]`):

```
g_β(W_i) = Z_i w(O_i) − ρ(O_i) − β[ X_i Z_i − μ_Z X_i − μ_X Z_i + E(XZ) ].   (6)
```

Plug-in: `ĝ_β(W_i) = (n−1)^{-1} Σ_{j≠i} φ_β(W_i,W_j)`, then
`ζ̂_1(β) = (n−1)^{-1} Σ_i (ĝ_β(W_i) − \bar{ĝ})^2`. The per-subject aggregates
`w(O_i)` (net wins/losses vs all others) and `ρ(O_i)` (the same, weighted by `Z_j`)
are exactly the outputs of `mrwin_subject_win_loss_cpp` (unweighted → `w`;
weights `= Z_j` → `ρ`). So the whole test is **implementable now**, O(n·m) with the
subject-subsample trick from the SE work, or exact O(n²) at moderate n.

---

## 4. Connection to Paper 01's Fieller (a strict generalization)

Fieller tests the linear contrast `ψ_h − β ψ_x = 0` using a *plug-in Gaussian*
covariance of `(ψ̂_h, ψ̂_x)`. (4) tests the same contrast with the *U-statistic*
covariance `4ζ_1(β0)`, which correctly accounts for the pairwise dependence the
plug-in ignores. In the strong-instrument, large-`n` limit the two coincide; in the
weak / small-`n` / heavy-tie regimes the AR variance is the right one. So Paper 03
*contains* Paper 01's CI as its Gaussian-approximation special case.

---

## 5. Over-identification: vector instrument and a weak-IV-robust pleiotropy test

`Z_i ∈ R^L` (e.g. `L` SNP dosages). Vector moment `ψ(β) ∈ R^L`,
`U_n(β) ∈ R^L`, asymptotic covariance `Σ(β) = 4 Cov(g_β(W_1))` (L×L). The
**over-identified AR**

```
AR_n(β0) = n · U_n(β0)' Σ̂(β0)^{-1} U_n(β0)  →d  χ²_L   under H0,            (7)
CI = { β0 : AR_n(β0) ≤ q_{χ²_L,0.95} }.
```

Two outputs for the price of one:

1. **Weak-IV-robust CI** for β via (7) — valid with `L` possibly-weak instruments,
   no first-stage F gate.
2. **Win-statistic over-identification / pleiotropy test.** `min_β AR_n(β) →d
   χ²_{L−1}` under instrument validity — a Sargan/Hansen-J / Cochran-Q analogue.
   Rejection signals a violated exclusion restriction (pleiotropy). This is the
   *identification-robust* counterpart of the SDPD diagnostic (§7.2 of Paper 01),
   and — unlike SDPD — it inherits the InSIDE geometry automatically: pleiotropy
   proportional to instrument strength shifts every moment by the same `β`, so it is
   absorbed into `β̂` and *not* flagged (correct), whereas InSIDE-satisfying
   pleiotropy inflates `min_β AR_n` (flagged). P3 thus unifies the weak-instrument
   and pleiotropy threads in one quadratic form.

---

## 6. The novel core: uniform validity across the identification continuum

Everything above is "standard non-degenerate U-statistic AR." The genuinely new,
hard contribution is the **boundary**.

`ζ_1(β0)` can vanish: when the first projection `g_β0(·)` is a.s. constant, the
Hájek term is gone and `U_n(β0)` is a **degenerate U-statistic** with a
*non-Gaussian* limit

```
n U_n(β0) →d Σ_k λ_k (χ²_{1,k} − 1),                                        (8)
```

the `λ_k` being eigenvalues of the second-projection operator
`φ̃_β0(w,w') = φ_β0 − g_β0(w) − g_β0(w') + ψ(β0)`. In this regime the Studentized
`AR_n = nU_n^2/(4ζ̂_1)` has a heavy-tailed, non-pivotal limit → **over-rejection**.
This happens precisely in the most MR-relevant cases: near-irrelevant instruments
and balanced win comparisons (low decided-pair fraction under heavy censoring).

The linear-IV AR literature (Anderson–Rubin 1949; Staiger–Stock; Moreira's CLR)
**never meets this**: their moment is a sample mean, asymptotically Gaussian with a
fixed variance, never degenerate. So the required theory is new.

**Target theorem (the paper's heart).** A statistic `T_n(β0)` and critical value
that are **uniformly size-correct over a neighborhood of the degeneracy**. Sketch of
the construction:

- Estimate both variance components: the Hájek `ζ̂_1(β0)` (O(n) projections) and the
  degenerate `ζ̂_2(β0)` (the Hilbert–Schmidt norm / spectrum of the empirical
  second-projection operator, via incomplete-U or random-projection sketching to stay
  subquadratic).
- Use a **self-normalized interpolating statistic** whose limit law is `χ²_1` when
  `ζ_1` dominates and the weighted-χ² (8) when `ζ_1 → 0`, with a continuous
  transition, e.g. a studentization that adds the `(2/n)·`trace term:
  `T_n = nU_n(β0)^2 / (4ζ̂_1(β0) + (2/n) ζ̂_2(β0))`, calibrated by a
  multiplier/permutation bootstrap on the *projected* kernel to get a single critical
  value valid across the regime.
- Prove `sup_{P ∈ neighborhood} | P(T_n > c_α) − α | → 0` (uniform/honest size),
  by a drifting-sequence argument `ζ_1 = ζ_1(n) → 0` at every rate, the standard
  device for uniform inference.

If achieved, this is the first **degeneracy-robust Anderson–Rubin test** — a
contribution to U-statistic inference beyond MR. Falsifiable cleanly (see §8).

---

## 7. Two-sample / summary-data version (immediate, extends Paper 02)

In two-sample MR `ψ_h` and `ψ_x` come from independent cohorts. The per-SNP
win-odds GWAS coefficients `δ̂_ℓ` (with the **closed-form SE built this session**)
estimate the instrument→win moments; `β̂_GX,ℓ` the instrument→exposure moments.
Under `H0: β = β0`, the vector `r̂(β0) = δ̂ − β0 β̂_GX` has mean 0 and covariance
`V̂(β0) = V̂_δ + β0² V̂_{GX}` (two independent cohorts; both covariances already
produced by `mrwin_win_gwas` + the exposure GWAS). Then

```
AR_2s(β0) = r̂(β0)' V̂(β0)^{-1} r̂(β0)  →d  χ²_L,                            (9)
```

invert for the weak-IV-robust CI, and `min_β AR_2s` for the pleiotropy test. **(9)
is implementable today on the shipped `mrwin_win_gwas` + SE machinery** — the
fastest route to a first Paper-03 result, and a strict upgrade of Paper 02's IVW
(which uses point inverse-variance weights and is not weak-IV-robust).

---

## 8. Decisive probe first (CLAUDE.md lesson A)

Before any prose, run the smallest experiment that could **falsify** the headline.
Reuse the validated `alpha_s` instrument-strength sweep and oracle `γ*/β*`:

1. **Coverage vs Fieller across instrument strength.** Compute, per cohort, the
   Fieller CI (Paper 01) and the U-statistic AR CI (5). Falsify if AR does **not**
   hold ≈95% coverage uniformly (incl. weak `alpha_s = 0.1`) **and** is no wider
   than Fieller in the strong regime. Prediction: AR maintains nominal coverage where
   Fieller is heuristic; AR ⊆ Fieller asymptotically (efficiency).
2. **Degeneracy stress.** Drive toward the boundary (near-irrelevant instrument +
   heavy censoring → low decided-pair fraction). Falsify the *need* for §6 if naive
   `χ²_1`-AR keeps nominal size; confirm it if naive-AR over-rejects while the
   degeneracy-robust `T_n` does not.
3. **Over-ID pleiotropy.** Reuse the InSIDE / score-proportional pleiotropy DGPs
   (Paper 01 §7.2): `min_β AR_n` should be calibrated under InSIDE-violating
   (instrument-proportional) pleiotropy and powered under InSIDE-satisfying — a
   weak-IV-robust echo of the SDPD result.

Each probe uses only the existing kernel, oracle, and per-subject C++ primitive.

---

## 9. Open technical questions (honest gaps)

- Exact form and *consistent subquadratic estimator* of the degenerate variance
  `ζ_2` and its spectrum (§6) — incomplete-U vs random sketching trade-off.
- Behavior under **ties/censoring**: `h` mixes decided/tie pairs; does the
  decided-pair fraction `→ 0` (heavy censoring) drive degeneracy, and is the test
  still honest? (Interacts with the censoring-invariance problem, Paper 04 candidate.)
- **CLR (Moreira) analogue:** is there a conditional-likelihood-ratio refinement of
  (7) for the U-statistic moment that is more powerful than AR while staying robust?
  Likely yes; needs the conditional distribution of the score given a sufficient
  statistic for the (nuisance) instrument-strength — non-trivial for U-statistics.
- **Many weak instruments** (`L` large, each weak): does (7) need an `L`-correction
  (à la Kleibergen / many-instrument asymptotics) when `L/n` is non-negligible?
- Scale link `γ = link(β)`: confirm the monotone map net-benefit → log-win-odds and
  that the AR set transports correctly (it does, by monotonicity, but the endpoints
  need care at the unbounded branch).

## 10. Why this is the right first frontier paper

Tied to **de-risked, shipped** machinery (Paper 02 estimator + closed-form SE +
per-subject C++ + oracle); the near-term deliverable (§7) is implementable now; it
**strictly generalizes** Paper 01 (Fieller) and Paper 02 (IVW); and it carries one
**genuinely novel theorem** (§6, degeneracy-robust AR) that contributes to
U-statistic inference beyond MR. Bounded scope, decisive probes, high payoff.

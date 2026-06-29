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

**The unified limit (derived).** Index identification strength by a drifting
sequence `ζ_1 = ζ_1(n)` with `n ζ_1(n) → c² ∈ [0, ∞]`. Hoeffding-decompose
`φ_β0 − θ = g(W_i) + g(W_j) + ψ̃(W_i,W_j)` (canonical, `E[ψ̃|W_i]=0`), with degenerate
second moment `δ_2 = E[ψ̃²] = Σ_k λ_k²` (`λ_k` the eigenvalues of the operator
`Af(w)=∫ψ̃(w,w')f dP`). Then **on the `nU_n` scale**, under H0,

```
n U_n(β0)  →d  R(c) := N(0, 4c²) + Σ_k λ_k (Z_k² − 1),   Z_k iid N(0,1),       (8)
```

the two terms independent — the Gaussian from the degree-1 (Hájek) part, the
weighted-χ² from the degenerate part. Endpoints: `c=∞` ⇒ `R/2c → N(0,1)`, recovering
the χ²_1 AR; `c=0` ⇒ `R = Σλ_k(Z_k²−1)`, a pure degenerate-U limit. The exact
variance is `Var(U_n) = 4(n−2)ζ_1/(n(n−1)) + 2δ_2/(n(n−1))`, i.e.
`n²Var(U_n) → 4c² + 2δ_2 = Var(R(c))`.

**The right pivot.** `nU_n` is `O_p(√n)` when `ζ_1` is fixed but `O_p(1)` when
`ζ_1 = O(1/n)`; only the *studentized* `S_n² = U_n²/Var(U_n)` is `O_p(1)` uniformly,
with `S_n² →d R(c)²/Var(R(c))` — interpolating χ²_1 (c=∞) and the degenerate law (c=0).
A single χ²_1 critical value is therefore **wrong** in the degenerate regime.

**What the probe established** (`tools/validation-scripts/p3-degeneracy-robust-probe.R`,
synthetic kernel `φ=ε(A_i+A_j)+A_iA_j`, `ζ_1=ε²`, rank-1 degenerate part; H0 size,
target 0.05, n=300):

| ζ_1 | naive χ²_1 AR | oracle (true c, λ) | plug-in (est. c, λ) |
|---|---|---|---|
| 1.00   | 0.059 | **0.050** | 0.059 |
| 0.01   | 0.183 | **0.057** | 0.242 |
| 0.0025 | 0.196 | **0.061** | 0.001 |
| 0.00   | 0.181 | **0.059** | 0.000 |

Three conclusions, each evidence-backed:

1. **The limit (8) is correct.** The *oracle* test — reject iff `(nU_n)² >`
   0.95-quantile of `R(c_true, λ_true)²` — holds ~0.05 across the **entire**
   strong→degenerate continuum. The theory is right.
2. **Naive AR really does fail under degeneracy** — size inflates to ~0.18–0.20
   (≈4×) in the pure-degenerate synthetic. (In the realistic win-MR DGP the inflation
   is mild, §8.2 coverage 0.94, because the win kernel is not purely degenerate.)
3. **The hard obstruction is precisely the boundary nuisance `c² = nζ_1`.** It is
   **not consistently estimable** there: `ζ_1` and its plug-in bias `δ_2/(n−1)` are
   both `O(1/n)`, so `ζ̂_1^{bc} = \widehat{Var}(ĝ) − δ̂_2/(n−1)` is a difference of two
   same-order noisy quantities. The plug-in test is consequently **unstable**
   (0.242 → 0.000 across a tiny ζ_1 range), strictly worse than naive. The spectrum
   `λ̂_k` *is* estimable (its eigenvalues are `O(1)`); only `c²` is not.

**Resolution (the remaining math, now correctly scoped).** Because `c²` is not
estimable under the null, plug-in calibration cannot work — the fix must be
**robust-to-non-estimable-nuisance inference**, two viable routes:

- *Andrews–Cheng least-favorable / sup-over-CI.* Build a `(1−α_1)` confidence
  interval `[c²_lo, c²_hi]` for `c²` (it includes 0 at the boundary); take the
  critical value `= max_{c²∈[c²_lo,c²_hi]}` of the `(1−α_2)` quantile of
  `R(c,λ̂)²/Var(R(c,λ̂))`, with `α_1+α_2=α` (Bonferroni). Uniformly valid by the
  drifting-sequence argument; needs only the estimable spectrum `λ̂` plus a CI for
  `c²` (not a point estimate).
- *Degenerate-U bootstrap* (Arcones–Giné): the ordinary nonparametric bootstrap is
  inconsistent for degenerate U-statistics; the corrected (canonical-kernel)
  multiplier bootstrap reproduces (8) and yields the critical value directly.

Either gives the first **degeneracy-robust Anderson–Rubin test** — a contribution to
U-statistic inference beyond MR. The oracle result certifies the target is reachable;
the open work is purely the non-estimable-nuisance handling. Subquadratic spectrum:
top-`r` `λ̂_k` via Nyström/randomized SVD on the centered kernel, tail as a trace
correction.

**Sup-over-CI: implemented and probed (2026-06-23, `tools/.../p3-degeneracy-supci-probe.R`).**
The least-favorable test was built and run. Key simplification: `q(c²) =`
0.95-quantile of `R(c)²` is **monotone increasing** in `c²` (the Gaussian variance
`4c²` grows), so the sup over the CI is attained at its **upper endpoint** `c²_hi`
— no grid needed. The CI for `c²` must be built by **bootstrap**, not the analytic
iid se: the empirical projections `ĝ_i` are *dependent* through the degenerate
common mode, so the sample-variance se underestimates `sd(ĉ²)` by ~10× near the
boundary (an analytic-se version over-rejected, 0.17 at `ζ_1=0.01`). Result (H0
size, target 0.05; naive / oracle / sup-CI):

| ζ_1 | naive | oracle | sup-CI |
|---|---|---|---|
| 1.00   | 0.062 | 0.049 | 0.016 |
| 0.09   | 0.086 | 0.052 | 0.008 |
| 0.01   | 0.186 | 0.054 | 0.000 |
| ≤0.0025| ~0.19 | ~0.06 | 0.000 |

**Verdict: the sup-CI test is uniformly VALID** (size ≤ 0.05 everywhere — it *fixes*
naive AR's 0.18–0.20 over-rejection under degeneracy) **but CONSERVATIVE, and this
is fundamental.** Near the boundary `c² = nζ_1` is *weakly identified* — its
estimation error is `O(1)` relative to its value (`sd(ĉ²) ≈ 3.5` when `c² ≈ 3`) — so
the valid upper limit `c²_hi` is necessarily large, inflating the critical value.
Tuning (analytic vs bootstrap se, percentile vs `+z·se`, the `α_1/α_2` split) trades
the over-rejection for conservativeness but cannot recover the oracle's exact size:
the oracle *knows* `c²`; no data-driven test can match it where `c²` is
unidentified. This is a genuine **validity–power frontier at the degenerate
boundary**, itself a reportable result (you pay, in CI width, exactly when you
cannot tell how strong the instrument is — the honest behaviour).

**Practical interim recommendation.** Since §8.2 shows the realistic-regime
degradation is mild (size ≈0.06), ship the basic AR test (§2) with a
**`weak_degenerate` flag** raised when `nζ̂_1` is small relative to the spectral
mass `Σλ̂_k²`; reserve the full sup-over-CI statistic for the flagged extreme tail.

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

### Probe results — EXECUTED 2026-06-23 (`tools/validation-scripts/p3-weak-iv-robust-probe.R`)

`beta*` (net-benefit gradient): interventional oracle **0.111** vs moment-ratio
plim **0.117** — agree, confirming the moment `beta* = Cov(Z,w(O))/Cov(Z,X)` IS the
causal net-benefit gradient (eq. 2).

**8.1 — AR vs Wald coverage (N=4000, R=300, Z=PRS).**

| alpha_s | AR cover | AR unbounded | AR med-width | Wald cover | Wald med-width |
|---|---|---|---|---|---|
| 0.10 | **0.947** | 0.73 | 1.27 | 1.000 | 1.48 |
| 0.20 | **0.960** | 0.26 | 0.73 | 0.993 | 0.63 |
| 0.30 | **0.943** | 0.01 | 0.49 | 0.967 | 0.41 |
| 0.40 | **0.943** | 0.00 | 0.33 | 0.977 | 0.30 |

AR is calibrated ~0.95 at **every** instrument strength (weak included, via the
unbounded branch); Wald **over-covers throughout** (1.00 → 0.97) — the heavy-tailed,
denominator-driven mis-calibration Paper 01 documented. Headline confirmed: the AR
moment test removes the ratio fragility.

**8.2 — degeneracy stress (naive chi2_1 AR coverage, R=400).**

| alpha_s | cens | coverage |
|---|---|---|
| 0.05 | 0.05 | 0.955 |
| 0.02 | 0.05 | 0.958 |
| 0.05 | 0.40 | 0.960 |
| 0.02 | 0.40 | 0.958 |
| 0.05 | 1.00 | **0.940** |

Degeneracy is **real but MILD**: even under extreme weak-instrument + very-heavy
censoring (low decided-pair fraction) naive-AR only sags to 0.940. So §6's
degeneracy-robust statistic is a **refinement for the tail, not a blocker** — the
basic AR test is robust in practice. (This *simplifies* the paper: §6 becomes a
completeness/extreme-regime section, not the critical path.)

**8.3 — over-ID min-beta AR pleiotropy test (L=40 SNPs, N=2500, R=120, chi2_39).**

| scenario | reject |
|---|---|
| null (no pleiotropy) | **0.033** (calibrated; median stat 39.4 ≈ df 39) |
| InSIDE-violating (gamma_direct=0.3, ∝ instrument) | **0.042** (correctly ~null) |
| InSIDE-satisfying (tau=0.06, 30% of SNPs) | **0.833** (power) |

The identification-robust over-ID test is calibrated, correctly **blind** to
instrument-proportional (InSIDE-violating) pleiotropy — absorbed into `beta-hat` —
and **powerful** against InSIDE-satisfying pleiotropy: the AR echo of Paper 01's
SDPD result, now weak-instrument-robust. (A normalization bug — an extra factor
`n`, giving type-I = 1.000 — was caught by the probe's own type-I check and fixed:
`Sigma = 4*Cov(g)`, not `/n`. Lesson H: a numerical path is not trusted until its
calibration probe passes.)

**Verdict.** All three headline claims of P3 are confirmed empirically. AR fixes
the weak-instrument ratio fragility (uniform calibration where Wald is
mis-calibrated); the U-statistic degeneracy is mild; the over-ID statistic delivers
a calibrated, InSIDE-aware, weak-IV-robust pleiotropy test. The theory is
de-risked; the genuine open math (the *uniform* degeneracy-robust statistic §6, the
many-instrument refinement §9) is now scoped as refinement, not blocker.

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

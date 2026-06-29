# Research frontier — the hard papers (brainstorm, do not lose)

Status: vision note (2026-06-23). Ambitious open-problem map for the `mrwin`
programme. Grounded by a retrieval-based literature scan (scite + PubMed); see the
honest caveat at the end. Companion deep-dive: `dev/p3-weak-iv-robust-winmr.md`.

## Grand frame

The win-statistics literature is mature on **estimands, regression adjustment,
censoring corrections, power** (Mao & Wang, proportional win-fractions regression,
*Biometrics* 2021, doi:10.1111/biom.13382; Dong et al., IPCW win ratio, *J Biopharm
Stat* 2020, doi:10.1080/10543406.2020.1757692). But the **causal-IV / MR layer is
essentially virgin territory** — a direct "win ratio Mendelian randomization"
search returns only surname collisions (author "Win"), zero methodology.

> **Unifying thesis for the whole programme:** *causal inference for pairwise
> (degree-2 U-statistic) functionals of potential-outcome distributions.* Means are
> degree-1; win statistics are degree-2. Identification, semiparametric efficiency,
> partial identification, mediation, and weak-instrument-robust inference have all
> been built for degree-1 and **never lifted to degree-2**. That lift is the
> programme.

**Unique asset:** the project already has a *validated interventional oracle* +
fast kernel + influence-function variance, so every new estimand below is
*testable immediately* against ground truth — infrastructure no competing group has
for win-MR.

## ⭐ Newly probed breakthrough (2026-06-23): within-family win-ratio MR

**The structural identity.** A win statistic is a PAIRWISE comparison `h(O_i,O_j)`;
the strongest MR design — within-family / sibling — is a PAIRWISE genetic contrast
`Z_i − Z_j` between sibs. **They are the same object.** A sibling pair *is* a pair in
the degree-2 U-statistic. So restricting the Paper-3 pairwise IV moment to sibling
pairs makes the contrast obey Mendelian segregation, and `dZ` becomes **independent
of population stratification, assortative mating, and dynastic (indirect parental)
effects** — the confounders Papers 1–3 cannot control. This is the robustness apex of
the chain, and it falls out *for free* from the pairwise structure.

**Estimator.** `β*_wf = Σ_{(i,j)∈sib} h(O_i,O_j)(Z_i−Z_j) / Σ_{(i,j)∈sib}(X_i−X_j)(Z_i−Z_j)`.
The win kernel handles the hierarchical/censored outcome; the sib contrast handles
confounding.

**Probe result — survived falsification** (`tools/validation-scripts/within-family-winmr-probe.R`;
population stratification: stratum raises trait-raising alleles → systematic PRS
shift, and worsens the outcome → a pure `Z↔stratum→Y` backdoor; stratum not in X):

| estimator | gradient | bias vs oracle |
|---|---|---|
| interventional oracle | 0.099 | — |
| **population** win-MR (all individuals) | **−0.794** | **−0.893** (sign-flipped by confounding) |
| **within-family** win-MR (sib pairs) | **0.098** | **−0.001** (recovers the causal effect) |

Population win-MR is catastrophically biased; the within-family restriction recovers
the truth.

**New algorithm / theory needed (Paper 04 candidate).**
1. A **design-restricted (incomplete) U-statistic**: the pairwise sum runs over sib
   pairs only — `O(N)` pairs, not `O(N²)`.
2. **Family-clustered influence-function variance** and a **family-clustered
   Anderson–Rubin** test (Paper 3 carries over with the clustering): pairs within a
   family are dependent, families independent.
3. Extension to **>2 sibs** (all within-family pairs, weighted), parent–offspring
   trios, and the **two-sample within-sibship summary-data** form (within-sibship
   win-odds GWAS → AR pooling).
4. Formal robustness theorem: within-family win-MR identifies the causal win-odds
   gradient under stratification + assortative mating + dynastic effects.

Position vs the scalar within-sibship MR literature (Brumpton/Davies/Howe et al.):
the novelty is the **win / pairwise-functional** outcome and the recognition that the
statistic and the design are the same pairwise object. Deep-scan that literature
before claiming first (see caveat).

## Paper numbering (reconciled)

- **01** `01-methods-scalable-cwr` — scalable + calibrated DS-CWR (one-sample). *Draft + full validation done.*
- **02** two-sample / summary-data win-ratio MR — `mrwin_win_snp/_win_gwas/_twosample_ivw` + closed-form SE. *Estimator shipped + validated this session.* (Was loosely tagged "Paper 03"; it is Paper 02.)
- **03** identification-robust (weak-instrument) inference for win-ratio MR — **P3 below**; deep-dive in `dev/p3-weak-iv-robust-winmr.md`.
- **04+** the frontier papers P1, P2, P4, P5 below.

## The five hardest branches (each needs NEW math)

| id | title | open? (lit) | new math | depends |
|----|-------|-------------|----------|---------|
| **P1** | Collapsible causal win effect via optimal transport on prioritized outcomes | PARTIAL→ causal slice OPEN | lexicographic OT on a latent partial order; collapsible "transport win effect" (maps compose); IV-identified under rank-preservation | keystone |
| **P2** | Semiparametric efficiency & double robustness for instrumented win ratios | OPEN | EIF of a *ratio of pairwise U-statistics under IV*; degree-2 higher-order IF + IV projection; one-step/TMLE doubly-robust in (win-regression, instrument-propensity) | — |
| **P3** | Identification-robust (weak-instrument) inference for win-ratio MR | OPEN | Anderson–Rubin on a **U-statistic moment**; uniform validity across the strong→weak→degenerate identification continuum (degenerate-U-statistic-robust AR) | extends P02; ⭐ start |
| **P4** | Sharp partial identification of win effects under IV without functional form | OPEN | identified set = sup/inf of a bilinear win functional over couplings consistent with IV-stratified marginals + monotonicity → infinite-dim OT/LP with Kantorovich dual | shares OT with P1 |
| **P5** | Win-mediation: direct/indirect causal win effects | OPEN | cross-world *pairwise* potential outcomes; copula coupling of cross-world mediators; exact 3-way split direct + indirect + **non-collapsibility interaction** (→0 iff collapsible, linking to P1) | needs P1 |

### Why each is hard and what the theorem is

- **P1 (keystone).** Win odds is non-collapsible: two studies with identical
  per-subject effects but different baseline spreads report different win ratios →
  non-transportable, non-poolable. Mao–Wang give *conditional* (still
  non-collapsible) regression; no *collapsible causal* target exists. Target =
  optimal-transport displacement between `do(x)` and `do(x')` latent-utility laws on
  the partial-order quotient (ties = semicompeting). Theorem: (i) collapsible
  (`T_{x→x''}=T_{x'→x''}∘T_{x→x'}`); (ii) → win-odds gradient under location shift;
  (iii) IV-identified. Position vs *probabilistic index models* (Thas) and
  *generalized pairwise comparisons* (Buyse).
- **P2.** No EIF / efficiency bound exists for the IV-identified win ratio; DS-CWR's
  GLS pooling is ad hoc. Estimand solves the pairwise moment
  `E[(h_ij − γ(X_i−X_j))(Z_i−Z_j)]=0`; derive γ's EIF via degree-2 higher-order
  influence functions (Robins–van der Vaart) + IV projection; build √n-efficient
  doubly-robust one-step/TMLE.
- **P3.** See `dev/p3-weak-iv-robust-winmr.md`. The new core: AR/CLR theory is built
  for sample-mean (Gaussian, never-degenerate) moments; the win moment is a degree-2
  U-statistic that can **degenerate** at a boundary. Uniformly-valid inference across
  the degeneracy is the novel theorem.
- **P4.** Drop Paper 01's local-linear (M): the win effect is only partially
  identified. Sharp set = sup/inf of the bilinear win functional over joint
  potential-outcome laws consistent with IV-stratified marginals, exclusion, and
  monotonicity → an infinite-dimensional OT/LP with a closed-form-ish Kantorovich
  dual; collapses to a point under (M). New class: partial-ID for pairwise
  non-collapsible functionals.
- **P5.** Mediation g-formula/product does not compose for win ratios
  (non-collapsibility). Cross-world *pairwise* outcomes `Y_i(x,M_i(x'))` need a copula
  coupling the two pair members' cross-world mediators; exact 3-way decomposition
  with a non-collapsibility interaction term that vanishes iff the functional is
  collapsible (unifies with P1).

## Sequencing

```
   P1 (collapsible/OT, keystone) ──► P4 (bounds, shared OT)
        └───────────────────────────► P5 (mediation = P1 + interaction)
   P3 (weak-IV-robust) ──◄ extends Paper 02 (shipped)   ⭐ start now
   P2 (efficiency/DR)  ──► upgrades the estimators of P01 & P02
```

Start **P3** (decisive, tied to de-risked shipped work, bounded math target) in
parallel with seeding **P1** (longest lead time). Then P2 → P4 → P5.

## Honest novelty caveat (do not skip before committing a paper)

The "OPEN" verdicts are **retrieval-based**: no relevant paper surfaced across
varied queries — strong but not dispositive. Before committing each paper, run a
**deep scan under adjacent terminology**: *probabilistic index models* (Thas, De
Neve), *Mann–Whitney / concordance causal effects*, *generalized pairwise
comparisons* (Buyse), *prioritized outcomes*, *net benefit*. P1 especially borders
the probabilistic-index literature — position, do not claim blank. The win+MR gap
(Paper 02/03 area) is the most robust finding (surname-collision-checked). Never
claim "first" without the deep scan; cite real DOIs only.

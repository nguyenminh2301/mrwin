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
the within-family / sibling MR design is a PAIRWISE genetic contrast `Z_i − Z_j`
between sibs. **They are the same degree-2 U-statistic object.** Restricting the
Paper-3 pairwise IV moment to sibling pairs makes `dZ` obey Mendelian segregation.
Estimator `β*_wf = Σ_{sib} h(O_i,O_j)(Z_i−Z_j) / Σ_{sib}(X_i−X_j)(Z_i−Z_j)`.

**RIGOROUS re-check (2026-06-23, corrected scope; `tools/.../within-family-winmr-probe.R`).**
The first probe compared within-family only to a *naive, PC-unadjusted* population
estimator — a **partial strawman**, since real MR adjusts for ancestry principal
components. Giving population MR its standard PC defence and adding a no-confounding
control and a dynastic scenario (bias vs oracle):

| scenario | naive-pop | PC-adj-pop | within-fam |
|---|---|---|---|
| none (control) | +0.01 | +0.01 | −0.01 |
| population stratification | **−0.89** | **−0.14** | **+0.02** |
| dynastic / indirect genetic | −1.03 | **−1.06** | **−0.01** |

**Honest conclusion (scope corrected).** PC adjustment removes *most* of population
stratification (−0.89 → −0.14), so within-family's advantage there is modest. The
**genuine, unique value is DYNASTIC / indirect-genetic confounding — which PC
adjustment cannot remove (−1.06) and within-family recovers exactly (−0.01)** — and,
by the same Mendelian logic, ASSORTATIVE MATING. The breakthrough survives, but its
headline is "robust to the confounders PCs *cannot* fix (dynastic, AM)", **not**
"robust to stratification" (PCs largely handle that). The no-confounding control
confirms all estimators are correct.

**Now also tested (`tools/.../within-family-ar-inference.R`):**
- **Assortative mating** — within-family recovers the causal gradient under AM (a
  continuous family-level latent correlating both parents + a dynastic path, which PCs
  cannot capture): point bias +0.005. **Robust to AM ✓.**
- **Family-clustered AR inference** — with 2 sibs/family each family contributes ONE
  independent moment `m_f(β)=(h_f−β·dX_f)·dZ_f`, so the within-family AR is a clean
  i.i.d.-across-families test `AR(β)=(Σ_f m_f)²/(Σ_f m_f²)→χ²_1`. 95% CI coverage of
  the causal gradient ≈ **0.925–0.950** across none/dyn/am — **calibrated ✓**.
- **Honest caveat — within-sibship weak-instrument finite-sample bias.** A small point
  bias (~+0.02) is present even with NO confounding (the within-family genetic variance
  is small → a weak instrument), and it **shrinks with N** (F=3000→9000: bias
  0.022→0.007, coverage 0.925→0.940), confirming consistency. This is the known power
  cost of within-family designs: lower precision, more finite-sample bias, larger N
  needed.

**>2 sibs — DONE (`tools/.../within-family-multisib-clustered-ar.R`).** Aggregating the
dependent C(s,2) within-family pairs to the family level `m_f` gives an
i.i.d.-across-families moment, so `AR=(Σ_f m_f)²/(Σ_f m_f²)→χ²_1`. Proven **necessary
AND sufficient**: the naive pair-independent AR under-covers as s grows (s=4: 0.85)
while the family-clustered AR holds ~0.95 (s=2/3/4: 0.950/0.925/0.950); point bias
stable (+0.013, the within-sibship weak-IV bias, not growing with s).

**Efficiency — DONE, the question dissolves (two probes).**
- *Within a family* (`tools/.../within-family-efficient-weighting.R`): full sibs are
  EXCHANGEABLE → the pair-moment covariance `Σ_g` is compound-symmetric and every pair
  loads equally on β → **optimal weights ARE equal**. CUE/optimal weighting gives no gain
  (s=2/3/4: sd_eq≈sd_eff to 3 digits); estimating the weight only adds finite-sample noise.
- *Between mixed-size families* (`tools/.../within-family-mixedsize-weighting.R`): pooling
  raw within-family pairs (`β̂=ΣA_f/ΣB_f`, each family weight 1) is **ALREADY efficient**
  even with mixed sizes — optimal stratified weighting `w_s∝b_s/v_s` buys nothing
  (rel-eff sd(PAIR)/sd(EFF)=**0.997**). Closed form: rel-eff
  `=[Σ C(s,2)g_s·Σ C(s,2)/g_s]/(Σ C(s,2))²`, `g_s=ν+2(s−2)c` (ν=per-pair var,
  c=cov of two pairs sharing a sib); =1 iff c=0, and the data say rel-eff≈1, i.e.
  **c≪ν empirically** (the per-pair win moment's idiosyncratic variance dominates the
  one-shared-sib covariance). The ONE mistake is normalizing per family (FAMEQ,
  `w=1/C(s,2)`): predicted sd penalty `sqrt((1/6)/(1/10))=1.29`, observed **1.285**.
  **⇒ Equal-pair-weight, family-clustered AR is BOTH the simplest and the efficient
  estimator; do not down-weight large families, do not normalize them.**

**Deep literature scan — DONE (`dev/findings-paper04-litscan.md`).** Retrieval-based scan
(scite + PubMed) over the win×{MR, IV, causal, family, sibling, twin} grid. **The three-way
intersection {win statistic} × {MR/genetic IV} × {within-family} is empty.** Each bordering
literature is occupied but distinct: within-family MR is **scalar-only** (Davies 2019
`10.1093/hmg/ddz204`; Brumpton 2020 `10.1038/s41467-020-17117-4`; Howe 2022
`10.1038/s41588-022-01062-7`); causal win ratio exists **only under no-unmeasured-confounding,
never IV/MR** — incl. the *clustered* "dependent subjects" estimator (arXiv 2212.06676),
whose clustering is a nuisance to balance, **not** a genetic sibling design; GPC/win
statistics (Buyse) are descriptive; and no one does IV/AR inference on a U-statistic outcome
moment. **Position, do not claim a bare "first":** extend within-family MR from scalar to
pairwise/win outcomes; identify under a within-family genetic instrument (targets the
PC-irreducible dynastic/AM confounders the causal-win-ratio designs cannot). Caveat carried:
two paywalled adjacencies (Mann–Whitney causal `insr.12326`; Biometrika net benefit) read at
abstract level only; finish the probabilistic-index deep-scan (Thas/De Neve) before submission.

**Parent–offspring trios — DONE (`tools/.../within-family-trio-winmr.R`).** A genuinely
DIFFERENT structure: the Mendelian instrument is the per-offspring **mid-parent residual**
`g_o = Z_o − (Z_f+Z_m)/2` (mean-zero in every family/subpopulation, since `E[Z_o|parents]`
= mid-parent PRS → clean of dynastic/AM/stratification). Cleaning the instrument at the
**individual** level makes the win comparison a **COMPLETE** U-statistic over all offspring
pairs (vs the sibship's restricted within-family-pairs U-statistic). Key algebra
(antisymmetric `h`): `Σ_{a<b} h(O_a,O_b)(g_a−g_b) = Σ_a g_a W_a`, so
`β̂_trio = Cov(g,w(O))/Cov(g,X)` in `O(N²K)` via the per-subject win-loss kernel.
- **Robust** (bias vs oracle, F=3000): none +0.007, **dyn +0.002**, **am +0.001**, **strat
  +0.009**; the naive offspring-PRS complete-U is catastrophically confounded (dyn −1.03,
  am −0.72, strat −0.35). (PCs would handle strat but NOT dyn/am — same PC-irreducible scope
  as the sibship design.)
- **Consistent** (convergence-in-N ladder, dyn): bias +0.013→+0.011→+0.007 for
  F=1500→3000→6000 (within-trio weak-IV finite-sample bias → 0).
- **Calibrated inference**: a closed-form **complete-U influence-function SE** — needs the
  g-weighted win-score `r(O_a)=Σ_b g_b h_ab` (a second kernel call with `weights=g`), giving
  `ψ_num,a = g_a w(O_a) − r(O_a)`. Wald-CI coverage 0.944–0.978; the IF SE slightly
  under-estimates in confounded scenarios (ratio ~0.90, the omitted higher-order ζ₂ term),
  coverage holds at nominal.

**Two-sample / summary-data within-family form — DONE (`tools/.../within-family-twosample-winmr.R`).**
The programme capstone: chains Papers 02 (IVW) + 03 (AR) into the family design. The family
structure changes ONLY how the per-SNP summary stats are computed — `β_GX,l^wf` = within-family
(sib-centered) exposure-on-dosage slope (what within-sibship GWAS already produce); `δ_l^wf` =
within-family slope of the marginal **win-score phenotype** `w(O_i)=E_j h(O_i,O_j)` on dosage —
then the EXISTING `mrwin_twosample_ivw` / `mrwin_winmr_ar` pool them. Deployable on within-sibship
GWAS summary stats (Howe 2022).
- **Removes dynastic/AM confounding**: naive (pop-centered) IVW is catastrophically biased
  (dyn −0.61, am −0.90); within-family is not.
- **But within-family instruments are intrinsically WEAK** (mean per-SNP F-stat 1.6–4.0 at
  5k–24k subjects), so the within-fam IVW POINT is attenuated and its Wald CI under-covers
  (0.82–0.90) — the textbook weak-IV failure of IVW. The **AR (Paper 03) is calibrated
  (0.91–0.98)** and is the recommended inference; over-ID Q type-I conservative (0.00–0.013).
- **Consistency**: F-ladder (none) — as F grows the per-SNP F-stat grows (1.58→4.02) and the IVW
  bias shrinks (−0.067→−0.040); AR coverage 0.96→1.00 (valid, conservative as the set widens
  honestly under weak IV). **Headline: report the AR confidence set, not the IVW point.** This is
  exactly why the programme needs AR — weak instruments are intrinsic to within-family designs.

**Still UNTESTED / proposed (do not claim):** realism (LD, sib–sib interaction, shared sib
environment); **mixed sibship+trio cohorts** (combine the restricted-U and complete-U moments
efficiently); a two-sample **within-trio** (transmitted/non-transmitted) summary form.

**New algorithm / theory for the paper (Paper 04 candidate):**
1. a **design-restricted (incomplete) U-statistic** over sib pairs (`O(N)` pairs);
2. **family-clustered** influence-function variance + clustered Anderson–Rubin
   (Paper 3 with clustering) — the untested inference piece;
3. extension to **>2 sibs / trios / two-sample within-sibship** summary data;
4. a robustness theorem scoped to dynastic + assortative mating (the PC-irreducible
   confounders) — stratification is a secondary benefit, not the headline.

Position vs scalar within-sibship MR (Brumpton/Davies/Howe et al.): the novelty is
the **win / pairwise-functional** outcome and the statistic-equals-design identity.
**Deep-scan that literature before claiming first** (see caveat).

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

# Findings — deep literature scan for Paper 04 (within-family win-ratio MR)

Status: 2026-06-30. Retrieval-based novelty scan run before any "first" claim, per the
roadmap's honest-novelty caveat (`research-frontier-roadmap.md`). Tools: scite
`search_literature` (210M-paper index, full-text + Smart Citations) and PubMed. Every
DOI below was returned by a tool and verified at abstract level; arXiv preprints are
cited by title + id only (authorship not independently verified — do not attribute).

## The question

Paper 04 sits at a **three-way intersection**:
`{win statistics / pairwise hierarchical-composite outcome}` × `{Mendelian randomization /
genetic IV}` × `{within-family (sibling) design}`. The claim to defend is that this
intersection is **empty** in the literature, and to position the paper against the three
mature single-/double-literatures it borders.

## Verdict: the three-way intersection is empty; each pair is occupied but distinct

**1. Within-family MR exists — but only for SCALAR exposures/outcomes.** This is the
literature to position against, not a strawman.
- Davies, Howe, Brumpton, Havdahl, Evans, Davey Smith (2019), *Hum Mol Genet* 28(R2):
  R170–R179, `10.1093/hmg/ddz204` — the methods review: sibling-pair / parent–offspring
  trio designs overcome population stratification, **dynastic effects, assortative
  mating**, selection/transmission-ratio distortion. Frames the SNP→phenotype association
  for **scalar** phenotypes.
- Brumpton, Sanderson, Heilbron, …, Davies (2020), *Nat Commun* 11:3519,
  `10.1038/s41467-020-17117-4` — within-family MR methods + simulation; the canonical
  source for exactly the three biases this project tested (dynastic, AM, stratification).
  Sib-difference instrument; **scalar** exposures (BMI, height) → **scalar** outcomes
  (diabetes, blood pressure, education). Windmeijer (weak-IV) is a co-author — they care
  about weak instruments, but for the linear ratio, not a U-statistic.
- Howe, Nivard, Morris, … (2022), *Nat Genet* 54:581–592, `10.1038/s41588-022-01062-7`
  — within-sibship GWAS for direct genetic effects; **scalar** marginal associations.
- Davey Smith et al., *Cold Spring Harb Perspect Med*, `10.1101/cshperspect.a039503`
  — "Integrating Family-Based and Mendelian Randomization Designs" (landscape review;
  no pairwise/win outcome proposed).

  → **None analyze a win statistic / pairwise functional / hierarchical composite
  endpoint.** Paper 04's statistic-equals-design identity (the win kernel `h(O_i,O_j)`
  IS the sib contrast) is absent. This is the genuine novelty axis.

**2. Causal win ratio exists — but only under no-unmeasured-confounding, never IV/MR,
never within-family.** Verified by reading abstracts/intros directly:
- "Rethinking the Win Ratio: A Causal Framework for Hierarchical Outcome Analysis",
  arXiv `10.48550/arXiv.2501.16933` — potential-outcomes / stratification-matching
  identification; **no IV, no MR, no family design** (confirmed from text).
- "Causal Inference on Win Ratio for Observational Data with Dependent Subjects",
  arXiv `10.48550/arXiv.2212.06676` — "dependent subjects" = **clustered multi-centre
  data** (patients nested in hospitals); identification by calibrated/propensity weights
  under no-unmeasured-confounding. The clustering is a **nuisance to balance, NOT a
  genetic sibling design**, and there is **no instrument**. Important to draw this
  distinction explicitly — it is the closest "clustered causal win ratio" but is not
  within-family win-MR.
- "Debiased learning of the causal net benefit with censored event time data",
  *Biometrika* `10.1093/biomet/asaf051` — semiparametric/debiased causal net benefit
  (relevant to the efficiency branch **P2**, not P4); **no IV**. (Abstract-level only.)

**3. Win statistics / generalized pairwise comparisons (Buyse) are descriptive
clinical-trial tools.** Overview `10.1002/bimj.202100354`; correlated prioritized
outcomes / net benefit `10.1002/sim.8788`; encyclopedia entry
`10.1002/9781118445112.stat08224`. Original net-benefit/GPC: Buyse (2010), *Stat Med*,
`10.1002/sim.3923` (already in the project bib). **No causal-IV layer.**

**4. Adjacent nonparametric causal effects — Mann–Whitney / probabilistic-index.**
"Estimating Mann–Whitney-type Causal Effects", *Int Stat Rev* `10.1111/insr.12326`
(full text not retrieved — paywall 403; cite as adjacency, do not assert its internals).
These target a Mann–Whitney/concordance causal estimand under covariate adjustment, not
under genetic IV, and not within-family. Position P1 against this literature; for P4 it
is a distant cousin.

**5. Degree-2 / pairwise IV machinery.** "Pairwise Valid Instruments" arXiv
`10.48550/arXiv.2203.08050` concerns **instrument validity assessed pairwise**, not an
IV-identified pairwise *outcome* functional. No retrieved work does Anderson–Rubin /
weak-IV-robust inference on a **U-statistic moment** (the Paper 03 / P3 novelty) nor an
efficiency bound for an IV-identified ratio of U-statistics (P2). This corroborates the
programme's unifying thesis: *IV identification/efficiency/robust-inference has been
built for degree-1 (mean) functionals and never lifted to degree-2 (pairwise) ones.*

## How to position Paper 04 (do NOT claim a bare "first")

- **vs within-family MR (Davies/Brumpton/Howe):** "We extend the within-family MR design
  from scalar exposures/outcomes to **pairwise hierarchical-composite (win) outcomes**,
  where the win kernel and the sibling contrast are the same degree-2 object." Inherit
  their bias taxonomy (dynastic, AM, stratification) — the project's probes already map
  onto it exactly.
- **vs causal win ratio (arXiv 2025/2022, Biometrika):** "Prior causal win-ratio work
  identifies under no-unmeasured-confounding (including clustered observational data);
  we identify under a **genetic instrument within families**, which targets the
  PC-irreducible confounders (dynastic, assortative mating) those designs cannot remove."
- **vs GPC / win statistics (Buyse):** descriptive estimands; we add the **causal-IV**
  layer.

## Honest caveat (carry into the manuscript)

- Retrieval-based: strong (no win-statistic MR or within-family win-MR surfaced across
  the win×{MR,IV,causal,family,sibling,twin} query grid; the "win ratio MR" string still
  returns only the author-surname "Win" collisions noted earlier), but not dispositive.
- Two adjacent full texts were **not read** (Mann–Whitney causal `insr.12326`; the
  Biometrika net-benefit) — paywalled. Do not assert their internals; both are
  abstract-/venue-level judged non-IV.
- Before submission, still run the deep-scan under: *probabilistic index models* (Thas,
  De Neve), *concordance/Mann–Whitney causal effects*, *prioritized outcomes*, and check
  whether any 2025–26 preprint has moved into the within-family win space.
- Cite real DOIs only; never attribute the arXiv preprints to unverified authors.

## Verified reference list (for the Paper 04 bib)

- Davies NM, Howe LJ, Brumpton B, Havdahl A, Evans DM, Davey Smith G. Within family
  Mendelian randomization studies. *Hum Mol Genet* 2019;28(R2):R170–R179.
  https://doi.org/10.1093/hmg/ddz204
- Brumpton B, Sanderson E, Heilbron K, et al. Avoiding dynastic, assortative mating, and
  population stratification biases in Mendelian randomization through within-family
  analyses. *Nat Commun* 2020;11:3519. https://doi.org/10.1038/s41467-020-17117-4
- Howe LJ, Nivard MG, Morris TT, et al. Within-sibship genome-wide association analyses
  decrease bias in estimates of direct genetic effects. *Nat Genet* 2022;54:581–592.
  https://doi.org/10.1038/s41588-022-01062-7
- Rethinking the Win Ratio: A Causal Framework for Hierarchical Outcome Analysis. arXiv
  2501.16933. https://doi.org/10.48550/arXiv.2501.16933
- Causal Inference on Win Ratio for Observational Data with Dependent Subjects. arXiv
  2212.06676. https://doi.org/10.48550/arXiv.2212.06676
- Debiased learning of the causal net benefit with censored event time data.
  *Biometrika* 2025. https://doi.org/10.1093/biomet/asaf051
- Buyse M. Generalized pairwise comparisons of prioritized outcomes in the two-sample
  problem. *Stat Med* 2010;29:3245–3257. https://doi.org/10.1002/sim.3923
- Péron J, Buyse M, et al. Generalized pairwise comparisons for censored data: an
  overview. *Biom J* 2021. https://doi.org/10.1002/bimj.202100354
- Estimating Mann–Whitney-type Causal Effects. *Int Stat Rev*.
  https://doi.org/10.1111/insr.12326  (adjacency; full text not retrieved)
- Pairwise Valid Instruments. arXiv 2203.08050.
  https://doi.org/10.48550/arXiv.2203.08050  (adjacency; instrument validity, not outcome)

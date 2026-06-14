# WP18 — Theory Foundation & Manuscript

Status block: continuous; updated whenever another WP lands a result.
Branch: `C-wp18` (from `C`); merges frequently and in small pieces.
Depends on: results from WP13–WP17 (no claim written before its test passes).
Blocks: external release / arXiv submission.

This work package turns the implementation into a **citable methodological
contribution**. It has two surfaces:

- **Public, committed** (`inst/spec/`): this file — the proof-obligation
  ledger and a public-safe statement of each result. Anyone reading the repo
  sees *what* is claimed and *that* it is tested, not the in-progress prose.
- **Private, gitignored** (`manuscript/`, `drafts/`, `*.tex`, `*.qmd`): the
  actual paper draft, figures, and derivations. Kept out of the public repo by
  `.gitignore` until submission.

Hard rule (from the roadmap §4.5.6): **no manuscript claim is written until the
corresponding implementation test passes.** Theory never gets ahead of code.

---

## 1. Manuscript skeleton (`manuscript/` — gitignored)

1. **Introduction** — hierarchical composite endpoints; causal win statistics;
   the scalability and `D`-sensitivity gaps Phase II closes.
2. **Estimand** — cCWR / DS-CWR, the continuous ISG (WP15), identification under
   the MR assumptions, non-collapsibility caveat.
3. **Computation** — the fast kernel (WP13), bootstrap reuse and incremental
   re-stratification (WP14); complexity theorems.
4. **Inference** — analytic influence-function variance + GWAS delta-method
   (WP17); relation to the multiplier bootstrap; Q statistic asymptotics
   (existing `q_statistic_asymptotics.py`).
5. **Stratification** — doubly-ranked strata (WP16) and the assumption
   trade-off vs residual ranking.
6. **Simulations** — coverage, type-I, power, scalability (WP19).
7. **Application** — template only unless a real cohort is authorised.
8. **Discussion** — limitations, pleiotropy sensitivity, scope of claims.

---

## 2. Proof-obligation ledger

Each row is discharged only when its test is green. Keep the status current.

| # | Claim | Source WP | Status |
|---|---|---|---|
| L-C1 | Single-endpoint win/loss = weighted dominance count; `Θ(N log N)` | WP13 | pending |
| T-C2 | Hierarchical win/loss = bounded-dimension orthogonal range counting; `Θ(N log^{K-1} N)`; exact vs dense kernel | WP13 | pending |
| P-C3 | Incremental re-stratification cost `O(I log N)`, `I` = inversions; exact vs full re-query | WP14 | pending |
| D-M1 | Continuous ISG: definition, identification, consistency; boxcar limit = decile estimator; U-process CLT | WP15 | pending |
| S-M2 | Doubly-ranked strata validity under rank-preservation; weaker than linearity/homogeneity | WP16 | pending |
| V-M3 | Influence function of stratified DS-CWR; delta-method through cutpoints; analytic = bootstrap in the limit | WP17 | pending |

---

## 3. Tasks

- **T1** — Stand up `manuscript/` (gitignored) with the skeleton above and a
  build (LaTeX or Quarto). Confirm `.gitignore` excludes it.
- **T2** — For each landed WP result, write the corresponding section and flip
  the ledger row to "discharged (test: <path>)".
- **T3** — Keep a public-safe one-paragraph statement of each result in this
  file (no unproven claims, no private data).
- **T4** — Pre-submission pass: confirm every ledger row is discharged and every
  number in the manuscript traces to a committed test or benchmark artifact.

---

## 4. Definition of Done

Manuscript builds; every ledger row discharged with a test reference; public
statements in this file match the implemented behaviour; citations verified
(see each WP's reference list — confirm exact bibliographic details rather than
relying on memory).

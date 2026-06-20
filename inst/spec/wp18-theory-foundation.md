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

Headline contributions (the validated three): (i) a subquadratic + compiled win/
loss kernel that takes the estimator to biobank scale; (ii) an analytic
influence-function variance that replaces the multiplier-bootstrap loop; (iii) a
calibration finding — the original delta-method ratio CI is mis-calibrated and
**Fieller** is the correct, now-primary interval. Plus two instructive **negative
findings** that carry the paper's methodological thesis: the continuous-ISG
generalisation (M1) and doubly-ranked stratification (M2) each pass every
*internal* consistency check yet fail *external* calibration — the second
(type-I ~1.0) because doubly-ranked balances the instrument across strata and so
is incompatible with the between-stratum DS-CWR contrast.

1. **Introduction** — hierarchical composite endpoints; causal win statistics;
   the scalability gap and the (subtler) inference-calibration gap Phase II
   closes.
2. **Estimand** — cCWR / DS-CWR, the instrument-standardised gradient (ISG),
   identification under the MR assumptions, non-collapsibility caveat. The
   continuous-ISG generalisation (WP15) is reported as an explored **negative
   result** (a ratio-variance obstruction), motivating the discrete decile +
   sensitivity-over-`D` recommendation.
3. **Computation** — the fast kernel (WP13), bootstrap reuse and incremental
   re-stratification (WP14); complexity theorems and the Rcpp wall-clock.
4. **Inference** — analytic influence-function variance + exact GWAS-only
   resample (WP17); relation to the multiplier bootstrap; **calibration: the
   bivariate-Delta vs Fieller ratio CI** (WP19, the R2 finding) — the headline
   correctness result.
5. **Stratification** — doubly-ranked strata (WP16) reported as a **negative
   finding**: it is wired end to end and mechanically correct, but balancing the
   instrument across strata makes the between-stratum DS-CWR contrast
   confounder-driven (external type-I ~1.0), so it is incompatible with this
   estimand (it suits a within-stratum LACE estimator instead).
6. **Simulations** — type-I, coverage, scalability, and the calibration grid
   across bootstrap/analytic/GWAS-uncertainty/IPTW/doubly-ranked paths (WP19).
7. **Application** — template only unless a real cohort is authorised.
8. **Discussion** — limitations (weak-instrument fragility of the ratio
   estimand; first-order IPTW correction), pleiotropy sensitivity, scope of
   claims; the methodological lesson that internal consistency (fast == dense,
   analytic == bootstrap) cannot detect a shared calibration defect.

---

## 2. Proof-obligation ledger

Each row is discharged only when its test is green. Keep the status current.

| # | Claim | Source WP | Status |
|---|---|---|---|
| L-C1 | Single-endpoint win/loss = weighted dominance count; `Θ(N log N)` | WP13 | discharged (test: `tests/python/test_kernel_fast.py`, `tests/testthat/test-kernel-fast.R`; benchmark: `benchmark-results.md`) |
| T-C2 | Hierarchical win/loss via tie-split + Fenwick 2D / CDQ 3D dominance counting; `Θ(N log^{K-1} N)`, exact vs dense kernel | WP13 | discharged for K=2,3 (tests: `test_kernel_fast.py::test_pair_2d_parity_random`/`test_pair_3d_parity_random`, `test-kernel-fast.R`); K≥4 open |
| P-C3 | Incremental re-stratification cost `O(I log N)`, `I` = inversions; exact vs full re-query | WP14 | pending |
| D-M1 | Continuous ISG: definition, identification, consistency; boxcar limit = decile estimator; U-process CLT | WP15 | **explored → not pursued (negative finding)**: boxcar-limit = decile reduction verified exactly (test: `test-continuous-isg.R`), but external validation shows the continuous estimator is 2.5–3.7× noisier than the decile and does **not** reduce `D`-sensitivity — the ISG is a ratio `logθ/ΔX` and finer smoothing shrinks `ΔX`, inflating the ratio variance. Decile + Fieller is the foundation. See `wp15-continuous-isg.md`. |
| S-M2 | Doubly-ranked strata validity under rank-preservation; weaker than linearity/homogeneity | WP16 | **explored → INCOMPATIBLE with the cCWR estimand (negative finding)**: the helper (Tian/Burgess 2023) is correct and wired end to end (fitted strata match across dense/sparse/fast/analytic, test: `test-doubly-ranked.R`), but external calibration gives **type-I ~1.000** because doubly-ranked *balances the instrument across strata* (its within-stratum-LACE selling point), which destroys the *between*-stratum IV contrast the DS-CWR estimand needs — the contrast becomes confounder-driven. Kept runnable with a `doubly_ranked_invalid` warning; default stays `prs_rank`. See `validation-findings.md`, `wp16-doubly-ranked-strata.md`. |
| V-M3 | Variance decomposition `cov_u = Σ_sampling + Σ_gwas`; Σ_sampling analytic (`CᵀC`), Σ_gwas exact GWAS-only resample (point estimate piecewise-constant in β → no pointwise gradient); IPTW via weighted IF (estimated-weights correction omitted, ~1–3%); analytic = bootstrap | WP17 | discharged: `adjustment="none"` exact (`se` ratio 0.997–1.009 incl. σ_β∈{0.05,0.15,0.30}); IPTW first-order (`se` ratio 0.985–1.026); wired as `mrwin(inference="analytic")` (test: `test-analytic-variance.R`) |
| C-R2 | Calibration: the package's original bivariate-Delta ratio CI is mis-calibrated (~2× too wide → structural type-I 0.000, coverage 0.995); the **Fieller** construction on the same `cov_u` is correctly calibrated (type-I 0.055, coverage 0.955) and is now the primary reported interval | WP19 | discharged: found by **external** calibration simulation, not internal consistency; Fieller made primary across `print`/`summary`/`tidy`/`mrwin_report` (test: `test-fieller-primary.R`); widened grid confirms no over-rejection across the **valid** paths — bootstrap/analytic, GWAS-uncertainty, and IPTW on both engines (bootstrap×IPTW 0.017, analytic×IPTW 0.060) — `prs_rank` stratification (`validation-findings.md`, harness `R/validate_calibration.R`). (The same external grid is what exposed the doubly-ranked incompatibility, S-M2.) |

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

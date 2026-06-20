# WP16 — Doubly-Ranked Stratification (M2)

Status block: `T1 [done] T2 [done] T3 [partial: unit + balance + end-to-end wiring tests done]`
Branch: developed on `C-wp13` alongside S1/S2 (independent of the kernel work).

Progress (2026-06-14, R-verified): `mrwin_doubly_ranked_strata(prs, exposure,
n_strata)` implemented in `R/strata.R`, exported, matching the
`mrwin_prs_strata` return contract. Tests `tests/testthat/test-doubly-ranked.R`:
hand-computed toy example, strictly-increasing mean exposure across strata with
balanced instrument means, exact balanced counts for full blocks, deterministic
remainder handling, input validation. Full suite 79 groups, 0 failures.

**T2 done (2026-06-20, R 4.3.3):** `stratification = c("prs_rank",
"doubly_ranked")` is threaded end to end. A single internal dispatch
`.mrwin_assign_strata(G, beta, X, n_strata, stratification)` (in `R/strata.R`)
returns the `mrwin_prs_strata` contract and is called at **every**
(re-)stratification point — `mrwin_estimate`, `mrwin_sparse_estimate`,
`mrwin_multiplier_bootstrap`, `mrwin_sparse_bootstrap`,
`mrwin_analytic_bootstrap`, `mrwin_analytic_inference`,
`mrwin_gwas_resample_covariance` — so the point estimate, each bootstrap
multiplier draw, and each GWAS-resample draw all honour the chosen scheme.
Doubly-ranked requires the exposure `X` at every re-stratification (the
incremental WP14 path therefore re-derives pre-strata from `X`; it does not
reuse a PRS-only bin map). `mrwin_controls(stratification = ...)` exposes it;
default stays `"prs_rank"` so existing behaviour is unchanged. End-to-end test
in `test-doubly-ranked.R` confirms the fitted strata equal the direct helper
across the dense / sparse / fast backends and the analytic inference path.
Depends on: WP4 stratification (`mrwin_prs_strata`).
Blocks: nothing hard; strengthens WP15 and WP19.

The current package stratifies on raw PRS rank — the "residual method" style,
which relies on **linearity and homogeneity** between instrument and exposure to
form valid strata. When those assumptions fail, the instrumental-variable
assumptions can fail *within strata* even if they hold in the population. This WP
adds **doubly-ranked stratification**, which forms strata with distinct exposure
levels without those parametric assumptions.

---

## 1. Method

Doubly-ranked stratification (Tian, Burgess, et al., 2023):

1. Rank individuals by the instrument (PRS).
2. Form pre-strata of size `k` by this ranking.
3. Within each pre-stratum, rank by the **exposure**; assign the `j`-th ranked
   individual to final stratum `j`.

This yields strata that differ in mean exposure while preserving the instrument
ordering, valid under a weaker "rank-preserving" assumption rather than
linearity/homogeneity.

---

## 2. Tasks

- **T1** — Implement `mrwin_doubly_ranked_strata(prs, exposure, n_strata)`
  returning the same contract as `mrwin_prs_strata` (1-indexed labels,
  deterministic tie handling). Add `stratification = c("prs_rank",
  "doubly_ranked")` to `mrwin_controls()`.
- **T2** — Wire it through `mrwin_estimate` / `mrwin_sparse_bootstrap` /
  `mrwin()`; ensure bootstrap re-stratification uses the chosen method
  consistently (and remains compatible with WP14 incremental updates — note:
  doubly-ranked changes how `S*` maps to bins, so the incremental path must
  re-derive pre-strata; document the interaction).
- **T3** — Scenario validation: on a non-linear exposure-instrument DGP
  (extend WP8), show doubly-ranked strata maintain within-stratum IV validity
  where residual ranking does not. Deterministic unit tests for the assignment
  on toy vectors.

---

## 3. Acceptance criteria

1. Assignment matches a hand-computed toy example exactly; ties deterministic.
2. On the linear DGP, doubly-ranked and PRS-rank give compatible DS-CWR
   (no regression); on a non-linear DGP, doubly-ranked shows reduced
   within-stratum IV violation.
3. PRS-rank remains the default; doubly-ranked is opt-in and documented.

## 4. Definition of Done

§3 satisfied; status block all `[done]`; WP18 "Stratification" section cites the
method and states the assumption trade-off; checkpoint updated.

## 5. Reference

- Tian, H., Mason, A. M., Liu, C., Burgess, S. (2023). Relaxing parametric
  assumptions for non-linear Mendelian randomization using a doubly-ranked
  stratification method. *PLOS Genetics* 19(6): e1010823.

# WP19 — Scalability Validation, Benchmarks & Release

Status block: `T1 [todo] T2 [todo] T3 [todo] T4 [todo] T5 [todo]`
Branch: `C-wp19` (from `C`); integrator — merges last.
Depends on: WP13, WP14 (hard); WP15, WP16, WP17 (for full validation).
Blocks: any external release / arXiv submission.

This is the integration and gate WP. It proves the speedups are real, proves the
fast paths match the reference, and flips defaults only when parity is green.

---

## 1. Parity gate (Track A correctness)

A single harness `mrwin_verify_fast_dense_parity()` (extend the existing
`mrwin_verify_sparse_dense_parity`) that asserts, across randomised cohorts and
seeds:

- `mrwin_fast_pair_win_loss` == dense reference (integer-exact unweighted;
  `1e-10` weighted), `K ∈ {1,2,3}`, including the adversarial tie / censoring
  cases from WP13 T5.
- Same-seed fast bootstrap == pre-WP14 bootstrap output within `1e-10`.
- Incremental re-stratification == full re-query within `1e-10`.

The fast backend may become a default **only after this gate is green in CI.**

---

## 2. Scaling benchmark (Track A performance)

Extend `R/benchmark.R::mrwin_benchmark`:

- Grid: `N ∈ {1e3, 3e3, 1e4, 3e4, 1e5}`, `K ∈ {1,3}`, `D ∈ {5,10}`,
  `B ∈ {0, 200}`.
- Compare: dense reference vs fast (WP13) vs fast+reuse (WP14).
- Record wall-clock and peak memory; fit empirical exponent `t ∝ N^p` and assert
  `p < 1.3` for the fast paths vs `p ≈ 2` for the reference.
- Output a committed summary table (numbers only; no environment-identifying
  details) under `inst/spec/benchmark-results.md`; raw artifacts stay gitignored.

---

## 3. Statistical validation (Track B)

Run the post-release validation plan from `project-checkpoint.md` against the
new estimators:

- Null type-I ≈ 0.05; valid-IV power; pleiotropy detection; weak-instrument
  Fieller behaviour; discordant-component warnings.
- Continuous-ISG coverage (WP15) and `D`-sensitivity reduction.
- Analytic-variance coverage and agreement (WP17).
- Doubly-ranked within-stratum validity (WP16).

---

## 4. Tasks

- **T1** — Parity harness + CI wiring (must pass before any default flip).
- **T2** — Scaling benchmark + committed summary table; assert exponent bounds.
- **T3** — Statistical validation grids for Track B estimators.
- **T4** — Flip defaults: once T1–T3 are green, make `backend = "fast"` and the
  bootstrap reuse the default; keep dense/quadratic path available as
  `backend = "dense"` for debugging and small-`N` exactness.
- **T5** — Release pass: `R CMD check --as-cran`, `pytest`, coverage ≥ 80%,
  update `project-checkpoint.md`, `work-package-roadmap.md`, and READMEs to
  reflect the shipped fast backend and continuous estimator.

---

## 5. Acceptance criteria

1. Parity gate green in CI.
2. Fast-path empirical exponent `< 1.3`; reference `≈ 2`; documented table.
3. All Track B validation targets met within MC error.
4. `R CMD check` and `pytest` clean; coverage ≥ 80%.

## 6. Definition of Done

§5 satisfied; defaults flipped; docs updated; checkpoint marks Phase II
complete and lists the shipped capabilities.

# Phase II Roadmap — Scalable & Continuous Causal Win Statistics

Status: canonical roadmap for the Phase II program (WP13–WP19).
Created: 2026-06-14.
Owner of record: package maintainer.
Audience: autonomous contributor agents and human reviewers.

This document is the single source of truth for the second program of work on
`mrwin`. Phase I (WP0–WP12, see `work-package-roadmap.md`) delivered a correct,
well-tested, but **computationally quadratic** package. Phase II makes the method
**biobank-scale** and **methodologically stronger**, and produces the **theory
manuscript** that turns the implementation into a citable contribution.

Read this file first. Then read the per-work-package spec it points to before
touching code.

> **Status (2026-06-18):** WP13 (fast kernel C1+C2, K=1/2/3) and its **Rcpp port**
> are DONE — the `O(N²)` goal is solved and biobank-scale. M2 core is done. The
> remaining directions are evaluated and re-prioritised in
> **`phase2-direction.md`**: recommended path **M3 → M1 → WP19 → release**;
> WP14 and K≥4 deferred. ~23 steps remain.

---

## 1. Why Phase II exists

### 1.1 The bottleneck, stated precisely

The estimator compares all cross-pairs between adjacent polygenic-risk-score
(PRS) strata. With `N` individuals, `K` priority endpoints, `D` strata, and a
multiplier bootstrap of `B` iterations, the dominant cost is:

```
Total ≈ B · [ N·M  (recompute PRS S* = G β*)
            + N log N  (re-sort into strata)
            + Σ_d N_d · N_{d-1} · K   (win/loss pair sweep)  ]

with  Σ_d N_d · N_{d-1} ≈ N² / D.
```

The pair sweep term `B · N²·K / D` dominates and is **quadratic in N**. The
existing "sparse" backend only avoids materialising the `N×N` matrix; it still
performs `Θ(N²/D)` comparisons. At `N = 337,000`, `B = 200`, `K = 3`, `D = 10`
this is `~7×10¹⁴` comparisons — infeasible. No constant-factor optimisation
(vectorisation, Rcpp, parallelism) changes the **order**; Phase II changes the
order.

### 1.2 The structural fact that makes acceleration possible

The hierarchical kernel `h(i,j) ∈ {-1,0,+1}` depends **only** on the observed
`(time, status)` of the pair. It does **not** depend on `β`, on the PRS, or on
stratum assignment. Across bootstrap iterations the kernel is invariant; only
(a) the **stratum labels** (rank-bands of `S*`) and (b) the per-individual
**weights** `ξ_i` change. Every acceleration in Phase II exploits this
invariance: build kernel-summary structures once on the fixed outcome
coordinates, then re-query them cheaply per iteration.

### 1.3 Relationship to prior art

A near-linear two-sample win-ratio computation already exists in the wild
(rank-counting via a Fenwick / binary-indexed tree; see WP13 references). That
prior art is **not** the contribution: it handles two fixed groups, no weights,
no instrument, no bootstrap, no stratification. Phase II's novelty is the
**dynamic-stratum, weighted, bootstrap-reused** version (WP13–WP14) plus a
**continuous estimator** that removes the arbitrary stratum count (WP15) and an
**analytic variance** that removes most of the bootstrap (WP17). These are the
publishable units.

---

## 2. The contribution ladder

Phase II is organised as two tracks that share infrastructure.

### Track A — Algorithmic acceleration (correctness-preserving)

| ID | Name | Outcome |
|---|---|---|
| C1 | Single-endpoint fast win/loss | `Θ(N²/D)` → `Θ(N log N)` per iteration |
| C2 | Hierarchical (`K`-priority) fast win/loss | `Θ(N log^{K-1} N)` per iteration via multidimensional dominance counting |
| C3 | Bootstrap reuse + incremental re-stratification | reuse fixed-kernel structures across `B`; exploit `β* ≈ β` to update by rank inversions, not from scratch |

Track A must produce **bit-for-bit / tolerance-exact parity** with the existing
dense estimator. It changes speed, never answers.

### Track B — Methodological strengthening (changes/adds estimands)

| ID | Name | Outcome |
|---|---|---|
| M1 | Continuous kernel-smoothed ISG | removes the discrete `D` choice; deciles become the boxcar special case; cleaner U-process asymptotics |
| M2 | Doubly-ranked stratification | strata valid without the linearity/homogeneity assumption of residual ranking |
| M3 | Analytic influence-function variance + GWAS delta-method | cut `B` by ~10× or eliminate it; faster CIs; connects to existing `q_statistic_asymptotics` work |

Track B changes results by design and therefore requires its own statistical
validation (coverage, calibration), not just parity.

### Work-package mapping

| WP | Track items | Spec file |
|---|---|---|
| WP13 | C1, C2 | `wp13-fast-kernel.md` |
| WP14 | C3 | `wp14-bootstrap-acceleration.md` |
| WP15 | M1 | `wp15-continuous-isg.md` |
| WP16 | M2 | `wp16-doubly-ranked-strata.md` |
| WP17 | M3 | `wp17-analytic-variance.md` |
| WP18 | Theory foundation / manuscript | `wp18-theory-foundation.md` |
| WP19 | Scalability validation, benchmarks, parity gates, release | `wp19-scalability-validation.md` |

---

## 3. Dependency graph (execution order)

```
WP13 (fast kernel)  ──► WP14 (bootstrap reuse) ──┐
        │                                        ├─► WP19 (scale validation + release)
        ├──────────────► WP15 (continuous ISG) ──┤
        │                                        │
WP16 (doubly-ranked) ────────────────────────────┤
WP17 (analytic variance) ────────────────────────┘

WP18 (theory) runs in parallel and is updated as each WP lands a result.
```

Rules:

1. **WP13 is the critical path.** Nothing in Track A or M1 is fast until WP13
   lands. Start here.
2. WP14 requires WP13's prebuilt structures.
3. WP15 (continuous ISG) reuses WP13's sweep but is statistically independent;
   it can be prototyped on the dense backend first, then accelerated.
4. WP16 and WP17 are independent of WP13 and may proceed in parallel by separate
   agents.
5. WP18 (theory) is continuous: every landed lemma/benchmark updates the
   manuscript. WP18 must never get ahead of an implemented, tested result —
   theory claims are only written once the corresponding test passes.
6. WP19 is the integrator and the release gate; it depends on everything.

---

## 4. Multi-agent execution protocol

This program is designed to be executed by multiple autonomous agents working
in parallel. Follow these rules exactly.

### 4.1 Branching

- The integration branch for Phase II is `C`.
- Each work package is developed on a child branch named `C-wp13`, `C-wp14`, …
  branched from `C`.
- Sub-tasks within a WP that warrant isolation use `C-wp13-<short-slug>`.
- Never push to `main`. Never force-push a shared branch.
- Open work merges back into `C` only when its WP "Definition of Done" (§4.5)
  is satisfied.

### 4.2 Task granularity

Each WP spec contains a numbered task list (`T1`, `T2`, …). A task is the unit
an agent claims. A task must be:

- independently testable,
- completable in a single focused session,
- accompanied by its own tests in the same commit.

### 4.3 Claiming and status

- Maintain a status block at the top of each WP spec: `T1 [done] T2 [in-progress: <branch>] …`.
- Refresh it in the same commit that changes task state. Do not rely on memory
  or external trackers; the spec file is the ledger.

### 4.4 Commit hygiene

- One logical change per commit; tests included.
- Commit messages: imperative subject, body explaining *why*. Reference the WP
  and task ID, e.g. `WP13 T2: linear-time single-endpoint win/loss sweep`.
- Do not commit generated artifacts (see `.gitignore`); do not commit data.
- Do not introduce tool-, vendor-, or assistant-identifying strings anywhere in
  the repository (code, comments, docs, commit messages). Use neutral wording.

### 4.5 Definition of Done (per WP)

A work package is done when **all** hold:

1. Every task `Tn` is `[done]`.
2. New code has unit tests; `R CMD check` and `pytest` pass locally.
3. **Track A only:** parity test passes — new fast path equals the WP3/WP4 dense
   reference within the tolerance the spec names (default: exact integer kernel
   equality; `1e-10` for weighted sums).
4. **Track B only:** statistical validation passes — the spec's coverage /
   calibration / consistency check meets its stated target.
5. A complexity/benchmark note is recorded (WP19 harness) showing the intended
   scaling on at least three `N`.
6. The corresponding manuscript section (WP18) is updated or explicitly marked
   "pending result".
7. The WP spec's status block reads all-done and the checkpoint
   (`project-checkpoint.md`) is updated.

### 4.6 Reuse-first rule

Before writing a new function, search `R/` and `python/p1_engine_v5/` for an
existing primitive. The fast paths must live behind the **existing public API**
(`mrwin_pair_win_loss`, `mrwin_sparse_bootstrap`, `mrwin()` with
`backend = "fast"`) so that swapping backends never changes user-visible
results, only runtime. Add new backends as opt-in; never silently change a
default until WP19 parity is green.

---

## 5. Theory ↔ code correspondence

Every implemented claim has a written counterpart; every written claim has a
test. WP18 maintains the manuscript; the table below is the contract.

| Implemented in | Proven / written in | Verified by |
|---|---|---|
| C1 linear sweep (WP13) | Lemma: win-count = weighted dominance count; complexity `Θ(N log N)` | parity test + scaling benchmark |
| C2 dominance counting (WP13) | Theorem: hierarchical win/loss reduces to bounded-dimension orthogonal range counting; complexity `Θ(N log^{K-1} N)` | parity test on `K=2,3` |
| C3 incremental update (WP14) | Proposition: per-iteration cost `O(I log N)` where `I` = rank inversions between `S` and `S*` | parity + inversion-count benchmark |
| M1 continuous ISG (WP15) | Definition + consistency of the kernel-smoothed gradient; deciles as boxcar limit; U-process CLT | coverage simulation; agreement with decile estimator as bandwidth → band width |
| M2 doubly-ranked (WP16) | Stratum-validity argument under rank-preservation; comparison to residual method | scenario simulation showing IV validity within strata |
| M3 analytic variance (WP17) | Influence function / Hájek projection of DS-CWR; delta-method through stratum cutpoints | variance agreement vs bootstrap; coverage |

---

## 6. Risk register

| Risk | Where | Mitigation |
|---|---|---|
| Dominance-counting constant blows up with `K` | C2 | Cap exact path at small `K`; provide thresholded/approximate fallback for large `K`; document the crossover |
| Censoring makes the order a *partial* order, breaking naive dominance reductions | C1/C2 | Treat tie/incomparable regions explicitly via the gating indicators `D_{·,k}`; validate against dense kernel on adversarial tie cases |
| Large GWAS uncertainty destroys the `β*≈β` assumption | C3 | Detect inversion count per iteration; fall back to full rebuild when it exceeds a threshold |
| Continuous estimator introduces bandwidth-selection instability | M1 | Cross-validated bandwidth + report sensitivity curve; keep decile estimator as default until coverage validated |
| Analytic variance under-covers in finite samples | M3 | Keep bootstrap available; use analytic variance as control variate first, replacement only after coverage passes |
| Theory outruns implementation | WP18 | Hard rule §4.5.6: no manuscript claim without a passing test |

---

## 7. Glossary (shared vocabulary for all agents)

- **cCWR / DS-CWR**: continuous / dose-standardised Causal Win Ratio, the
  package estimand. `DS-CWR = exp(δ_GLS)`.
- **Kernel `h(i,j)`**: hierarchical pairwise comparison in `{-1,0,+1}`,
  antisymmetric, fixed per cohort.
- **ISG**: instrument-standardised gradient, `log θ_d / ΔX_d`.
- **Win/loss sweep**: the cross-pair aggregation that Phase II accelerates.
- **Dominance counting**: counting pairs where one point exceeds another in all
  coordinates of a partial order; the geometric form of the win count.
- **Fenwick tree / BIT**: prefix-sum structure supporting `O(log N)` point
  update and prefix query; the workhorse of C1/C2.
- **Multiplier bootstrap**: resampling by drawing `β* ~ N(β, Σ_GWAS)` and
  `ξ_i ~ Exp(1)` and recomputing the pipeline.
- **Parity**: numerical equality between a fast path and the dense reference.

---

## 8. Where to start (for a fresh agent)

The **ordered, sprint-by-sprint execution plan** with exact function signatures,
algorithm pseudocode, parity tests, and the benchmark protocol is
`inst/spec/implementation-plan.md`. It sequences the work as
S1=C1 → S2=C2 → S3=M2 → S4=(M3+M1). Start there.

1. Read this file, then `implementation-plan.md`.
2. Read `wp13-fast-kernel.md`. Claim sprint S1 / task `T1`.
3. Branch `C-wp13` from `C`.
4. Implement on the dense reference first to prove correctness, then optimise.
5. Land with parity test green and a scaling benchmark recorded.
6. Update the WP status block and `project-checkpoint.md`.

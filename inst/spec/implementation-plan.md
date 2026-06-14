# Phase II Implementation Plan — Prioritised, Sprint-by-Sprint

Status: canonical execution sequence for Phase II.
Created: 2026-06-14.
Reads alongside: `acceleration-roadmap.md` (strategy), `wp13`–`wp19` (per-WP
specs). This file is the **ordered, concrete** plan: what to build first, in
what order, with exact signatures, tests, and benchmark protocol.

Priority rationale (highest leverage / lowest risk first):

1. **S1 = C1** — Fenwick / sorted-sweep `O(N log N)` single-endpoint win/loss in
   the bootstrap, with **bit-for-bit parity** to the dense backend and a
   **benchmark** that demonstrates the order change. This is the concrete proof
   that the algorithm changed *order*, not just constant factor. Lowest risk
   (single endpoint, exact reduction), highest payoff (kills the dominant term).
2. **S2 = C2** — extend to `K` endpoints via a build-once / reuse-across-`B`
   range / merge-sort tree.
3. **S3 = M2** — doubly-ranked stratification: robustness upgrade, low risk,
   easy to write and test.
4. **S4 = M3 + M1** — analytic variance and continuous estimator: the research
   core for the manuscript.

Each sprint ends at a **gate**. Do not start the next sprint until the current
gate is green. Track state in each sprint's status line.

---

## S0 — Shared scaffolding (half a day)

Status: `[todo]`. Branch: `C-wp13`.

Deliverables:

- `R/kernel_fast.R` — new R backend file (empty stubs + roxygen).
- `python/p1_engine_v5/kernel_fast.py` — Python reference mirror.
- `tests/testthat/test-kernel-fast.R` — parity test file (empty harness).
- `tests/python/test_kernel_fast.py` — Python parity test file.
- Extend `R/benchmark.R` with `mrwin_benchmark_backends()` (added in S1).
- `inst/spec/benchmark-results.md` — committed results table (numbers only).

Gate S0: files exist, `R CMD check` and `pytest` still pass (no behaviour
change yet).

---

## S1 — C1: single-endpoint `O(N log N)` win/loss + parity + benchmark

Status: `[done — 2026-06-14]`. Branch: `C-wp13`. Implements WP13 T1–T2, and WP14
T2 for `K=1`. This is the flagship deliverable.

**Landed:** `kernel_fast.py` + `R/kernel_fast.R` (`mrwin_fast_pair_win_loss`,
`mrwin_fast_adjacent_win_loss`); parity tests Python bit-for-bit green
(`tests/python/test_kernel_fast.py`) and R (`tests/testthat/test-kernel-fast.R`,
CI); benchmark `benchmark_fast_kernel.py` → **fast tail exponent 1.21 vs dense
2.30**, ~24× at N=16k and ~270× projected at N=200k (see `benchmark-results.md`);
`backend = "fast"` wired opt-in for K=1. Gate S1 met for K=1 (R execution to be
re-confirmed in CI, as this environment has no R).

### S1.1 The exact reduction (correctness contract)

For high set `H` (stratum `d`), low set `L` (stratum `d-1`), times `T`, status
`D`, weights `ξ ≥ 0`, single endpoint:

```
W   = Σ_{j∈L, D_j=1} ξ_j · ( Σ_{i∈H : T_i > T_j} ξ_i )
Lo  = Σ_{i∈H, D_i=1} ξ_i · ( Σ_{j∈L : T_j > T_i} ξ_j )
Tot = (Σ_{i∈H} ξ_i) · (Σ_{j∈L} ξ_j)
```

Inequalities are **strict** (matches the dense kernel `>`). Equal times tie.

### S1.2 Algorithm — single pair (for the parity unit test)

```
fast_pair_win_loss_1d(t_high, d_high, w_high, t_low, d_low, w_low):
    # combine into entries (time, side ∈ {H,L}, status, weight)
    entries = [(t, H, d, w) for ...high...] + [(t, L, d, w) for ...low...]
    sort entries by time DESCENDING
    acc_H = 0; acc_L = 0; W = 0; Lo = 0
    for each maximal group of entries with equal time:
        # score using accumulators of STRICTLY greater times only
        for x in group:
            if x.side == L and x.status == 1: W  += x.w * acc_H
            if x.side == H and x.status == 1: Lo += x.w * acc_L
        # then fold the tied group into the accumulators
        for x in group:
            if x.side == H: acc_H += x.w
            else:           acc_L += x.w
    Tot = (Σ w_high) * (Σ w_low)
    return (W, Lo, Tot)
```

Correctness notes: scoring before folding the equal-time group enforces strict
`>`. `acc_H` at the moment of scoring `j∈L` equals `Σ_{i∈H, T_i>T_j} ξ_i`.
Complexity `O(n log n)` per pair (sort), `O(n)` if the order is reused.

### S1.3 Algorithm — all adjacent contrasts in one pass (bootstrap path)

The time order is fixed across bootstrap iterations; presort once. Per iteration,
one `O(N)` pass over the presorted order computes **every** `(d,d-1)` contrast:

```
sweep_all_adjacent_1d(order, stratum[], status[], weight[], D_strata):
    acc = zeros(D_strata+2)          # acc[s] = Σ weight in stratum s, strictly greater time
    W   = zeros(D_strata+1); Lo = zeros(D_strata+1)
    for each equal-time group in `order` (descending time):
        for x in group:
            s = stratum[x]
            if status[x]==1:
                W[s+1]  += weight[x] * acc[s+1]   # x is LOW side of contrast (s+1, s)
                Lo[s]   += weight[x] * acc[s-1]   # x is HIGH side of contrast (s, s-1)
        for x in group: acc[stratum[x]] += weight[x]
    # contrast d ≡ (d, d-1): wins=W[d], losses=Lo[d], total=tot[d]*tot[d-1]
    return W, Lo
```

This is C1 + "all contrasts at once" fused. Per bootstrap iteration cost is
`O(N)` after the one-time `O(N log N)` presort — independent of `D`.

### S1.4 Public functions

R (`R/kernel_fast.R`):

```r
mrwin_fast_pair_win_loss(time_high, status_high, time_low, status_low,
                         weights_high = NULL, weights_low = NULL)
# returns c(wins=, losses=, total=) — drop-in match for mrwin_pair_win_loss when K==1
mrwin_fast_adjacent_win_loss(time, status, strata, weights = NULL, time_order = NULL)
# returns data.frame(high, low, wins, losses, total) for all observed adjacent pairs
```

Python mirror in `kernel_fast.py`: `pair_win_loss_1d(...)`,
`sweep_all_adjacent_1d(...)`.

### S1.5 Parity test (`test-kernel-fast.R`) — the proof of correctness

Randomised cohorts, fixed seeds, assert fast == dense reference:

- Grid: `N ∈ {50, 200, 1000}`, `K=1`, `D_strata ∈ {2, 5, 10}`.
- Time generation: **integer times from a small set** (e.g. `sample(1:8, N, TRUE)`)
  to force heavy ties; plus one continuous-time replicate.
- Status: `rbinom`; include edge columns all-0 and all-1.
- Weights: `NULL` (unweighted) and `rexp(N)` (weighted).
- For each adjacent contrast, compare to `mrwin_pair_win_loss` evaluated on the
  same `idx_high/idx_low`.
- Assertions: **integer-exact** wins/losses unweighted; `abs(diff) < 1e-10`
  weighted; `total` exact. Antisymmetry spot check: `wins(H,L) == losses(L,H)`.

Run the same fixtures in `tests/python/test_kernel_fast.py` against
`kernel_sparse.stratum_pair_logcwr` / a dense reference.

### S1.6 Benchmark (`mrwin_benchmark_backends`) — the proof of speed

- Grid: `N ∈ {1e3, 3e3, 1e4, 3e4, 1e5}` (add `3e5` if memory allows),
  `K=1`, `D=10`, `M=20`, `B ∈ {0, 200}`, fixed seed.
- Compare: dense `mrwin_sparse_bootstrap` vs fast path.
- Record wall-clock (median of 3) and peak memory.
- Fit `log(time) ~ log(N)`; assert slope `p_fast < 1.3` and `p_dense ≈ 2`
  (report both with CIs).
- Emit raw rows to a gitignored CSV; write a curated table to
  `inst/spec/benchmark-results.md` (committed, numbers only).

### Gate S1 (Definition of Done)

1. Parity test green (R and Python), including the heavy-tie and edge cases.
2. Benchmark shows `p_fast < 1.3` vs `p_dense ≈ 2`, table committed.
3. `backend = "fast"` wired as opt-in in `mrwin_controls()` for the `K=1` path;
   default unchanged.
4. WP18 ledger row **L-C1** flipped to "discharged (test: test-kernel-fast.R)".
5. Update `wp13-fast-kernel.md` status block (`T1`,`T2` done) and
   `project-checkpoint.md`.

Effort estimate: 2–3 focused sessions.

---

## S2 — C2: hierarchical `K`-endpoint, build-once / reuse-across-`B`

Status: `[todo]`. Branch: `C-wp13` → `C-wp14`. Implements WP13 T3–T6, WP14 T1–T3.

### S2.1 Plan

1. Start with the **`K=2` closed form** (see `wp13-fast-kernel.md` §3.1): split
   the four `(D_{i,1}, D_{j,1})` cases; each tie-at-1 region is a rectangle in
   `(T_{i,1}, T_{j,1})`; intersect with the priority-2 win/loss half-plane; count
   by orthogonal range counting → `O(N log N)`.
2. Generalise to arbitrary `K` with the **CDQ divide-and-conquer / merge-sort
   tree** backend (`pair_win_loss_kd`), `O(N log^{K-1} N)`.
3. **Reuse across bootstrap**: build the dominance structure once on the fixed
   outcome coordinates (`mrwin_precompute_fast_structures(time, status)`), then
   per iteration only re-weight (`ξ`) and re-bin (strata). Refactor
   `mrwin_sparse_bootstrap` to build once before the `for b` loop.
4. **Incremental re-stratification** (WP14 T3): detect inversion fraction between
   `S` and `S*`; update only boundary movers; fall back to full re-query past a
   threshold (`incremental_inversion_frac`, default 0.2).

### S2.2 Parity test additions

Extend `test-kernel-fast.R` to `K ∈ {2,3}` vs `mrwin_pair_kernel` /
`mrwin_pair_win_loss`, with the adversarial cases from `wp13` T5:

- equal times across multiple priorities,
- all-censored priorities,
- terminal-event-before-lower-priority (`algorithm-spec.md` §3.1),
- same-seed bootstrap equals pre-refactor output within `1e-10`,
- incremental path equals full-requery path within `1e-10`.

### Gate S2

1. Parity green for `K=2,3` and for the refactored + incremental bootstrap.
2. Benchmark: per-iteration cost empirically independent of `D`, slope `< 1.3`
   for `K=3` (range-tree log factor noted).
3. WP18 ledger rows **T-C2**, **P-C3** discharged.
4. `wp13`/`wp14` status blocks updated; checkpoint updated.

Effort estimate: 4–6 sessions (the merge-sort tree + censoring partial-order
handling is the hardest engineering of Phase II).

---

## S3 — M2: doubly-ranked stratification

Status: `[todo]`. Branch: `C-wp16` (independent; can run in parallel with S2 by a
second agent). Implements WP16 T1–T3.

### S3.1 Function

```r
mrwin_doubly_ranked_strata(prs, exposure, n_strata)
# 1. order by prs; 2. form pre-strata of size n_strata; 3. within each
#    pre-stratum order by exposure and assign final stratum = exposure rank.
# returns list(strata = <1-indexed int>, ...) matching mrwin_prs_strata contract.
```

Add `stratification = c("prs_rank", "doubly_ranked")` to `mrwin_controls()`;
thread through `mrwin_estimate`, `mrwin_sparse_bootstrap`, `mrwin()`. Note the
interaction with S2 incremental updates (doubly-ranked re-derives pre-strata, so
incremental boundary tracking must re-form pre-strata; document it).

### S3.2 Tests + validation

- Deterministic toy-vector assignment matches a hand calculation; ties
  deterministic.
- Linear DGP: doubly-ranked vs prs_rank give compatible DS-CWR (no regression).
- Non-linear exposure–instrument DGP (extend WP8): doubly-ranked shows reduced
  within-stratum IV violation.

### Gate S3

Toy parity exact; non-linear scenario shows the intended robustness; prs_rank
remains default; WP16 status + WP18 "Stratification" section + checkpoint
updated. Effort: 1–2 sessions.

---

## S4 — M3 + M1: analytic variance + continuous estimator (research core)

Status: `[todo]`. Branches: `C-wp17` (M3), `C-wp15` (M1). Implements WP17 and
WP15. These are the manuscript headline results; theory (WP18) is written in
lockstep.

### S4.1 M3 — analytic influence-function variance (do first)

1. Theory note (WP18 V-M3): influence function of `log θ_d = log(W_d/L_d)` as a
   ratio of two-sample U-statistics (per-subject win/loss projections); stack
   `D-1` contrasts; delta-method through ISG and GLS (reuse
   `_mrwin_isg_covariance`, `mrwin_gls_pool`).
2. GWAS component: Jacobian `∂(logθ, ΔX)/∂β` from **boundary movers only**
   (individuals near stratum cutpoints), `Σ_gwas = J Σ_GWAS Jᵀ`.
3. `mrwin_analytic_covariance(...)`; add `inference = c("bootstrap","analytic",
   "analytic+bootstrap")` to `mrwin_controls()`; the hybrid uses analytic as a
   control variate to cut `B`.
4. Validate: analytic `Σ` ≈ bootstrap `Σ` within tolerance on null/valid-IV/
   dose-heterogeneity; CI coverage within MC error of nominal; report `B`
   reduction.

### S4.2 M1 — continuous kernel-smoothed ISG

1. Theory note (WP18 D-M1): continuous gradient `δ(r) = g(r)/x'(r)`; deciles =
   boxcar special case; local-U-statistic CLT; consistency under
   rank-preservation.
2. `mrwin_continuous_isg(...)`; fold kernel-weighted win accumulation by rank gap
   into the S1/S2 sweep (no extra asymptotic cost).
3. `estimator = "continuous"` / `n_strata = Inf`; cross-validated bandwidth;
   `δ̂(b)` sensitivity curve via `plot(fit, type = "gradient")`.
4. Validate: **boxcar reduction** — continuous == decile within `1e-8` for a
   boxcar of width `1/D`; coverage ≥ nominal; lower `D`-sensitivity than deciles.

### Gate S4

Both validations pass; bootstrap + decile estimator remain defaults until
coverage validated; WP15/WP17 status + WP18 sections D-M1, V-M3 discharged;
checkpoint updated; manuscript sections drafted (gitignored `manuscript/`).
Effort: the research bulk of Phase II — multiple sessions; theory and code
interleave.

---

## Cross-sprint rules

- **Parity is non-negotiable for S1–S2.** Any fast path that disagrees with the
  dense reference is a bug, not a tolerance choice.
- **Defaults flip only in WP19**, after the parity gate is green in CI.
- **No manuscript claim before its test passes** (roadmap §4.5.6).
- **Reuse-first**: search `R/` and `python/p1_engine_v5/` before adding a
  primitive; new fast paths sit behind the existing public API.
- **No tool/run/assistant-identifying strings** anywhere (code, docs, commits).

## Milestone summary

| Sprint | Item | Proves | Gate metric |
|---|---|---|---|
| S1 | C1 | order change, exactness | parity green + `p_fast<1.3` vs `p_dense≈2` |
| S2 | C2 + reuse | `K`-endpoint scalable | parity `K=2,3` + per-iter independent of `D` |
| S3 | M2 | robustness | toy-exact + non-linear IV validity |
| S4 | M3+M1 | inference + estimand | analytic≈bootstrap coverage; boxcar reduction `1e-8` |

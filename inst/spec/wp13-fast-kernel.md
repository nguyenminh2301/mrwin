# WP13 — Fast Hierarchical Win/Loss Kernel (C1 + C2)

Status block: `T1 [done] T2 [done] T3 [done: K=2] T4 [done: K=2] T5 [done: K<=2] T6 [done: K<=2 backend wired]  | K>=3 [open]`
Branch: `C-wp13` (from `C`).

Progress note (2026-06-14), R-verified under R 4.3.3:
- **S1 (K=1)** and **S2 (K=2)** landed. `python/p1_engine_v5/kernel_fast.py` +
  R `R/kernel_fast.R` (`mrwin_fast_pair_win_loss`, `mrwin_fast_adjacent_win_loss`,
  Fenwick-based 2D dominance counter for K=2).
- Parity: Python differential tests vs the brute-force oracle — 20k random K=2
  cohorts + 8k max-tie (support=2) cohorts, 0 mismatches; R testthat
  `test-kernel-fast.R` (K=1 and K=2) and `test-backend-fast.R` (end-to-end
  `mrwin(backend="fast")` for K=1/K=2 fast paths and K=3 fallback). Full suite 74
  groups, 0 failures.
- Scaling: K=1 fast tail exponent 1.21, K=2 exponent 1.18 (vs dense 2.30); see
  `benchmark-results.md`.
- `backend = "fast"` wired opt-in for K∈{1,2} via `.mrwin_pair_win_loss_backend`;
  K≥3 transparently falls back to the dense pair kernel (identical results).

### Complexity reality for K≥3 (honest note, supersedes the optimistic §3.1 claim)

The clean subquadratic fast paths are K=1 and K=2, both `Θ(N log N)`. For the
**generic** kernel at K≥3 the first-separation decomposition does not reduce to a
subproblem on a subset of *individuals* (the "tied at level 1" set is a relation
over *pairs*), and a direct orthogonal-range-counting formulation needs
dimension `2K`, so the log-power grows and the constant becomes unattractive by
K=3. A genuinely fast K≥3 path therefore needs to **exploit the nested
time-to-event structure** of the v5 endpoints (death terminal, censoring lower
priorities) rather than the generic per-priority `(t, status)` contract — this is
the open algorithmic problem and the next research checkpoint. Until then, K≥3
(the v5 flagship) uses the dense fallback and is correct but quadratic.
Depends on: WP3 (`mrwin_pair_win_loss`) and WP4 (`mrwin_estimate`) as the
correctness reference.
Blocks: WP14, WP15, WP19.

This work package replaces the `Θ(N²/D)` pair sweep with a subquadratic
algorithm that returns **identical** win/loss/total sums. It is the critical
path of Phase II.

---

## 1. Exact problem statement

Given two index sets `H` (high stratum `d`) and `L` (low stratum `d-1`), with
per-priority observed times `T ∈ ℝ^{N×K}`, event indicators `D ∈ {0,1}^{N×K}`,
and non-negative weights `ξ ∈ ℝ^N`, compute:

```
W = Σ_{i∈H} Σ_{j∈L} ξ_i ξ_j · 1{ h(i,j) = +1 }
Lo = Σ_{i∈H} Σ_{j∈L} ξ_i ξ_j · 1{ h(i,j) = -1 }
Tot = (Σ_{i∈H} ξ_i)(Σ_{j∈L} ξ_j)
```

where `h` is the hierarchical kernel defined in `algorithm-spec.md` §4.1 and
implemented in `R/kernel.R::mrwin_pair_kernel`. The reference semantics, per
priority `k` (first informative priority wins):

```
i wins  at k  ⇔  D_{j,k}=1 ∧ T_{i,k} > T_{j,k}     (and pair not yet decided)
i loses at k  ⇔  D_{i,k}=1 ∧ T_{j,k} > T_{i,k}     (and pair not yet decided)
otherwise the pair ties at k and proceeds to k+1.
```

The output must match `mrwin_pair_win_loss(...)` exactly (integer counts when
unweighted; `≤1e-10` relative for weighted sums).

---

## 2. C1 — Single endpoint (`K = 1`)

### 2.1 The reduction

Drop the priority index. The win sum factorises:

```
W = Σ_{j∈L, D_j=1} ξ_j · ( Σ_{i∈H : T_i > T_j} ξ_i )
Lo = Σ_{i∈H, D_i=1} ξ_i · ( Σ_{j∈L : T_j > T_i} ξ_j )
```

Each inner sum is a **suffix-weighted sum over a fixed time order**. Therefore:

1. Sort the union `H ∪ L` by `T` **once** (this order is reusable across all
   bootstrap iterations because `T` is fixed — exploited fully in WP14).
2. Sweep in decreasing `T`, maintaining a running weighted sum of `H`-members
   seen so far (`acc_H`). When the sweep reaches a member `j ∈ L` with `D_j=1`,
   its contribution to `W` is `ξ_j · acc_H`.
3. Symmetrically maintain `acc_L` for the loss sum.

Ties in `T` (equal times) must be handled so that `T_i > T_j` is strict: process
all equal-time entries as a group and add them to the accumulator **after**
scoring the events at that exact time. Validate against the dense kernel on
exact-tie inputs (the dense kernel uses strict `>`).

Complexity: one sort `Θ(N log N)` (amortised away in WP14) + one linear sweep
`Θ(N)`. Per stratum pair: `Θ((|H|+|L|))` after the global sort.

### 2.2 Tasks

- **T1** — Implement `mrwin_fast_pair_win_loss_1d(time, status, idx_high,
  idx_low, weights)` in a new file `R/kernel_fast.R`. Pure R first
  (vectorised cumulative sums over the sorted order). Mirror the Python
  reference in `python/p1_engine_v5/kernel_fast.py::pair_win_loss_1d`.
- **T2** — Parity test `tests/testthat/test-kernel-fast.R`: for random cohorts
  (`N` up to 2000, `K=1`, including heavy ties and all-censored columns), assert
  exact equality with `mrwin_pair_win_loss` for unweighted and weighted cases.

---

## 3. C2 — Hierarchical (`K ≥ 2`)

### 3.1 The decomposition

Win total = wins decided at priority 1 + (tied at 1) ∧ (win at 2) + … The first
term is a C1 problem. The coupling term is the hard part:

```
W = Σ_k  Σ_{(i,j)}  1{ tied at 1..k-1 } · 1{ i wins at k }.
```

For fixed `k`, the event `{tied at 1..k-1} ∧ {i wins at k}` is a conjunction of
a bounded number (`≤ 2(k-1)+2`) of threshold conditions on the coordinates
`(T_{i,1},T_{j,1},…,T_{i,k},T_{j,k})`, branched over the relevant event
indicators `D_{·,·}`. Counting weighted pairs satisfying a conjunction of
coordinate-threshold conditions is **orthogonal range counting**, solvable with
a `(k-1)`-level merge-sort tree / range tree.

Worked case `K = 2`:

```
tie at 1  ⇔  ¬(i wins 1) ∧ ¬(i loses 1)
          = ( D_{j,1}=0 ∨ T_{i,1} ≤ T_{j,1} ) ∧ ( D_{i,1}=0 ∨ T_{j,1} ≤ T_{i,1} ).
win at 2  ⇔  D_{j,2}=1 ∧ T_{i,2} > T_{j,2}.
```

Split the four `(D_{i,1},D_{j,1})` cases; in each the tie-at-1 region is a fixed
rectangle/half-plane in `(T_{i,1},T_{j,1})`. Combined with the `T_{i,2}>T_{j,2}`
condition this is a constant-dimension range count → `Θ(N log N)` per case →
`Θ(N log N)` total for `K=2`, `Θ(N log² N)` for `K=3`.

### 3.2 Implementation strategy

Two acceptable implementations; pick per the cost note in WP19:

- **(a) Offline divide-and-conquer (CDQ) over priorities.** Recursively split,
  sort by one coordinate, sweep with a Fenwick tree keyed by the next
  coordinate. Cleanest for general `K`; `Θ(N log^{K-1} N)`.
- **(b) Sequential "still-tied set" sweep.** Process priorities in order; carry
  the set of still-tied pairs implicitly as a collection of canonical rectangles
  in time-rank space; intersect with the next priority's win/loss half-planes.
  More intuitive for `K=3`; matches the existing `decided`-mask semantics.

Both must reproduce the dense kernel exactly. Implement (a) as the general
backend; (b) optionally as a `K=3` fast path if it benchmarks faster.

### 3.3 Tasks

- **T3** — Implement the general dominance-counting backend
  `pair_win_loss_kd(...)` in `python/p1_engine_v5/kernel_fast.py` for arbitrary
  `K`, with the `K=2` closed-form path as the first correctness milestone.
- **T4** — Port to R: `mrwin_fast_pair_win_loss(time, status, idx_high,
  idx_low, weights)` dispatching to the 1d path when `K=1`.
- **T5** — Parity test extending `test-kernel-fast.R` to `K=2,3` against
  `mrwin_pair_kernel` / `mrwin_pair_win_loss`, including:
  - adversarial ties (equal times across multiple priorities),
  - all-censored priorities,
  - terminal-event-before-lower-priority cases (per `algorithm-spec.md` §3.1),
  - antisymmetry spot check (`W(H,L) == Lo(L,H)`).
- **T6** — Wire a new opt-in backend `backend = "fast"` in `mrwin_controls()`
  and route `mrwin_sparse_estimate` / `mrwin_sparse_bootstrap` to the fast pair
  function when selected. Default backend unchanged until WP19 parity is green.

---

## 4. Acceptance criteria

1. `mrwin_fast_pair_win_loss` returns integer-exact `wins`/`losses` and
   `≤1e-10` weighted sums vs the dense reference across all parity cases.
2. Antisymmetry and zero-diagonal properties preserved.
3. Scaling benchmark (WP19 harness) shows empirical near-`N log N` growth for
   `K=1` and near-`N log^{K-1} N` for `K=2,3`, across `N ∈ {1e3, 1e4, 1e5}`,
   versus the quadratic reference.
4. No change to any existing default output (the fast path is opt-in).

## 5. Definition of Done

All of §4 plus: status block all `[done]`; `wp18-theory-foundation.md` Lemma C1
and Theorem C2 marked "implemented + verified"; `project-checkpoint.md` updated.

## 6. References (rank-counting / win-ratio computation)

- Fenwick, P. M. (1994). A new data structure for cumulative frequency tables.
  *Software: Practice and Experience*, 24(3), 327–336.
- Pocock, S. J., et al. (2012). The win ratio. *Eur Heart J*, 33(14),
  1744–1749. (estimand definition)
- Bebu, I., & Lachin, J. M. (2016). Large-sample inference for a win ratio.
  *Biostatistics*, 17(1), 178–191. (variance baseline)
- Standard computational-geometry range-counting (merge-sort tree / CDQ
  divide-and-conquer); confirm exact citation in WP18 before publishing.

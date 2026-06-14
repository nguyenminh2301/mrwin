# WP14 — Bootstrap Acceleration via Fixed-Kernel Reuse (C3)

Status block: `T1 [todo] T2 [todo] T3 [todo] T4 [todo]`
Branch: `C-wp14` (from `C`, after WP13 merges).
Depends on: WP13 (fast pair structures).
Blocks: WP19.

Phase II's largest absolute saving: the multiplier bootstrap currently redoes
the entire win/loss sweep every iteration because re-stratification changes
group membership. This WP makes each iteration cheap by **reusing structures
built once on the fixed outcome coordinates** and by **updating across the small
perturbation `β* ≈ β`** instead of rebuilding.

---

## 1. The two levers

### 1.1 Lever A — Reuse fixed-order structures across `B`

The time order (C1) and the dominance structures (C2) depend only on `(time,
status)`, which are invariant across bootstrap iterations. Build them **once**,
before the bootstrap loop. Per iteration, only:

- recompute `S* = G β*` and the stratum labels (rank-bands),
- recompute weights `ξ_i` (and IPTW `ω_i` if adjustment is on),
- re-query the prebuilt structures with the new labels/weights.

This alone turns `B · Θ(N²/D)` into `B · Θ(N log^{K-1} N)` and removes per-iter
structure construction.

Implementation: refactor `mrwin_sparse_bootstrap` so that the fixed-order index
and any merge-sort-tree nodes are computed before the `for (b in 1:B)` loop and
passed into the per-iteration win/loss call. The per-iteration call takes
`(labels, weights)` and returns all `D-1` adjacent contrasts in one sweep
(tag each individual by stratum; a single pass over the fixed order accumulates
per-stratum prefix sums).

### 1.2 Lever B — Incremental re-stratification under `β* ≈ β`

`S* = G β*` with `β* = β + δβ`, `δβ ~ N(0, Σ_GWAS)`. When `Σ_GWAS` is small and
SNPs are LD-clumped, the rank order of `S*` differs from that of `S` by a small
number of adjacent inversions `I`. Maintain the stratum boundaries incrementally:

- Precompute the baseline rank order of `S`.
- Per iteration, compute `S*`, then update the order by counting/applying only
  the `I` inversions (e.g. via an adaptive sort that is `O(N + I log N)` on
  nearly-sorted input), and move only the individuals that cross a stratum
  boundary.
- Recompute win/loss **incrementally** for the contrasts whose membership
  changed, not all of them.

Cost per iteration becomes `O(I log N + (boundary movers) · log N)`, which is
sublinear when perturbations are small.

Robustness gate: compute `I` (or a cheap proxy: Spearman footrule between `S`
and `S*`). If it exceeds a configurable threshold (default: order changes for
> 20% of individuals), **fall back** to the Lever-A full re-query for that
iteration. Record fallback rate as a diagnostic.

---

## 2. Tasks

- **T1** — Refactor `mrwin_sparse_bootstrap` (and the Python
  `multiplier_bootstrap`) to build fixed-order structures once and pass them
  into the per-iteration win/loss. Parity: identical bootstrap output (same
  seed) to the current implementation within `1e-10`.
- **T2** — Implement the single-pass "all adjacent contrasts at once" sweep so
  one ordered pass yields every `(d, d-1)` win/loss. Replaces the per-contrast
  loop.
- **T3** — Implement Lever B incremental re-stratification with the inversion
  threshold fallback. Add `bootstrap_incremental = TRUE/FALSE` and
  `incremental_inversion_frac` to `mrwin_controls()`. Parity vs Lever A.
- **T4** — Diagnostics: record per-iteration inversion fraction and fallback
  rate in the bootstrap result object; expose in `summary()`.

---

## 3. Acceptance criteria

1. Same-seed bootstrap output equals the pre-WP14 result within `1e-10`
   (parity is non-negotiable; this is a speed change only).
2. Per-iteration cost is empirically independent of `D` and sub-quadratic in `N`
   (WP19 benchmark across `N ∈ {1e3,1e4,1e5}`, `B=200`).
3. Incremental path matches full-requery path within `1e-10` and degrades
   gracefully (correct results) when forced past the inversion threshold.

## 4. Definition of Done

§3 satisfied; status block all `[done]`; WP18 Proposition C3 marked
"implemented + verified"; `project-checkpoint.md` updated.

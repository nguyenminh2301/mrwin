# Phase II — Direction Evaluation & Remaining Roadmap

Created: 2026-06-18. Audience: maintainer + contributor agents.
This is the strategic decision document: what is left, which direction is most
valuable, the step-by-step path of each to its final result, and the
recommended critical path. It re-prioritises the WP14–WP19 specs given what has
actually landed.

---

## 1. Where the project stands

The original motivation was two-fold: (a) break the `O(N²)` win/loss bottleneck,
and (b) build a *publishable* contribution that goes beyond a plain fast
win-ratio calculator (e.g. jaFS).

**Goal (a) is solved.** The fast hierarchical kernel (K=1,2,3) + the Rcpp port
make `mrwin(backend="fast")` biobank-scale: K=3 is 72–300× faster than the dense
backend, K=3 at N=80k runs in 0.78 s, end-to-end `mrwin()` is 16× faster at
N=2k/B=100 with identical results, and the empirical exponent is ≈1.1–1.13. The
remaining performance cost is the `B`-iteration bootstrap loop (each iteration
re-stratifies + runs the now-cheap kernel).

**Goal (b) is the remaining value.** With speed solved, the marginal value is now
*methodological* — the pieces that make this a citable, differentiated method
rather than "a fast stratified win ratio".

Landed: WP13 (fast kernel + Rcpp), M2 core (`mrwin_doubly_ranked_strata`).
Remaining: M1, M3, M2-wiring, WP14, WP19, WP18, K≥4.

---

## 2. Direction scorecard

Value = contribution to the end goal (publishable, validated, biobank-scale
method). Effort in focused steps. Risk = chance of a hard/uncertain outcome.

| Dir | What | Value | Effort | Risk | Note |
|---|---|---|---|---|---|
| **M1** | Continuous kernel-smoothed ISG (removes `D`) | ★★★★★ | 6 | med-high | the *new method*; headline novelty; differentiates from all prior work |
| **M3** | Analytic influence-function variance | ★★★★☆ | 5 | med | rigorous inference **and** removes the `B`-loop (the remaining perf cost) |
| **WP19** | Simulation validation + release | ★★★★☆ | 6 | low-med | without coverage/type-I evidence the method is not credible/publishable |
| **WP18** | Theory manuscript | ★★★★☆ | 4* | low | turns the work into a citable artifact; interleaves with M1/M3/WP19 |
| **M2-wire** | Thread doubly-ranked through pipeline | ★★★☆☆ | 2 | low | quick robustness completion; core already done |
| **WP14** | Bootstrap build-once/reuse + incremental | ★★☆☆☆ | 4 | med | constant-factor only now; partly moot if M3 removes the `B`-loop |
| **K≥4** | Extend fast kernel beyond 3 priorities | ★☆☆☆☆ | 4 | med | v5 endpoint is K=3; low demand |

\*WP18 is continuous, not a single block.

---

## 3. The most valuable direction

**M1 (continuous kernel-smoothed ISG) is the single most valuable direction.**
Reasoning: the speed problem is solved, so value now lives in the science, and M1
is the genuine *new estimator* — it replaces the arbitrary decile count `D` with
a bandwidth-controlled local gradient for which deciles are the boxcar special
case, with cleaner U-process asymptotics. It is exactly the "new method" the
project set out to find, it closes the package's own flagged `D`-sensitivity gap,
and it is what distinguishes this from existing win-ratio and non-linear-MR work.

**M3 (analytic variance) is the highest-leverage companion** because it is both a
rigorous inference contribution *and* the fix for the remaining performance cost
(it removes the `B`-iteration bootstrap), and its influence-function machinery is
reused by M1's inference. For these reasons the recommended *execution order*
puts M3 first even though M1 is the headline (see §5).

---

## 4. Roadmap of each direction to its final result

Each step is an independently testable unit (commit + tests), following the
same discipline used for WP13 (differential testing against a reference, R
verification, honest benchmarking).

### M1 — continuous ISG (final result: a `estimator="continuous"` mode with validated coverage and a manuscript section). 6 steps.
1. Theory note: define `δ(r) = g(r)/x'(r)`; local-U-statistic estimator;
   consistency under rank-preservation; **boxcar limit = decile estimator**;
   CLT sketch. (→ WP18)
2. Prototype `mrwin_continuous_isg()` on the dense path; **verify** it reproduces
   the decile DS-CWR exactly when the kernel is a boxcar of width `1/D`.
3. Fold kernel-weighted win accumulation by rank-gap into the fast sweep
   (R + C++), so the continuous estimator inherits the subquadratic speed.
4. Cross-validated bandwidth + `δ̂(b)` sensitivity curve + `plot(type="gradient")`.
5. Coverage simulation on the WP8 DGP; show lower `D`-sensitivity than deciles.
6. Wire `estimator = c("decile","continuous")` into `mrwin_controls()`; docs.

### M3 — analytic variance (final result: `inference="analytic"` matching the bootstrap). 5 steps. **[steps 1–4 DONE 2026-06-18; only IPTW + manuscript remain]**
Done: influence-function `Σ_sampling` (closed form) + exact GWAS-only resample
`Σ_gwas`, validated vs the full bootstrap (`se` ratio 0.997–1.009), wired as
`mrwin(inference="analytic")` (7× at σ_β=0). Remaining: IPTW influence terms
(`adjustment="ordinal_iptw"`) and the manuscript section.
1. Theory note: influence function of `log θ_d` (ratio of two-sample U-stats →
   per-subject win/loss projections); Jacobian of `(logθ,ΔX)` through the stratum
   cutpoints for GWAS-uncertainty propagation. (→ WP18)
2. Implement `mrwin_analytic_covariance()` (R + Python ref): `Σ_sampling` from the
   IF, `Σ_gwas = J Σ_GWAS Jᵀ` from boundary movers; reuse `_mrwin_isg_covariance`,
   `mrwin_gls_pool`, `.mrwin_fieller_ci`.
3. Wire `inference = c("bootstrap","analytic","analytic+bootstrap")`.
4. **Validate**: analytic `Σ` ≈ bootstrap `Σ` on null/valid-IV/dose-heterogeneity;
   CI coverage within MC error; quantify the speedup from dropping the `B`-loop.
5. Docs + manuscript section; keep bootstrap as default until coverage passes.

### WP19 — validation + release (final result: defaults flipped, CI green, release-ready). 6 steps.
1. Parity-gate harness `mrwin_verify_fast_dense_parity()` + CI wiring.
2. Scaling benchmark in `R/benchmark.R` + committed `benchmark-results.md` table;
   assert exponent bounds.
3. Statistical-validation grids: type-I (null), power (valid IV), pleiotropy
   (SDPD vs γ), weak-IV (Fieller), discordant-components.
4. Continuous-ISG coverage (M1) + analytic-variance coverage (M3) grids.
5. Flip defaults: `backend="fast"` default once the parity gate is green in CI;
   keep `dense` for debugging.
6. `R CMD check --as-cran`, coverage ≥ 80%, release checklist.

### M2-wiring — 2 steps.
1. Thread `stratification=c("prs_rank","doubly_ranked")` through
   `mrwin_controls` / `mrwin_estimate` / `mrwin_sparse_bootstrap` / `mrwin()`
   (exposure available at every re-stratification).
2. Validate on a non-linear exposure–instrument DGP (IV validity within strata).

### WP14 — bootstrap build-once/reuse — 4 steps (deferred; do only if M3 does not remove the `B`-loop for a given workflow).
1. Precompute fixed time-order / outcome structures once before the `B`-loop.
2. Single-pass all-adjacent-contrasts sweep (fuse contrasts).
3. Incremental re-stratification with inversion-threshold fallback.
4. Parity (same-seed == pre-WP14 within 1e-10) + benchmark.

### K≥4 — 4 steps (low priority).
1. General `d`-dim weighted dominance counter (CDQ recursion to a base 2D).
2. General-K tie-split assembly over `4^{K-1}` regimes.
3. Differential test vs `dense_pair_win_loss_kd` for K=4,5.
4. C++ port + dispatch.

### WP18 — manuscript — continuous (~4 effective steps): one section per landed
result, discharge each proof-obligation ledger row when its test passes, final
pre-submission consistency pass.

---

## 5. Recommended critical path (to a publishable + released package)

```
M3 (5)  ──►  M1 (6)  ──►  WP19 validation (6)  ──►  release
   │            │                 ▲
   └── M2-wire (2, parallel) ─────┘
WP18 manuscript (4): interleaved, one section per landed result.
```

Order rationale: **M3 first** — its influence-function machinery is the
foundation M1's inference reuses, it removes the remaining `B`-loop cost (so the
whole pipeline gets fast end-to-end, not just the kernel), and it is lower-risk
than M1. **M1 second** — the headline estimator, built on the IF machinery.
**M2-wiring** runs in parallel (independent, 2 steps). **WP19** validates M1+M3
and flips defaults. **WP18** is written as each result lands.

**Remaining steps to publishable + released:** **~23** (M3 5 + M1 6 + M2-wire 2 +
WP19 6 + WP18 4). Optional/deferred: WP14 (4) + K≥4 (4) = +8.

Deferred deliberately: **WP14** (constant-factor; largely superseded by M3
removing the `B`-loop) and **K≥4** (v5 is K=3). Revisit WP14 only for workflows
where analytic variance is not applicable and very large `B` is required.

---

## 6. One-line recommendation

The speed goal is done; invest the remaining effort in the **methodological
track — M3 then M1, validated by WP19 and written up in WP18** — because that is
what converts a fast win-ratio engine into a citable, differentiated method. Do
**M3 next**: highest combined value (rigorous inference + removes the last
performance bottleneck + foundation for M1).

# Benchmark Results

Curated scaling numbers for Phase II. Raw artifacts are gitignored; reproduce
with the harnesses named below. Numbers are wall-clock medians on a single core;
absolute values are environment-dependent, the **empirical exponent** is the
portable signal.

## S1 — single-endpoint (K=1) win/loss: dense O(N²/D) vs fast O(N log N)

Harness: `python -m p1_engine_v5.benchmark_fast_kernel` (D=10, ties present).
`dense` = numpy boolean blocks (current backend style); `fast` =
`kernel_fast.sweep_all_adjacent_1d` (pure Python). Comparing pure-Python `fast`
against numpy-vectorised `dense` is conservative — a compiled `fast` widens the
gap further.

| N | dense (s) | fast (s) | speedup |
|---:|---:|---:|---:|
| 1,000 | 0.0008 | 0.0003 | ~3× |
| 4,000 | 0.0093 | 0.0011 | ~8× |
| 8,000 | 0.0353 | — | |
| 16,000 | 0.1377 | 0.0057 | ~24× |
| 32,000 | 0.8608 | — | |
| 64,000 | — | 0.0330 | |
| 200,000 | ~33 (∝N² extrapolation) | 0.1208 | ~270× |

Empirical log-log exponent:

- dense: full-fit `p = 1.85`, large-N **tail `p = 2.30`** (super-quadratic; the
  full fit is dragged down by numpy block overhead at small N).
- fast: full-fit `p = 1.17`, tail `p = 1.21`.

Gate (`implementation-plan.md` S1.6): `p_fast < 1.3` ✅ and `p_dense_tail > 1.6`
✅. This is the concrete evidence that the algorithm changes **order**, not just
the constant factor.

Parity: `tests/python/test_kernel_fast.py` (400 random pairs + edge columns +
antisymmetry + 120 multi-stratum sweeps) passes bit-for-bit (integer-exact
unweighted, `<1e-10` weighted). R parity verified under **R 4.3.3**:
`tests/testthat/test-kernel-fast.R` and `tests/testthat/test-backend-fast.R`
(end-to-end `mrwin(backend="fast")`) pass; full suite 71 groups, 0 failures.

## S2 — hierarchical K=2 win/loss fast path

`fast_pair_win_loss_2d` (level-1 sweep + four-regime tie-split with a Fenwick 2D
dominance counter), one balanced pair, half high / half low:

| N (per side ×2) | fast_2d (s) |
|---:|---:|
| 2,000 | 0.011 |
| 8,000 | 0.055 |
| 32,000 | 0.263 |
| 128,000 | 1.51 |

Empirical exponent `p = 1.18` (gate `< 1.3` ✅). Parity: 20k random K=2 cohorts
+ 8k max-tie (support=2) cohorts vs the brute-force oracle, **0 mismatches**;
R `test-kernel-fast.R` (K=2) and `test-backend-fast.R` (K=2 end-to-end) green.

K≥3: no fast path yet (dense fallback). See `wp13-fast-kernel.md` honesty note.

_Last updated: 2026-06-14 (S1 + S2/K=2 landed, R-verified)._

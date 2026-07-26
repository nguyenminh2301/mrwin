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

Gate (`dev/implementation-plan.md` S1.6): `p_fast < 1.3` ✅ and `p_dense_tail > 1.6`
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

## K=3 hierarchical win/loss fast path (the v5 flagship)

`fast_pair_win_loss_3d` = K=2 fast on levels 1–2 + a level-3 term over 16
regimes using a **CDQ divide-and-conquer 3D dominance counter** (O(N) memory):

| N (per side ×2) | fast_3d (s) |
|---:|---:|
| 2,000 | 0.068 |
| 8,000 | 0.334 |
| 32,000 | 1.71 |
| 128,000 | 9.17 |

Empirical exponent `p = 1.18` (≈ N log² N). An earlier dense-2D-BIT version was
O(N²) memory (exponent 1.58, 18 s at N=64k); the CDQ rewrite fixed it.

Parity: K=3 fast vs dense oracle — 35k+ random cohorts incl. max-tie
(support=2), 0 mismatches; the 3D counter vs a dense-2D-BIT oracle — 160k
comparisons over all comparison-direction combos, 0 mismatches. R
`test-kernel-fast.R` (K=3) + `test-backend-fast.R` (K=3 end-to-end) green.

K≥4: dense fallback (correct).

### Rcpp port — the kernel now realises the order advantage in wall-clock

`mrwin_fast_pair_win_loss` dispatches to the compiled `mrwin_fast_pair_cpp`
(`src/fast_kernel.cpp`), an exact port of the validated R/Python kernel. The
earlier pure-R CDQ recursion was correct but had a large interpreted constant
(K=3 only overtook dense around n≈2000). With the C++ kernel the order advantage
is fully realised:

| case | fast(C++) vs dense | match |
|---|---|---|
| K=3, n=2000 | **72× faster** | exact |
| K=3, n=5000 | **300× faster** (dense ≈ 10 s) | exact |

C++ fast-path scaling (one balanced pair, fast-only; dense is infeasible here):

| N (per side ×2) | K=3 (s) | K=1 (s) |
|---:|---:|---:|
| 5,000 | 0.034 | — |
| 20,000 | 0.163 | 0.006 |
| 80,000 | 0.78 | — |
| 200,000 | — | 0.083 |

K=3 at N=80,000 in 0.78 s vs a dense `O(N²)` extrapolation of ~40 min (and OOM)
— biobank-scale feasibility. Empirical exponent ≈ 1.1 (K=1) / ≈ 1.13 (K=3).
Validation: the C++ kernel matches the dense oracle over 16k random K=1/2/3
cohorts (0 mismatches) and equals the pure-R fast path (triangle test
`test-kernel-fast.R`). A differential test caught one C++ bug (double-`eq` regime
collapsing the level-2 equality), now fixed. Full R 4.3.3 suite: 82 groups, 0
failures.

_Last updated: 2026-06-18 (Rcpp port landed, R-verified)._

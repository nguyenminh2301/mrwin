"""
benchmark_fast_kernel.py — WP13 / Phase II S1 scaling proof.

Times the current-style dense adjacent win/loss (O(N^2 / D), numpy boolean
blocks, mirroring kernel_sparse.py) against the fast sorted sweep
(O(N log N), kernel_fast.sweep_all_adjacent_1d) and fits the empirical
log-log exponent. The gate (implementation-plan.md S1.6) is:
    p_fast < 1.3   and   p_dense ~ 2.

Run: python -m p1_engine_v5.benchmark_fast_kernel
Raw CSV output is gitignored; the curated table lives in
inst/spec/benchmark-results.md.
"""
from __future__ import annotations
import time as _time
import numpy as np

from .kernel_fast import sweep_all_adjacent_1d


def _make(n, n_strata, seed):
    rng = np.random.default_rng(seed)
    t = rng.integers(1, max(2, n // 4), size=n).astype(float)  # ties present
    s = rng.integers(0, 2, size=n)
    strat = rng.integers(1, n_strata + 1, size=n)
    w = rng.exponential(1.0, size=n)
    return t, s, strat, w


def dense_adjacent_np(t, s, strat, w, n_strata):
    """O(N^2 / D) adjacent win/loss using numpy boolean blocks (backend style)."""
    out = {}
    for d in range(2, n_strata + 1):
        hi = np.where(strat == d)[0]
        lo = np.where(strat == d - 1)[0]
        if hi.size == 0 or lo.size == 0:
            continue
        th = t[hi][:, None]; tl = t[lo][None, :]
        sh = s[hi][:, None]; sl = s[lo][None, :]
        wh = w[hi][:, None]; wl = w[lo][None, :]
        win = (sl == 1) & (th > tl)
        loss = (sh == 1) & (tl > th)
        W = float(np.sum(wh * wl * win))
        Lo = float(np.sum(wh * wl * loss))
        out[d] = (W, Lo, float(w[hi].sum() * w[lo].sum()))
    return out


def _median_time(fn, repeats=3):
    times = []
    for _ in range(repeats):
        t0 = _time.perf_counter()
        fn()
        times.append(_time.perf_counter() - t0)
    return float(np.median(times))


def _slope(ns, ts):
    return float(np.polyfit(np.log(np.array(ns)), np.log(np.array(ts)), 1)[0])


def _tail_slope(ns, ts, k=3):
    return _slope(ns[-k:], ts[-k:])


def main(n_strata=10, seed=20260614):
    dense_ns = [500, 1000, 2000, 4000, 8000, 16000, 32000]
    fast_ns = [1000, 4000, 16000, 64000, 200000]

    print(f"{'method':6s} {'N':>8s} {'time_s':>12s}")
    dense_t = []
    for n in dense_ns:
        t, s, strat, w = _make(n, n_strata, seed)
        tt = _median_time(lambda: dense_adjacent_np(t, s, strat, w, n_strata))
        dense_t.append(tt)
        print(f"{'dense':6s} {n:>8d} {tt:>12.5f}")
    fast_t = []
    for n in fast_ns:
        t, s, strat, w = _make(n, n_strata, seed)
        tl = t.tolist(); sl = s.tolist(); strl = strat.tolist(); wl = w.tolist()
        order = sorted(range(n), key=lambda x: tl[x], reverse=True)
        tt = _median_time(
            lambda: sweep_all_adjacent_1d(tl, sl, strl, wl, n_strata, order)
        )
        fast_t.append(tt)
        print(f"{'fast':6s} {n:>8d} {tt:>12.5f}")

    p_dense = _slope(dense_ns, dense_t)
    p_fast = _slope(fast_ns, fast_t)
    p_dense_tail = _tail_slope(dense_ns, dense_t)
    p_fast_tail = _tail_slope(fast_ns, fast_t)
    print()
    print(f"empirical exponent  dense p = {p_dense:.3f}  (tail {p_dense_tail:.3f}, expect ~2)")
    print(f"empirical exponent  fast  p = {p_fast:.3f}  (tail {p_fast_tail:.3f}, gate: < 1.3)")
    print(f"GATE p_fast<1.3 : {'PASS' if p_fast < 1.3 else 'FAIL'}")
    print(f"GATE p_dense~2  : {'PASS' if p_dense_tail > 1.6 else 'FAIL'}")
    return {"p_dense": p_dense, "p_fast": p_fast,
            "p_dense_tail": p_dense_tail, "p_fast_tail": p_fast_tail,
            "dense": list(zip(dense_ns, dense_t)),
            "fast": list(zip(fast_ns, fast_t))}


if __name__ == "__main__":
    main()

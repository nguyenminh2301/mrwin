"""
fast_kernel_benchmark.py — Reproducible scaling benchmark for the
sub-quadratic hierarchical win/loss algorithm.

Measures wall-clock runtime of ``fast_pair_win_loss`` (O(N log^{K-1} N)) against
the dense O(n_i n_j) kernel, across two-group sizes spanning into the 100k+
range, and fits the empirical growth exponent. Produces a CSV and a printed
table suitable for inclusion in a manuscript / preprint.

Reproducibility:
  - Fixed RNG seed derived from each size; same machine class -> same numbers.
  - Exactness is re-verified at every size where the brute path is run.

Usage:
  python -m p1_engine_v5.fast_kernel_benchmark
  python -m p1_engine_v5.fast_kernel_benchmark --max-fast 262144 --reps 3 \
        --out results/fast_kernel_scaling.csv

Notes:
  - This times the PURE-PYTHON reference. A compiled (Rust/extendr) port shares
    the asymptotics with a far smaller constant; treat absolute seconds as an
    upper bound and the log-log slope as the portable result.
  - Cost is per stratum pair per bootstrap iteration. The full pipeline runs
    this (D-1) * B times.
"""
from __future__ import annotations

import argparse
import csv
import math
import time
from pathlib import Path

import numpy as np

from .fast_kernel import fast_pair_win_loss
from .kernel_sparse import stratum_pair_kernel


def _make_groups(n_high, n_low, K, seed):
    """Continuous times + right-censoring + exponential multiplier weights."""
    rng = np.random.default_rng(seed)
    TH = rng.exponential(1.0, size=(n_high, K))
    DH = rng.integers(0, 2, size=(n_high, K))
    wH = rng.exponential(1.0, size=n_high)
    TL = rng.exponential(1.0, size=(n_low, K))
    DL = rng.integers(0, 2, size=(n_low, K))
    wL = rng.exponential(1.0, size=n_low)
    return TH, DH, wH, TL, DL, wL


def _brute(TH, DH, wH, TL, DL, wL):
    H = stratum_pair_kernel(TH, DH, TL, DL)
    W = np.outer(wH, wL)
    return float(np.sum(W * (H == 1))), float(np.sum(W * (H == -1)))


def _time(fn, reps):
    best = math.inf
    out = None
    for _ in range(reps):
        t0 = time.perf_counter()
        out = fn()
        best = min(best, time.perf_counter() - t0)
    return best, out


def run(sizes, K, max_brute, reps, seed0):
    rows = []
    prev_n = prev_t = None
    print(f"Hierarchical win/loss scaling  (K={K}, weighted, censored, "
          f"best-of-{reps})")
    print(f"{'N/group':>9} {'fast(s)':>10} {'brute(s)':>10} {'speedup':>8} "
          f"{'slope':>6} {'exact':>6}")
    for i, n in enumerate(sizes):
        TH, DH, wH, TL, DL, wL = _make_groups(n, n, K, seed0 + i)
        ft, (fw, fl, _) = _time(
            lambda: fast_pair_win_loss(TH, DH, TL, DL, wH, wL), reps)

        bt = None
        exact = ""
        if n <= max_brute:
            bt, (bw, bl) = _time(lambda: _brute(TH, DH, wH, TL, DL, wL), 1)
            err = max(abs(fw - bw), abs(fl - bl))
            exact = "ok" if err < 1e-6 else f"ERR{err:.1e}"

        slope = ""
        if prev_t is not None:
            slope = f"{math.log(ft / prev_t) / math.log(n / prev_n):.2f}"
        prev_n, prev_t = n, ft

        sp = f"{bt / ft:.2f}" if bt else "--"
        bts = f"{bt:10.3f}" if bt is not None else f"{'(skip)':>10}"
        print(f"{n:>9} {ft:10.3f} {bts} {sp:>8} {slope:>6} {exact:>6}")
        rows.append({
            "n_per_group": n, "K": K,
            "fast_seconds": ft,
            "brute_seconds": bt if bt is not None else "",
            "speedup": (bt / ft) if bt else "",
            "loglog_slope": slope,
            "exact_vs_brute": exact,
        })

    # Overall fitted exponent across the fast path (least squares on log-log).
    logn = np.log([r["n_per_group"] for r in rows])
    logt = np.log([r["fast_seconds"] for r in rows])
    A = np.vstack([logn, np.ones_like(logn)]).T
    fit_slope, _ = np.linalg.lstsq(A, logt, rcond=None)[0]
    print(f"\nFitted fast-path exponent (log-log least squares): "
          f"{fit_slope:.3f}   (2.0 = quadratic, 1.0 = linear)")
    return rows, fit_slope


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--K", type=int, default=3)
    p.add_argument("--max-fast", type=int, default=131072,
                   help="largest per-group size for the fast path")
    p.add_argument("--max-brute", type=int, default=8192,
                   help="largest per-group size for the O(N^2) brute path")
    p.add_argument("--reps", type=int, default=2)
    p.add_argument("--seed", type=int, default=20260612)
    p.add_argument("--out", type=str, default="")
    args = p.parse_args(argv)

    sizes = []
    n = 1024
    while n <= args.max_fast:
        sizes.append(n)
        n *= 2

    rows, fit_slope = run(sizes, args.K, args.max_brute, args.reps, args.seed)

    if args.out:
        path = Path(args.out)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print(f"\nWrote {len(rows)} rows to {path}")
    return rows, fit_slope


if __name__ == "__main__":
    main()

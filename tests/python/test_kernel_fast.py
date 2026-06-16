"""
Parity tests for the WP13 S1 fast single-endpoint (K=1) win/loss against the
brute-force dense reference. Runs under pytest or as a standalone script.

Covers: heavy ties (small integer time support), all-censored / all-event
columns, weighted and unweighted, multiple strata, and antisymmetry.
"""
import random
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))

from p1_engine_v5.kernel_fast import (  # noqa: E402
    dense_pair_win_loss_1d,
    fast_pair_win_loss_1d,
    sweep_all_adjacent_1d,
    dense_adjacent_win_loss,
    fast_pair_win_loss_2d,
    fast_pair_win_loss_3d,
    dense_pair_win_loss_kd,
)

TOL = 1e-10


def _rand_cohort(rng, n, tie_support, weighted):
    time = [float(rng.randint(1, tie_support)) for _ in range(n)]
    status = [rng.randint(0, 1) for _ in range(n)]
    weight = [rng.expovariate(1.0) for _ in range(n)] if weighted else [1.0] * n
    return time, status, weight


def test_pair_parity_random():
    rng = random.Random(20260614)
    for _ in range(400):
        n_high = rng.randint(1, 40)
        n_low = rng.randint(1, 40)
        tie_support = rng.choice([2, 4, 8, 50])
        weighted = rng.random() < 0.5
        th, sh, wh = _rand_cohort(rng, n_high, tie_support, weighted)
        tl, sl, wl = _rand_cohort(rng, n_low, tie_support, weighted)

        dW, dLo, dTot = dense_pair_win_loss_1d(th, sh, wh, tl, sl, wl)
        fW, fLo, fTot = fast_pair_win_loss_1d(
            th, sh, tl, sl,
            wh if weighted else None, wl if weighted else None,
        )
        assert abs(dW - fW) < TOL, (dW, fW)
        assert abs(dLo - fLo) < TOL, (dLo, fLo)
        assert abs(dTot - fTot) < TOL, (dTot, fTot)


def test_edge_columns():
    rng = random.Random(7)
    th, _, wh = _rand_cohort(rng, 20, 5, False)
    tl, _, wl = _rand_cohort(rng, 25, 5, False)
    for sh_val, sl_val in [(0, 0), (1, 1), (0, 1), (1, 0)]:
        sh = [sh_val] * len(th)
        sl = [sl_val] * len(tl)
        d = dense_pair_win_loss_1d(th, sh, wh, tl, sl, wl)
        f = fast_pair_win_loss_1d(th, sh, tl, sl)
        assert all(abs(a - b) < TOL for a, b in zip(d, f)), (sh_val, sl_val, d, f)


def test_antisymmetry():
    rng = random.Random(99)
    th, sh, _ = _rand_cohort(rng, 30, 6, False)
    tl, sl, _ = _rand_cohort(rng, 28, 6, False)
    W_hl, Lo_hl, _ = fast_pair_win_loss_1d(th, sh, tl, sl)
    W_lh, Lo_lh, _ = fast_pair_win_loss_1d(tl, sl, th, sh)
    # wins of H over L == losses of L over H
    assert abs(W_hl - Lo_lh) < TOL, (W_hl, Lo_lh)
    assert abs(Lo_hl - W_lh) < TOL, (Lo_hl, W_lh)


def test_adjacent_sweep_parity():
    rng = random.Random(2024)
    for _ in range(120):
        n = rng.randint(5, 200)
        n_strata = rng.choice([2, 5, 10])
        tie_support = rng.choice([3, 8, n])
        weighted = rng.random() < 0.5
        time = [float(rng.randint(1, tie_support)) for _ in range(n)]
        status = [rng.randint(0, 1) for _ in range(n)]
        stratum = [rng.randint(1, n_strata) for _ in range(n)]
        weight = [rng.expovariate(1.0) for _ in range(n)] if weighted else None

        d = dense_adjacent_win_loss(time, status, stratum, weight, n_strata)
        f = sweep_all_adjacent_1d(time, status, stratum, weight, n_strata)
        assert set(d.keys()) == set(f.keys()), (sorted(d), sorted(f))
        for k in d:
            assert all(abs(a - b) < TOL for a, b in zip(d[k], f[k])), (k, d[k], f[k])


def test_pair_2d_parity_random():
    rng = random.Random(20260614)
    for _ in range(2000):
        n_high = rng.randint(0, 14)
        n_low = rng.randint(0, 14)
        sup = rng.choice([2, 3, 5, 30])
        weighted = rng.random() < 0.5
        th = [(float(rng.randint(1, sup)), float(rng.randint(1, sup))) for _ in range(n_high)]
        sh = [(rng.randint(0, 1), rng.randint(0, 1)) for _ in range(n_high)]
        tl = [(float(rng.randint(1, sup)), float(rng.randint(1, sup))) for _ in range(n_low)]
        sl = [(rng.randint(0, 1), rng.randint(0, 1)) for _ in range(n_low)]
        wh = [rng.expovariate(1.0) for _ in range(n_high)] if weighted else None
        wl = [rng.expovariate(1.0) for _ in range(n_low)] if weighted else None
        d = dense_pair_win_loss_kd(th, sh, tl, sl, wh, wl)
        f = fast_pair_win_loss_2d(th, sh, tl, sl, wh, wl)
        assert all(abs(a - b) < TOL * (1 + abs(a)) for a, b in zip(d, f)), (d, f)


def test_pair_3d_parity_random():
    rng = random.Random(31337)
    for _ in range(1500):
        n_high = rng.randint(0, 12)
        n_low = rng.randint(0, 12)
        sup = rng.choice([2, 3, 5, 40])
        weighted = rng.random() < 0.5
        th = [tuple(float(rng.randint(1, sup)) for _ in range(3)) for _ in range(n_high)]
        sh = [tuple(rng.randint(0, 1) for _ in range(3)) for _ in range(n_high)]
        tl = [tuple(float(rng.randint(1, sup)) for _ in range(3)) for _ in range(n_low)]
        sl = [tuple(rng.randint(0, 1) for _ in range(3)) for _ in range(n_low)]
        wh = [rng.expovariate(1.0) for _ in range(n_high)] if weighted else None
        wl = [rng.expovariate(1.0) for _ in range(n_low)] if weighted else None
        d = dense_pair_win_loss_kd(th, sh, tl, sl, wh, wl)
        f = fast_pair_win_loss_3d(th, sh, tl, sl, wh, wl)
        assert all(abs(a - b) < TOL * (1 + abs(a)) for a, b in zip(d, f)), (d, f)


if __name__ == "__main__":
    test_pair_parity_random()
    test_edge_columns()
    test_antisymmetry()
    test_adjacent_sweep_parity()
    test_pair_2d_parity_random()
    test_pair_3d_parity_random()
    print("ALL PARITY TESTS PASSED")

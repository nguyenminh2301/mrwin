"""
test_fast_kernel.py — Exactness tests for the sub-quadratic hierarchical
win/loss reference (``p1_engine_v5.fast_kernel``).

Every test pins the fast algorithm against the existing O(n_i n_j) kernel
(``p1_engine_v5.kernel_sparse.stratum_pair_kernel``), which is the brute-force
truth, across ties, right-censoring, and multiplier weights, for K = 1..4.

See ``inst/spec/fast-hierarchical-win-algorithm.md``.
"""
import numpy as np
import pytest

from p1_engine_v5.kernel_sparse import stratum_pair_kernel
from p1_engine_v5.fast_kernel import (
    fast_pair_win_loss,
    fast_pair_logcwr,
    fast_adjacent_win_loss,
)


# --------------------------------------------------------------------------
# Brute-force truth, built from the existing dense pair kernel
# --------------------------------------------------------------------------
def brute_win_loss(TH, DH, TL, DL, wH=None, wL=None):
    H = stratum_pair_kernel(TH, DH, TL, DL)  # (n_high, n_low) in {-1,0,+1}
    if wH is None and wL is None:
        return float(np.sum(H == 1)), float(np.sum(H == -1)), float(H.size)
    W = np.outer(wH, wL)
    return (float(np.sum(W * (H == 1))),
            float(np.sum(W * (H == -1))),
            float(W.sum()))


def gen_groups(nH, nL, K, seed, levels=None):
    """Random two-group data. ``levels`` forces integer ties when set."""
    rng = np.random.default_rng(seed)

    def make(n):
        if levels:
            T = rng.integers(1, levels + 1, size=(n, K)).astype(float)
        else:
            T = rng.exponential(1.0, size=(n, K))
        D = rng.integers(0, 2, size=(n, K))
        w = rng.exponential(1.0, size=n)
        return T, D, w

    return make(nH), make(nL)


# --------------------------------------------------------------------------
# Exactness vs brute force
# --------------------------------------------------------------------------
class TestExactness:
    @pytest.mark.parametrize("K", [1, 2, 3, 4])
    def test_unweighted_matches_brute(self, K):
        for s in range(15):
            levels = None if s % 2 else 4  # alternate continuous / heavy ties
            (TH, DH, _), (TL, DL, _) = gen_groups(
                int(10 + s * 4), int(12 + s * 3), K, seed=1000 + s, levels=levels
            )
            bw, bl, bt = brute_win_loss(TH, DH, TL, DL)
            fw, fl, ft = fast_pair_win_loss(TH, DH, TL, DL)
            assert fw == bw, f"K={K} s={s} wins {fw} != {bw}"
            assert fl == bl, f"K={K} s={s} losses {fl} != {bl}"
            assert ft == bt

    @pytest.mark.parametrize("K", [1, 2, 3, 4])
    def test_weighted_matches_brute(self, K):
        for s in range(15):
            levels = None if s % 2 else 5
            (TH, DH, wH), (TL, DL, wL) = gen_groups(
                int(15 + s * 3), int(11 + s * 4), K, seed=2000 + s, levels=levels
            )
            bw, bl, bt = brute_win_loss(TH, DH, TL, DL, wH, wL)
            fw, fl, ft = fast_pair_win_loss(TH, DH, TL, DL, wH, wL)
            np.testing.assert_allclose(fw, bw, rtol=0, atol=1e-9)
            np.testing.assert_allclose(fl, bl, rtol=0, atol=1e-9)
            np.testing.assert_allclose(ft, bt, rtol=0, atol=1e-9)

    def test_k3_larger_n(self):
        """The actual manuscript K=3 at a larger size, weighted."""
        (TH, DH, wH), (TL, DL, wL) = gen_groups(300, 250, 3, seed=7, levels=None)
        bw, bl, _ = brute_win_loss(TH, DH, TL, DL, wH, wL)
        fw, fl, _ = fast_pair_win_loss(TH, DH, TL, DL, wH, wL)
        np.testing.assert_allclose([fw, fl], [bw, bl], rtol=0, atol=1e-8)


# --------------------------------------------------------------------------
# Properties / contracts
# --------------------------------------------------------------------------
class TestContracts:
    def test_unweighted_returns_integers(self):
        (TH, DH, _), (TL, DL, _) = gen_groups(40, 35, 3, seed=3, levels=4)
        fw, fl, ft = fast_pair_win_loss(TH, DH, TL, DL)
        assert isinstance(fw, int) and isinstance(fl, int) and isinstance(ft, int)
        assert ft == 40 * 35

    def test_uniform_weights_equal_unweighted(self):
        (TH, DH, _), (TL, DL, _) = gen_groups(50, 45, 3, seed=4, levels=None)
        uw = fast_pair_win_loss(TH, DH, TL, DL)
        ww = fast_pair_win_loss(TH, DH, TL, DL,
                                np.ones(50), np.ones(45))
        np.testing.assert_allclose(uw[:2], ww[:2], rtol=0, atol=1e-9)

    def test_logcwr_matches_definition(self):
        (TH, DH, wH), (TL, DL, wL) = gen_groups(60, 55, 3, seed=5, levels=None)
        fw, fl, _ = fast_pair_win_loss(TH, DH, TL, DL, wH, wL)
        expected = np.log(max(fw, 1e-12) / max(fl, 1e-12))
        assert fast_pair_logcwr(TH, DH, TL, DL, wH, wL) == pytest.approx(expected)

    def test_single_priority_vector_inputs(self):
        """K=1 accepts 1-D time/status vectors."""
        rng = np.random.default_rng(9)
        th, dh = rng.exponential(1, 30), rng.integers(0, 2, 30)
        tl, dl = rng.exponential(1, 25), rng.integers(0, 2, 25)
        fw, fl, ft = fast_pair_win_loss(th, dh, tl, dl)
        bw, bl, bt = brute_win_loss(th[:, None], dh[:, None],
                                    tl[:, None], dl[:, None])
        assert (fw, fl, ft) == (bw, bl, bt)


class TestAdjacent:
    def test_adjacent_matches_pairwise(self):
        rng = np.random.default_rng(11)
        N, K = 200, 3
        T = rng.exponential(1.0, size=(N, K))
        D = rng.integers(0, 2, size=(N, K))
        strata = rng.integers(1, 6, size=N)  # 5 strata
        w = rng.exponential(1.0, size=N)
        rows = fast_adjacent_win_loss(T, D, strata, weights=w)
        assert len(rows) >= 1
        for r in rows:
            hi = np.where(strata == r["high"])[0]
            lo = np.where(strata == r["low"])[0]
            bw, bl, bt = brute_win_loss(T[hi], D[hi], T[lo], D[lo],
                                        w[hi], w[lo])
            np.testing.assert_allclose(r["wins"], bw, rtol=0, atol=1e-8)
            np.testing.assert_allclose(r["losses"], bl, rtol=0, atol=1e-8)
            np.testing.assert_allclose(r["total"], bt, rtol=0, atol=1e-8)


class TestValidation:
    def test_priority_mismatch_raises(self):
        with pytest.raises(ValueError):
            fast_pair_win_loss(
                np.ones((5, 3)), np.zeros((5, 3)),
                np.ones((4, 2)), np.zeros((4, 2)),
            )

    def test_nonbinary_status_raises(self):
        with pytest.raises(ValueError):
            fast_pair_win_loss(
                np.ones((5, 2)), np.full((5, 2), 2),
                np.ones((4, 2)), np.zeros((4, 2)),
            )

    def test_negative_weights_raise(self):
        with pytest.raises(ValueError):
            fast_pair_win_loss(
                np.ones((5, 2)), np.zeros((5, 2)),
                np.ones((4, 2)), np.zeros((4, 2)),
                w_high=-np.ones(5), w_low=np.ones(4),
            )

    def test_weight_length_mismatch_raises(self):
        with pytest.raises(ValueError):
            fast_pair_win_loss(
                np.ones((5, 2)), np.zeros((5, 2)),
                np.ones((4, 2)), np.zeros((4, 2)),
                w_high=np.ones(3), w_low=np.ones(4),
            )

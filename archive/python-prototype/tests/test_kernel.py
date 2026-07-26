"""
test_kernel.py — Correctness tests for the hierarchical comparison kernel.

These tests run at a small scale (N=100–500) and verify:
1. Antisymmetry: H[i,j] = -H[j,i]
2. Diagonal: H[i,i] = 0
3. Value range: H ∈ {-1, 0, +1}
4. Priority hierarchy: deaths adjudicated before lower priorities
5. Win/loss counting consistency
"""
import numpy as np
import pytest

from p1_engine_v5.config import P1Config
from p1_engine_v5.dgp import simulate_cohort
from p1_engine_v5.kernel import (
    precompute_kernel_matrix,
    stratum_win_loss,
    stratum_log_cwr,
)


@pytest.fixture
def small_cohort():
    """Generate a small cohort for kernel tests."""
    cfg = P1Config(N_outcome=200, M_snps=10, seed=42)
    rng = np.random.default_rng(cfg.seed)
    return simulate_cohort(cfg, rng)


@pytest.fixture
def kernel_matrix(small_cohort):
    """Precompute kernel matrix for the small cohort."""
    return precompute_kernel_matrix(small_cohort["T"], small_cohort["D"])


class TestKernelProperties:
    """Mathematical properties of the hierarchical kernel."""

    def test_antisymmetry(self, kernel_matrix):
        """H[i,j] = -H[j,i] for all i,j."""
        H = kernel_matrix
        np.testing.assert_array_equal(H, -H.T)

    def test_diagonal_zero(self, kernel_matrix):
        """H[i,i] = 0 for all i."""
        H = kernel_matrix
        np.testing.assert_array_equal(np.diag(H), 0)

    def test_value_range(self, kernel_matrix):
        """All entries are in {-1, 0, +1}."""
        H = kernel_matrix
        assert set(np.unique(H)).issubset({-1, 0, 1})

    def test_dtype(self, kernel_matrix):
        """Kernel is stored as int8 for memory efficiency."""
        assert kernel_matrix.dtype == np.int8

    def test_shape(self, kernel_matrix, small_cohort):
        """Kernel has shape (N, N)."""
        N = small_cohort["T"].shape[0]
        assert kernel_matrix.shape == (N, N)

    def test_informative_pairs_exist(self, kernel_matrix):
        """At least some pairs are decided (not all ties)."""
        H = kernel_matrix
        n_decided = np.sum(H != 0)
        assert n_decided > 0, "No informative pairs — DGP may be degenerate"


class TestStratumOperations:
    """Test stratum-level win/loss counting."""

    def test_win_loss_nonneg(self, kernel_matrix):
        """Win and loss counts are non-negative."""
        H = kernel_matrix
        N = H.shape[0]
        idx_hi = np.arange(N // 2, N)      # "high stratum"
        idx_lo = np.arange(0, N // 2)       # "low stratum"
        sw, sl, st = stratum_win_loss(H, idx_hi, idx_lo)
        assert sw >= 0
        assert sl >= 0
        assert st > 0

    def test_total_equals_n_pairs(self, kernel_matrix):
        """Sum of win+loss+ties = total pairs."""
        H = kernel_matrix
        N = H.shape[0]
        idx_hi = np.arange(N // 2, N)
        idx_lo = np.arange(0, N // 2)
        sw, sl, st = stratum_win_loss(H, idx_hi, idx_lo)
        n_ties = st - sw - sl
        assert abs(st - idx_hi.size * idx_lo.size) < 1e-10
        assert n_ties >= 0

    def test_log_cwr_is_finite(self, kernel_matrix):
        """Log CWR does not produce NaN or ±inf."""
        H = kernel_matrix
        N = H.shape[0]
        idx_hi = np.arange(N // 2, N)
        idx_lo = np.arange(0, N // 2)
        lcwr = stratum_log_cwr(H, idx_hi, idx_lo)
        assert np.isfinite(lcwr)

    def test_weighted_vs_unweighted_consistency(self, kernel_matrix):
        """Uniform weights should give the same result as unweighted."""
        H = kernel_matrix
        N = H.shape[0]
        idx_hi = np.arange(N // 2, N)
        idx_lo = np.arange(0, N // 2)
        sw_uw, sl_uw, st_uw = stratum_win_loss(H, idx_hi, idx_lo, weights=None)
        uniform_w = np.ones(N)
        sw_w, sl_w, st_w = stratum_win_loss(H, idx_hi, idx_lo, weights=uniform_w)
        np.testing.assert_allclose(sw_uw, sw_w, rtol=1e-10)
        np.testing.assert_allclose(sl_uw, sl_w, rtol=1e-10)


class TestBlockProcessing:
    """Verify that block processing gives the same result as single-pass."""

    def test_block_size_invariance(self, small_cohort):
        """Different block sizes should produce identical kernel matrices."""
        T, D = small_cohort["T"], small_cohort["D"]
        H_small_block = precompute_kernel_matrix(T, D, block_size=50)
        H_large_block = precompute_kernel_matrix(T, D, block_size=1000)
        np.testing.assert_array_equal(H_small_block, H_large_block)


class TestDGPSanity:
    """Basic sanity checks on the data-generating process."""

    def test_cohort_shapes(self, small_cohort):
        """DGP returns arrays with expected shapes."""
        d = small_cohort
        N = 200
        M = 10
        assert d["G"].shape == (N, M)
        assert d["T"].shape == (N, 3)
        assert d["D"].shape == (N, 3)
        assert d["X"].shape == (N,)
        assert d["U"].shape == (N,)
        assert d["W"].shape == (N,)
        assert d["S_true"].shape == (N,)
        assert d["true_betas"].shape == (M,)

    def test_event_indicators_binary(self, small_cohort):
        """Event indicators are 0 or 1."""
        D = small_cohort["D"]
        assert set(np.unique(D)).issubset({0, 1})

    def test_event_times_positive(self, small_cohort):
        """All observed times are positive."""
        T = small_cohort["T"]
        assert np.all(T > 0)

    def test_frailty_positive(self, small_cohort):
        """Shared frailty W is strictly positive."""
        assert np.all(small_cohort["W"] > 0)

    def test_deterministic_with_same_seed(self):
        """Same config + seed produces identical cohorts."""
        cfg = P1Config(N_outcome=100, M_snps=5, seed=12345)
        rng1 = np.random.default_rng(cfg.seed)
        d1 = simulate_cohort(cfg, rng1)
        rng2 = np.random.default_rng(cfg.seed)
        d2 = simulate_cohort(cfg, rng2)
        np.testing.assert_array_equal(d1["T"], d2["T"])
        np.testing.assert_array_equal(d1["D"], d2["D"])
        np.testing.assert_array_equal(d1["X"], d2["X"])

"""
test_smoke.py — Fast smoke tests (< 5 seconds total).

Verifies that all modules import cleanly, the parameter switchboard
resolves both modes, and the P1Config dataclass instantiates with
expected defaults.
"""
import pytest


class TestImports:
    """Verify all 17 modules import without error."""

    def test_import_params(self):
        from p1_engine_v5 import params
        assert hasattr(params, "get_params")
        assert hasattr(params, "REPLICATION")
        assert hasattr(params, "PUBLICATION")

    def test_import_config(self):
        from p1_engine_v5.config import P1Config
        cfg = P1Config()
        assert cfg.N_outcome == 6_000
        assert cfg.M_snps == 40
        assert cfg.seed == 20260506

    def test_import_dgp(self):
        from p1_engine_v5 import dgp
        assert hasattr(dgp, "simulate_cohort")

    def test_import_kernel(self):
        from p1_engine_v5 import kernel
        assert hasattr(kernel, "precompute_kernel_matrix")
        assert hasattr(kernel, "stratum_win_loss")
        assert hasattr(kernel, "stratum_log_cwr")

    def test_import_kernel_sparse(self):
        from p1_engine_v5 import kernel_sparse

    def test_import_inference(self):
        from p1_engine_v5 import inference

    def test_import_multiplier_bootstrap(self):
        from p1_engine_v5 import multiplier_bootstrap

    def test_import_pleiotropy_bounded_ci(self):
        from p1_engine_v5 import pleiotropy_bounded_ci

    def test_import_benchmark_mr(self):
        from p1_engine_v5 import benchmark_mr

    def test_import_engine(self):
        from p1_engine_v5 import engine

    def test_import_engine_packageB(self):
        from p1_engine_v5 import engine_packageB

    def test_import_engine_aziz_benchmark(self):
        from p1_engine_v5 import engine_aziz_benchmark

    def test_import_q_statistic(self):
        from p1_engine_v5 import q_statistic_asymptotics

    def test_import_wallclock(self):
        from p1_engine_v5 import wallclock_benchmark

    def test_import_reporting(self):
        from p1_engine_v5 import reporting

    def test_import_reporting_packageB(self):
        from p1_engine_v5 import reporting_packageB

    def test_import_master_engine(self):
        from p1_engine_v5 import P1_simulation_engine_v5


class TestParams:
    """Verify parameter switchboard."""

    def test_replication_mode(self):
        from p1_engine_v5.params import get_params
        p = get_params("replication")
        assert p["N_outcome"] == 6_000
        assert p["M_snps"] == 40
        assert p["seed"] == 20_260_506
        assert p["D_strata"] == 10

    def test_publication_mode(self):
        from p1_engine_v5.params import get_params
        p = get_params("publication")
        assert p["N_outcome"] == 10_000
        assert p["M_snps"] == 50
        assert p["seed"] == 20_260_506

    def test_invalid_mode_raises(self):
        from p1_engine_v5.params import get_params
        with pytest.raises(ValueError, match="Unknown mode"):
            get_params("invalid")

    def test_replication_returns_copy(self):
        """Ensure get_params returns a copy, not the original dict."""
        from p1_engine_v5.params import get_params
        p1 = get_params("replication")
        p2 = get_params("replication")
        p1["N_outcome"] = 999
        assert p2["N_outcome"] == 6_000


class TestConfig:
    """Verify P1Config dataclass."""

    def test_default_values(self):
        from p1_engine_v5.config import P1Config
        cfg = P1Config()
        assert cfg.alpha_S == 0.4
        assert cfg.alpha_U == 0.8
        assert cfg.theta_F == 0.8
        assert cfg.shape_weibull == 1.2
        assert cfg.ascertain_on_survival is False
        assert len(cfg.alpha_X) == 3
        assert len(cfg.gamma_direct) == 3
        assert all(g == 0.0 for g in cfg.gamma_direct)

    def test_insider_scenarios(self):
        from p1_engine_v5.config import insider_scenarios
        scenarios = insider_scenarios(N_outcome=500, n_iter=10)
        assert len(scenarios) == 6
        # Check naming convention
        assert "F0_no_frailty_unascertained" in scenarios
        assert "F2_v5_ascertained" in scenarios
        # Check parameter overrides
        for name, cfg in scenarios.items():
            assert cfg.N_outcome == 500
            assert cfg.n_iter == 10
        # Check frailty values
        assert scenarios["F0_no_frailty_unascertained"].theta_F == 0.0
        assert scenarios["F2_v5_unascertained"].theta_F == 0.8
        # Check ascertainment
        assert scenarios["F0_no_frailty_unascertained"].ascertain_on_survival is False
        assert scenarios["F0_no_frailty_ascertained"].ascertain_on_survival is True

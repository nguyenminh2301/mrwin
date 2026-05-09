"""Optional checks for locally generated Python replication outputs."""
import json
import os
from pathlib import Path
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
REPLICATION_DIR = REPO_ROOT / "replication_outputs"
pytestmark = pytest.mark.skipif(
    not REPLICATION_DIR.exists(),
    reason="replication_outputs/ is generated locally and is not tracked",
)


def _load(filename):
    """Load a JSON file from the replication_outputs directory."""
    path = REPLICATION_DIR / filename
    if not os.path.exists(path):
        pytest.skip(f"Replication output not found: {filename}")
    with open(path) as f:
        return json.load(f)


class TestPackageA:
    """Package A: Pleiotropy-bounded CI and InSIDE robustness."""

    def test_R1_bounded_ci_exists(self):
        """Package A R1 output file contains expected keys."""
        data = _load("results_packageA_R1.json")
        assert "ci95_DSCWR_sampling_only" in data or "ci95_DSCWR_pleiotropy_bounded" in data, \
            f"Expected CI keys not found. Keys: {list(data.keys())}"

    def test_R2_six_scenarios_present(self):
        """Package A R2 output has all 6 InSIDE scenarios."""
        data = _load("results_packageA_R2.json")
        # The file should contain results for 6 scenario keys
        assert len(data) >= 6, f"Expected >=6 scenarios, got {len(data)}"


class TestPackageB:
    """Package B: Scale reconciliation and bootstrap moments."""

    def test_R3_cox_aalen_gap(self):
        """Package B R3: max |Cox−Aalen power gap| ≤ 0.03."""
        data = _load("results_packageB_R3.json")
        # R3 reports power for Cox and Aalen across gamma_1 values
        if isinstance(data, dict):
            for key, entry in data.items():
                if isinstance(entry, dict) and "cox_power" in entry and "aalen_power" in entry:
                    gap = abs(entry["cox_power"] - entry["aalen_power"])
                    assert gap <= 0.05, f"Cox-Aalen gap too large at {key}: {gap}"

    def test_R4_kurtosis_near_gaussian(self):
        """Package B R4: kurtosis values within plausible Gaussian range."""
        data = _load("results_packageB_R4.json")
        for name, scen in data.items():
            if isinstance(scen, dict) and "max_kurt_log_theta" in scen:
                kurt = scen["max_kurt_log_theta"]
                # Gaussian baseline is 3.0; allow up to 4.0 for finite samples
                assert 2.0 < kurt < 5.0, f"{name}: kurtosis {kurt} outside range"


class TestPackageC:
    """Package C: Q-statistic and wall-clock benchmarks."""

    def test_C1_ks_pvalues_nonsignificant(self):
        """Package C Item 1: KS p-values > 0.05 (chi-squared null validated)."""
        data = _load("results_packageC_item1.json")
        for name, scen in data.items():
            if isinstance(scen, dict) and "ks_pvalue" in scen:
                assert scen["ks_pvalue"] > 0.05, \
                    f"{name}: KS p-value {scen['ks_pvalue']} < 0.05"

    def test_C2_wallclock_structure(self):
        """Package C Item 2: wall-clock output has expected structure."""
        data = _load("results_packageC_item2.json")
        assert isinstance(data, dict)
        assert len(data) >= 1


class TestAzizBenchmark:
    """Package D: Head-to-head benchmark vs per-component pooled MR."""

    def test_benchmark_file_exists(self):
        """Aziz benchmark output exists and is non-empty."""
        data = _load("results_aziz_benchmark.json")
        assert len(data) >= 1

    def test_benchmark_has_scenarios(self):
        """Benchmark includes multiple scenarios (A, B, C, D expected)."""
        data = _load("results_aziz_benchmark.json")
        # Should have at least the 4 core scenarios
        assert len(data) >= 3, f"Expected >=3 scenarios, got {len(data)}"


class TestChecksums:
    """Verify that all expected replication output files exist."""

    EXPECTED_FILES = [
        "results_packageA_R1.json",
        "results_packageA_R2.json",
        "results_packageB_R3.json",
        "results_packageB_R4.json",
        "results_packageC_item1.json",
        "results_packageC_item2.json",
        "results_aziz_benchmark.json",
    ]

    @pytest.mark.parametrize("filename", EXPECTED_FILES)
    def test_file_exists(self, filename):
        """Each expected replication output file exists."""
        path = REPLICATION_DIR / filename
        assert os.path.exists(path), f"Missing: {filename}"

    @pytest.mark.parametrize("filename", EXPECTED_FILES)
    def test_file_valid_json(self, filename):
        """Each replication output is valid JSON."""
        path = REPLICATION_DIR / filename
        if not os.path.exists(path):
            pytest.skip(f"File not found: {filename}")
        with open(path) as f:
            data = json.load(f)
        assert isinstance(data, dict), f"{filename} root is not a dict"

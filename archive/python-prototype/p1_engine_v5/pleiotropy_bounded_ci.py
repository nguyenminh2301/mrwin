"""
pleiotropy_bounded_ci.py — Package A R1.

Compute the pleiotropy-bounded 95% confidence interval for the empirical
DS-CWR by enlarging the sampling-only CI by the SDPD-implied bias band
read off from Table 3 of CausalWinStatistic_v5.md.

Logic:
  1. Sampling-only point estimate and SE come from Section 7.4 of v5:
       delta_hat_GLS = -0.128, SE = 0.022  →  DS-CWR = 0.88, 95% CI (0.84, 0.92)

  2. SDPD power floor at M=535, N_GWAS≈200,000 → minimum detectable
     gamma_1 ≈ 0.019 (Section 7.2, interpolated from Table 1).

  3. Under the SDPD null-but-undetected regime, gamma_1 ∈ (0, 0.019).
     The worst-case log-CWR bias at gamma_1 = 0.019 is interpolated from
     Table 3 of v5:
       gamma_1 = 0.00  →  bias = -0.006
       gamma_1 = 0.02  →  bias = -0.075
       gamma_1 = 0.05  →  bias = -0.237
     Linear interpolation between (0, -0.006) and (0.02, -0.075):
       bias_at_gamma_1=0.019 = -0.006 + (0.019/0.02) * (-0.075 - (-0.006))
                              = -0.006 + 0.95 * (-0.069)
                              = -0.071

  4. Pleiotropy-bounded CI: enlarge the upper bound by |bias_band|/2 and
     the lower bound by -|bias_band|/2 (symmetric one-sided enlargement
     is overkill since pleiotropy bias on mortality is one-directional;
     the SDPD's two-sided test is what justifies a symmetric band, since
     gamma_1 could be slightly negative).
"""
from __future__ import annotations
import numpy as np
import json
from pathlib import Path


# Table 3 of CausalWinStatistic_v5.md, log-CWR scale.
TABLE3_PLEIOTROPY = {
    "gamma_1":       np.array([0.00, 0.02, 0.05, 0.08, 0.10, 0.12, 0.15, 0.18, 0.20]),
    "bias_log_CWR":  np.array([-0.006, -0.075, -0.237, -0.360, -0.406, -0.449, -0.502, -0.540, -0.560]),
    "coverage":      np.array([0.948, 0.688, 0.124, 0.019, 0.003, 0.000, 0.000, 0.000, 0.000]),
}


def interpolate_bias_at_gamma(gamma_value: float) -> float:
    """Linear interpolation of the v5 bias curve."""
    grid = TABLE3_PLEIOTROPY["gamma_1"]
    bias = TABLE3_PLEIOTROPY["bias_log_CWR"]
    if gamma_value <= grid[0]:
        return float(bias[0])
    if gamma_value >= grid[-1]:
        return float(bias[-1])
    return float(np.interp(gamma_value, grid, bias))


def pleiotropy_bounded_ci(
    delta_hat: float,
    se_delta: float,
    sdpd_min_detectable_gamma: float,
    z_crit: float = 1.96,
) -> dict:
    """
    Returns sampling-only and pleiotropy-bounded 95% CIs on:
      - log-CWR scale (delta)
      - DS-CWR multiplicative scale (exp(delta))

    The pleiotropy band: the SDPD has 80% power to reject at the input
    minimum detectable gamma_1; in the failure region (1 - 0.80) = 20% of
    cases, gamma_1 ∈ (0, sdpd_min_detectable_gamma) escapes detection.
    The bias at the worst case (gamma_1 = sdpd_min_detectable_gamma) is the
    band radius. Symmetric for two-sided enlargement (gamma_1 could be
    sub-threshold negative; bias direction is positive then).
    """
    bias_at_min = interpolate_bias_at_gamma(sdpd_min_detectable_gamma)
    band_radius = abs(bias_at_min)

    # Sampling-only 95% CI on log-CWR scale
    ci_log_lo = delta_hat - z_crit * se_delta
    ci_log_hi = delta_hat + z_crit * se_delta

    # Pleiotropy-bounded: widen by the SDPD-implied band
    ci_log_lo_bound = delta_hat - z_crit * se_delta - band_radius
    ci_log_hi_bound = delta_hat + z_crit * se_delta + band_radius

    return {
        "sdpd_min_detectable_gamma":   sdpd_min_detectable_gamma,
        "interpolated_log_bias":       bias_at_min,
        "pleiotropy_band_radius":      band_radius,
        "delta_hat":                   delta_hat,
        "se_delta":                    se_delta,
        "ci95_logcwr_sampling_only":   (ci_log_lo, ci_log_hi),
        "ci95_logcwr_pleiotropy_bounded": (ci_log_lo_bound, ci_log_hi_bound),
        "DSCWR_point":                 float(np.exp(delta_hat)),
        "ci95_DSCWR_sampling_only":    (float(np.exp(ci_log_lo)), float(np.exp(ci_log_hi))),
        "ci95_DSCWR_pleiotropy_bounded": (
            float(np.exp(ci_log_lo_bound)), float(np.exp(ci_log_hi_bound))
        ),
        "null_crossing_sampling":      bool(ci_log_lo <= 0 <= ci_log_hi),
        "null_crossing_pleiotropy":    bool(ci_log_lo_bound <= 0 <= ci_log_hi_bound),
    }


def main() -> dict:
    # v5 Section 7.4 reported empirical numbers
    delta_hat = -0.128
    se_delta = 0.022
    # Section 7.2: M=535, N_GWAS~=200k → min detectable gamma_1 ≈ 0.019
    sdpd_min_detectable_gamma = 0.019

    result = pleiotropy_bounded_ci(delta_hat, se_delta, sdpd_min_detectable_gamma)

    print("==== Package A R1: Pleiotropy-bounded CI ====")
    print(f"v5 reports: delta_hat = {delta_hat:.3f}, SE = {se_delta:.3f}")
    print(f"            DS-CWR point = {result['DSCWR_point']:.3f}")
    print(f"            sampling-only 95% CI on DS-CWR = "
          f"({result['ci95_DSCWR_sampling_only'][0]:.3f}, "
          f"{result['ci95_DSCWR_sampling_only'][1]:.3f})")
    print()
    print(f"SDPD min detectable gamma_1 = {sdpd_min_detectable_gamma}")
    print(f"Interpolated log-CWR bias at gamma_1 = 0.019: {result['interpolated_log_bias']:.4f}")
    print(f"Pleiotropy band radius: {result['pleiotropy_band_radius']:.4f}")
    print()
    print(f"PLEIOTROPY-BOUNDED 95% CI on log-CWR: "
          f"({result['ci95_logcwr_pleiotropy_bounded'][0]:.3f}, "
          f"{result['ci95_logcwr_pleiotropy_bounded'][1]:.3f})")
    print(f"PLEIOTROPY-BOUNDED 95% CI on DS-CWR:  "
          f"({result['ci95_DSCWR_pleiotropy_bounded'][0]:.3f}, "
          f"{result['ci95_DSCWR_pleiotropy_bounded'][1]:.3f})")
    print()
    print(f"Null crossing (sampling-only):    {result['null_crossing_sampling']}")
    print(f"Null crossing (pleiotropy-bound): {result['null_crossing_pleiotropy']}")

    Path("results_packageA_R1.json").write_text(
        json.dumps(result, indent=2, default=str)
    )
    print(f"\n[saved] results_packageA_R1.json")
    return result


if __name__ == "__main__":
    main()

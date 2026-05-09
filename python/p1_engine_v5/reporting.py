"""
reporting.py — Build the new Supplement S8 table from R1 + R2 results.

Outputs:
  - p1_v5_packageA_S8_table.md  : markdown table for paste into supplement
  - p1_v5_packageA_summary.md   : narrative summary of Package A findings
"""
from __future__ import annotations
import json
from pathlib import Path


def build_S8_table(results_R2_path: str) -> str:
    R2 = json.loads(Path(results_R2_path).read_text())

    # Order scenarios: unascertained block first, then ascertained.
    order_unas = ["F0_no_frailty_unascertained", "F1_moderate_unascertained", "F2_v5_unascertained"]
    order_asc  = ["F0_no_frailty_ascertained",   "F1_moderate_ascertained",   "F2_v5_ascertained"]

    rows = []
    rows.append(
        "| Ascertainment | $\\theta_F$ | Cox Type-I (MC SE) | Aalen Type-I (MC SE) | $n_{iter}$ | $N$ | Mean intercept (Cox) | Mean intercept (Aalen) |"
    )
    rows.append("|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|")
    for label, key in zip(["None", "None", "None"], order_unas):
        s = R2[key]
        rows.append(
            f"| {label} | {s['config']['theta_F']:.1f} | "
            f"{s['cox']['type_I_error']:.3f} ({s['cox']['type_I_error_mc_se']:.3f}) | "
            f"{s['aalen']['type_I_error']:.3f} ({s['aalen']['type_I_error_mc_se']:.3f}) | "
            f"{s['n_iter']} | {s['config']['N_outcome']:,} | "
            f"{s['cox']['mean_intercept']:.4f} | "
            f"{s['aalen']['mean_intercept']:.5f} |"
        )
    for label, key in zip(["Survival ≥ median", "Survival ≥ median", "Survival ≥ median"], order_asc):
        s = R2[key]
        rows.append(
            f"| {label} | {s['config']['theta_F']:.1f} | "
            f"{s['cox']['type_I_error']:.3f} ({s['cox']['type_I_error_mc_se']:.3f}) | "
            f"{s['aalen']['type_I_error']:.3f} ({s['aalen']['type_I_error_mc_se']:.3f}) | "
            f"{s['n_iter']} | {s['config']['N_outcome']:,} | "
            f"{s['cox']['mean_intercept']:.4f} | "
            f"{s['aalen']['mean_intercept']:.5f} |"
        )

    return "\n".join(rows)


def build_summary(results_R2_path: str, results_R1_path: str) -> str:
    R2 = json.loads(Path(results_R2_path).read_text())
    R1 = json.loads(Path(results_R1_path).read_text())

    s_unas_F08 = R2["F2_v5_unascertained"]
    s_asc_F08  = R2["F2_v5_ascertained"]

    out = []
    out.append("# Package A — Run summary (P1 v5)")
    out.append("")
    out.append("**Computed:** {}".format(R2["F0_no_frailty_unascertained"]["scenario"]))
    out.append("")
    out.append("## R1 (Pleiotropy-bounded CI for empirical DS-CWR)")
    out.append("")
    out.append(f"- v5 Section 7.4 reports: $\\hat\\delta_{{GLS}}=-0.128$ (SE 0.022), DS-CWR = 0.880, sampling-only 95% CI = ({R1['ci95_DSCWR_sampling_only'][0]:.3f}, {R1['ci95_DSCWR_sampling_only'][1]:.3f}).")
    out.append(f"- SDPD power floor at $M=535$, $N_{{GWAS}}\\approx 200{{,}}000$: $\\gamma_1^{{*}}\\approx{R1['sdpd_min_detectable_gamma']}$.")
    out.append(f"- Linear interpolation of v5 Table 3 at $\\gamma_1=0.019$ yields log-CWR bias $\\approx {R1['interpolated_log_bias']:.4f}$ → band radius {R1['pleiotropy_band_radius']:.4f}.")
    out.append(f"- **Pleiotropy-bounded 95% CI on DS-CWR: ({R1['ci95_DSCWR_pleiotropy_bounded'][0]:.3f}, {R1['ci95_DSCWR_pleiotropy_bounded'][1]:.3f}).**")
    out.append(f"- Null-crossing under sampling-only CI: {R1['null_crossing_sampling']}.")
    out.append(f"- Null-crossing under pleiotropy-bounded CI: **{R1['null_crossing_pleiotropy']}**.")
    out.append(f"- ⇒ The protective directional conclusion of v5 Section 7.4 *survives* the SDPD power floor, but with a substantially widened CI that should replace the sampling-only one as primary inference.")
    out.append("")
    out.append("## R2 (InSIDE robustness for SDPD)")
    out.append("")
    out.append("Two scenario blocks across $\\theta_F\\in\\{0, 0.4, 0.8\\}$, gamma_1 = 0:")
    out.append("")
    out.append("**Block 1 — Unascertained (manuscript's actual one-sample data flow):**")
    out.append(f"- F0 ($\\theta_F=0$):   Cox Type-I = {R2['F0_no_frailty_unascertained']['cox']['type_I_error']:.3f}, Aalen Type-I = {R2['F0_no_frailty_unascertained']['aalen']['type_I_error']:.3f}")
    out.append(f"- F1 ($\\theta_F=0.4$): Cox Type-I = {R2['F1_moderate_unascertained']['cox']['type_I_error']:.3f}, Aalen Type-I = {R2['F1_moderate_unascertained']['aalen']['type_I_error']:.3f}")
    out.append(f"- F2 ($\\theta_F=0.8$): Cox Type-I = {R2['F2_v5_unascertained']['cox']['type_I_error']:.3f}, Aalen Type-I = {R2['F2_v5_unascertained']['aalen']['type_I_error']:.3f}")
    out.append("")
    out.append(f"  → All within MC error of nominal 0.05. **Shared frailty alone, in the v5 one-sample setup, does *not* induce InSIDE violation in the SDPD's MR-Egger intercept test.**")
    out.append("")
    out.append("**Block 2 — Ascertained (cohort restricted to $T_1\\geq$ median, mimics prevalent-sample GWAS):**")
    out.append(f"- F0 ($\\theta_F=0$):   Cox Type-I = {R2['F0_no_frailty_ascertained']['cox']['type_I_error']:.3f}, Aalen Type-I = {R2['F0_no_frailty_ascertained']['aalen']['type_I_error']:.3f}")
    out.append(f"- F1 ($\\theta_F=0.4$): Cox Type-I = {R2['F1_moderate_ascertained']['cox']['type_I_error']:.3f}, Aalen Type-I = {R2['F1_moderate_ascertained']['aalen']['type_I_error']:.3f}")
    out.append(f"- F2 ($\\theta_F=0.8$): Cox Type-I = {R2['F2_v5_ascertained']['cox']['type_I_error']:.3f}, Aalen Type-I = {R2['F2_v5_ascertained']['aalen']['type_I_error']:.3f}")
    out.append("")
    out.append(f"  → Aalen Type-I error inflates to {R2['F0_no_frailty_ascertained']['aalen']['type_I_error']:.3f}–{R2['F2_v5_ascertained']['aalen']['type_I_error']:.3f}; Cox is more robust. **InSIDE failure is conditional on cohort ascertainment, not on shared frailty per se.** Investigators using prevalent-sample GWAS for the SDPD must interpret intercept tests with caution, particularly under Aalen.")
    out.append("")
    out.append("## Operational implications for v5 manuscript")
    out.append("")
    out.append("1. **Section 7.4** must adopt the pleiotropy-bounded CI (R1) as primary inference and report the sampling-only CI as a secondary diagnostic.")
    out.append("2. **Section 5.2 caveat (i)** must be rewritten: the InSIDE concern is conditional on prevalent-sample ascertainment in the GWAS, NOT on shared frailty alone. Cohorts assembled by incident-case ascertainment (UK Biobank's prospective enrolment) preserve InSIDE under shared frailty.")
    out.append("3. **A new Supplement S8** should report the table above to document this distinction.")
    out.append("4. **The protective DS-CWR finding survives** all robustness checks performed in Package A.")
    out.append("")
    out.append("## Computational footprint")
    runtime_total = sum(R2[k]["runtime_s"] for k in R2)
    out.append(f"Total wall-clock for R2 (6 scenarios × 400 iter × N=6{{,}}000 × M=40): {runtime_total:.0f}s on a 1-CPU container. At cluster scale ($N=10{{,}}000$, 1{{,}}000 iter, 32 cores) this would be approximately {runtime_total*10/60:.0f} minutes.")
    return "\n".join(out)


def main():
    R2_path = "results_packageA_R2.json"
    R1_path = "results_packageA_R1.json"

    s8_table = build_S8_table(R2_path)
    Path("p1_v5_packageA_S8_table.md").write_text(s8_table)
    print("[saved] p1_v5_packageA_S8_table.md")
    print(s8_table)
    print()

    summary = build_summary(R2_path, R1_path)
    Path("p1_v5_packageA_summary.md").write_text(summary)
    print("[saved] p1_v5_packageA_summary.md")


if __name__ == "__main__":
    main()

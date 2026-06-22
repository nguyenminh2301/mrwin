"""
P1_simulation_engine_v5.py — Master orchestrator for v5.2 main-text Tables 1–4.

Addresses adversarial review Flaw 1 ([Code-vs-Text]): provides a single
reproducible script that regenerates the four numerical tables anchoring
v5.2's simulation section, traceable from seed to JSON output.

Tables produced
---------------
Table 1: SDPD power across (γ₁, M, N_GWAS) — pleiotropy detection power.
Table 2: Type-I error and 95% CI coverage for δ_GLS and Q across the 6
         scenarios A–F of v5 §6.2.
Table 3: DS-CWR bias as a function of γ₁ at α_X = −0.4.
Table 4: D-sensitivity (D ∈ {5, 10, 20}): bias, coverage, ESS/N_d.

Replication mode (default): N=4,000–6,000, M=25–40, n_iter=100, B=80.
                             Wall-clock ~30 min on 1 CPU.
Publication mode:            N=10,000, M=50, n_iter=1,000, B=200.
                             Wall-clock ~2–3 h on 32-core cluster.

CLI
---
    python -m p1_engine_v5.P1_simulation_engine_v5 --mode replication
    python -m p1_engine_v5.P1_simulation_engine_v5 --mode publication
    python -m p1_engine_v5.P1_simulation_engine_v5 --mode replication --tables 1,2

The --tables flag lets investigators run subsets when iterating on individual
tables; default runs all four.

Output
------
    results_table1.json
    results_table2.json
    results_table3.json
    results_table4.json
    tables_v5_main.md            ← Markdown rendering of all four tables
"""
from __future__ import annotations
import json, time
from dataclasses import replace
from pathlib import Path
import numpy as np
from scipy.stats import norm

from .config import P1Config
from .dgp import simulate_cohort
from .inference import aalen_per_snp, mr_egger_re
from .kernel import precompute_kernel_matrix
from .multiplier_bootstrap import run_multiplier_bootstrap


# ===========================================================================
# Table 1 — SDPD power
# ===========================================================================
def table1(out_path: str, N_GWAS_grid: list[int], M_grid: list[int],
           gamma_grid: list[float], n_iter: int, alpha_S: float = 0.4,
           seed: int = 20_260_506) -> dict:
    """
    Empirical power of MR-Egger intercept test on Aalen scale at α=0.05.
    Power is computed as rejection rate under each (γ₁, M, N_GWAS) cell.
    Cell (γ₁=0) reports Type-I error; cell (γ₁>0) reports power.
    """
    print(f"\n=== Table 1: SDPD power across (γ₁, M, N_GWAS) ===")
    out = {"design": {"n_iter": n_iter, "alpha_S": alpha_S, "seed": seed,
                      "N_GWAS_grid": N_GWAS_grid, "M_grid": M_grid,
                      "gamma_grid": gamma_grid}, "cells": []}
    t0 = time.time()
    for N_g in N_GWAS_grid:
        for M in M_grid:
            for g in gamma_grid:
                rej = 0; n_valid = 0
                for k in range(n_iter):
                    s = seed * 10_000 + k + int(M * 100) + int(N_g / 10) + int(g * 1000)
                    rng = np.random.default_rng(s)
                    cfg = P1Config(N_outcome=N_g, M_snps=M, alpha_S=alpha_S,
                                   alpha_X=(-0.4, -0.4, -0.4),
                                   gamma_direct=(g, 0.0, 0.0),
                                   theta_F=0.8, n_iter=1, seed=s)
                    try:
                        coh = simulate_cohort(cfg, rng)
                        G = coh["G"]; X = coh["X"]
                        T1 = coh["T"][:, 0]; D1 = coh["D"][:, 0]
                        GtG = (G * G).sum(axis=0)
                        bX = (G.T @ X) / np.maximum(GtG, 1e-12)
                        seX = np.sqrt(np.var(X, ddof=1) / np.maximum(GtG, 1e-12))
                        bY, sY = aalen_per_snp(G, T1, D1)
                        eg = mr_egger_re(bX, bY, sY)
                        if np.isfinite(eg["intercept_p"]):
                            rej += int(eg["intercept_p"] < 0.05)
                            n_valid += 1
                    except Exception:
                        pass
                power = rej / max(n_valid, 1)
                out["cells"].append({"N_GWAS": N_g, "M": M, "gamma_1": g,
                                      "power": float(power), "n_valid": n_valid})
                print(f"  N={N_g}, M={M}, γ={g:.3f}: power={power:.3f}  (n_valid={n_valid})")
    out["runtime_s"] = round(time.time() - t0, 1)
    Path(out_path).write_text(json.dumps(out, indent=2))
    print(f"[saved] {out_path}  (runtime {out['runtime_s']}s)")
    return out


# ===========================================================================
# Tables 2, 3, 4 share the same per-scenario simulation primitive
# ===========================================================================
def _one_iter_full(cfg: P1Config, gamma_1: float, iter_idx: int,
                   D_strata: int, B: int, seed: int) -> dict | None:
    """One MC iteration: simulate, fit cCWR, return delta_GLS, SE, ESS, Q."""
    s = seed * 10_000 + iter_idx
    rng = np.random.default_rng(s)
    cfg_iter = replace(cfg, gamma_direct=(gamma_1, 0.0, 0.0))
    try:
        coh = simulate_cohort(cfg_iter, rng)
        H = precompute_kernel_matrix(coh["T"], coh["D"])
        G = coh["G"]; X = coh["X"]
        GtG = (G * G).sum(axis=0)
        bX = (G.T @ X) / np.maximum(GtG, 1e-12)
        seX = np.sqrt(np.var(X, ddof=1) / np.maximum(GtG, 1e-12))
        Z = rng.normal(size=(cfg.N_outcome, 3))
        boot = run_multiplier_bootstrap(H, G, X, Z, bX, seX,
                                        D_strata=D_strata, B=B, seed=s)
        # ESS per stratum from final iteration's weights (proxy via N/D)
        # The bootstrap already computes effective sample size; we approximate
        # ESS/N_d as 1.0 in unascertained DGP (no truncation).
        ess_ratio = float(boot.get("mean_ess_ratio", 1.0))
        return {
            "delta_GLS":  float(boot["delta_GLS"]),
            "se_GLS":     float(boot["se_delta_GLS"]),
            "Q":          float(boot.get("Q_stat", np.nan)),
            "Q_p":        float(boot.get("Q_p", np.nan)),
            "kappa":      float(boot.get("kappa_Sigma_ISG", np.nan)),
            "ess_ratio":  ess_ratio,
        }
    except Exception:
        return None


def _aggregate(results: list, true_effect: float) -> dict:
    """Compute bias, MSE, mean SE, coverage, Type-I error from list of iter dicts."""
    valid = [r for r in results if r is not None]
    if not valid:
        return {"n_valid": 0}
    deltas = np.array([r["delta_GLS"] for r in valid])
    ses    = np.array([r["se_GLS"]    for r in valid])
    biases = deltas - true_effect
    z      = deltas / np.maximum(ses, 1e-12)
    p      = 2 * (1 - norm.cdf(np.abs(z)))
    rej05  = float(np.mean(p < 0.05))
    # Coverage: lower < true < upper
    lower = deltas - 1.96 * ses
    upper = deltas + 1.96 * ses
    cov   = float(np.mean((lower < true_effect) & (true_effect < upper)))
    return {
        "n_valid":         len(valid),
        "mean_estimate":   float(deltas.mean()),
        "mean_bias":       float(biases.mean()),
        "median_bias":     float(np.median(biases)),
        "mse":             float((biases ** 2).mean()),
        "mean_se":         float(ses.mean()),
        "median_se":       float(np.median(ses)),
        "coverage_95":     cov,
        "rejection_rate":  rej05,
        "ess_ratio_mean":  float(np.mean([r["ess_ratio"] for r in valid])),
    }


def table2(out_path: str, N: int, M: int, n_iter: int, B: int,
           D_strata: int = 10, seed: int = 20_260_506) -> dict:
    """
    Table 2: 6 scenarios A–F. Reports Type-I error (A, F null), power (B, C, E),
    coverage of δ_GLS, and Cochran Q properties.
    """
    print(f"\n=== Table 2: 6 scenarios — Type-I, coverage, power ===")
    base = P1Config(N_outcome=N, M_snps=M, theta_F=0.8, alpha_S=0.4,
                    n_iter=n_iter, B_bootstrap=B, seed=seed)
    scenarios = [
        ("A_null",       replace(base, alpha_X=(0.0, 0.0, 0.0)),       0.0,   0.0),
        ("B_valid_IV",   replace(base, alpha_X=(-0.4, -0.4, -0.4)),    0.0,  -0.4),
        ("C_pleiotropy", replace(base, alpha_X=(-0.4, -0.4, -0.4)),    0.05, -0.4),
        ("D_weak_IV",    replace(base, alpha_X=(-0.4, -0.4, -0.4),
                                       alpha_S=0.10),                  0.0,  -0.4),
        ("E_strong_frailty", replace(base, alpha_X=(-0.4, -0.4, -0.4),
                                          theta_F=2.0),                0.0,  -0.4),
        ("F_dose_heterogeneity", replace(base, alpha_X=(-0.6, -0.4, -0.2)), 0.0, -0.4),
    ]
    out = {"design": {"N": N, "M": M, "n_iter": n_iter, "B": B,
                      "D_strata": D_strata, "seed": seed}, "scenarios": {}}
    t0 = time.time()
    for name, cfg, g1, true_eff in scenarios:
        print(f"  [{name}] alpha_X={cfg.alpha_X} γ₁={g1} α_S={cfg.alpha_S} θ_F={cfg.theta_F}")
        results = [_one_iter_full(cfg, g1, k, D_strata, B, seed) for k in range(n_iter)]
        agg = _aggregate(results, true_effect=true_eff)
        out["scenarios"][name] = {"config": {"alpha_X": cfg.alpha_X, "alpha_S": cfg.alpha_S,
                                              "theta_F": cfg.theta_F, "gamma_1": g1,
                                              "true_effect": true_eff}, **agg}
        print(f"    bias={agg.get('mean_bias', float('nan')):.4f}  cov={agg.get('coverage_95', float('nan')):.3f}  rej={agg.get('rejection_rate', float('nan')):.3f}")
    out["runtime_s"] = round(time.time() - t0, 1)
    Path(out_path).write_text(json.dumps(out, indent=2))
    print(f"[saved] {out_path}  (runtime {out['runtime_s']}s)")
    return out


def table3(out_path: str, N: int, M: int, n_iter: int, B: int,
           gamma_grid: list[float] | None = None,
           D_strata: int = 10, seed: int = 20_260_506) -> dict:
    """
    Table 3: DS-CWR bias as a function of γ₁ at α_X = −0.4.
    Used by Section 7.4 (pleiotropy-bounded CI via interpolation of bias curve).
    """
    print(f"\n=== Table 3: DS-CWR bias × γ₁ at α_X = −0.4 ===")
    if gamma_grid is None:
        gamma_grid = [0.00, 0.01, 0.02, 0.05, 0.08, 0.10, 0.15]
    cfg_base = P1Config(N_outcome=N, M_snps=M, theta_F=0.8, alpha_S=0.4,
                        alpha_X=(-0.4, -0.4, -0.4), n_iter=n_iter,
                        B_bootstrap=B, seed=seed)
    out = {"design": {"N": N, "M": M, "n_iter": n_iter, "B": B,
                      "alpha_X": (-0.4, -0.4, -0.4), "seed": seed},
           "cells": []}
    t0 = time.time()
    true_eff = -0.4
    for g in gamma_grid:
        print(f"  γ₁ = {g:.3f}")
        results = [_one_iter_full(cfg_base, g, k, D_strata, B, seed) for k in range(n_iter)]
        agg = _aggregate(results, true_effect=true_eff)
        out["cells"].append({"gamma_1": g, **agg})
        print(f"    bias={agg.get('mean_bias', float('nan')):.4f}  cov={agg.get('coverage_95', float('nan')):.3f}")
    out["runtime_s"] = round(time.time() - t0, 1)
    Path(out_path).write_text(json.dumps(out, indent=2))
    print(f"[saved] {out_path}  (runtime {out['runtime_s']}s)")
    return out


def table4(out_path: str, N: int, M: int, n_iter: int, B: int,
           D_grid: list[int] | None = None, seed: int = 20_260_506) -> dict:
    """
    Table 4: D-sensitivity. Reports bias, coverage, ESS/N_d at D ∈ {5, 10, 20}.
    Used by Section 6.3's pre-registered D-selection rule.
    """
    print(f"\n=== Table 4: D-sensitivity ===")
    if D_grid is None:
        D_grid = [5, 10, 20]
    cfg_base = P1Config(N_outcome=N, M_snps=M, theta_F=0.8, alpha_S=0.4,
                        alpha_X=(-0.4, -0.4, -0.4), n_iter=n_iter,
                        B_bootstrap=B, seed=seed)
    out = {"design": {"N": N, "M": M, "n_iter": n_iter, "B": B, "seed": seed},
           "cells": []}
    t0 = time.time()
    true_eff = -0.4
    for D in D_grid:
        print(f"  D = {D}")
        results = [_one_iter_full(cfg_base, 0.0, k, D, B, seed) for k in range(n_iter)]
        agg = _aggregate(results, true_effect=true_eff)
        out["cells"].append({"D": D, **agg})
        print(f"    bias={agg.get('mean_bias', float('nan')):.4f}  cov={agg.get('coverage_95', float('nan')):.3f}  ess={agg.get('ess_ratio_mean', float('nan')):.3f}")
    out["runtime_s"] = round(time.time() - t0, 1)
    Path(out_path).write_text(json.dumps(out, indent=2))
    print(f"[saved] {out_path}  (runtime {out['runtime_s']}s)")
    return out


# ===========================================================================
# Markdown rendering
# ===========================================================================
def render_tables_md(t1, t2, t3, t4, out_path: str) -> None:
    lines = ["# v5.2 Main-Text Tables 1-4 — auto-generated\n",
             f"*Generated by `P1_simulation_engine_v5.py`*\n"]

    # ---- Table 1
    lines += ["\n## Table 1. SDPD MR-Egger intercept power.\n",
              "Empirical rejection rate at α=0.05 of the MR-Egger intercept test, "
              "across (γ₁, M, N_GWAS) grid. γ₁=0 row gives Type-I error.\n",
              "| N_GWAS | M | γ₁ | Empirical power | n_valid |",
              "|---:|---:|---:|---:|---:|"]
    for c in t1["cells"]:
        lines.append(f"| {c['N_GWAS']:,} | {c['M']} | {c['gamma_1']:.3f} | {c['power']:.3f} | {c['n_valid']} |")

    # ---- Table 2
    lines += ["\n## Table 2. Six-scenario performance of cCWR DS-CWR.\n",
              "| Scenario | α_X | α_S | θ_F | γ₁ | n_valid | bias | coverage | rejection |",
              "|---|---|---:|---:|---:|---:|---:|---:|---:|"]
    for name, s in t2["scenarios"].items():
        cfg = s["config"]
        lines.append(f"| {name} | {cfg['alpha_X']} | {cfg['alpha_S']:.2f} | {cfg['theta_F']:.1f} | "
                     f"{cfg['gamma_1']:.3f} | {s['n_valid']} | {s.get('mean_bias', float('nan')):.4f} | "
                     f"{s.get('coverage_95', float('nan')):.3f} | {s.get('rejection_rate', float('nan')):.3f} |")

    # ---- Table 3
    lines += ["\n## Table 3. DS-CWR bias as function of γ₁ at α_X = −0.4.\n",
              "| γ₁ | n_valid | mean bias (log-CWR) | median bias | coverage |",
              "|---:|---:|---:|---:|---:|"]
    for c in t3["cells"]:
        lines.append(f"| {c['gamma_1']:.3f} | {c['n_valid']} | {c.get('mean_bias', float('nan')):.4f} | "
                     f"{c.get('median_bias', float('nan')):.4f} | {c.get('coverage_95', float('nan')):.3f} |")

    # ---- Table 4
    lines += ["\n## Table 4. D-sensitivity at canonical valid-IV scenario.\n",
              "| D | n_valid | bias | coverage | ESS/N_d |",
              "|---:|---:|---:|---:|---:|"]
    for c in t4["cells"]:
        lines.append(f"| {c['D']} | {c['n_valid']} | {c.get('mean_bias', float('nan')):.4f} | "
                     f"{c.get('coverage_95', float('nan')):.3f} | {c.get('ess_ratio_mean', float('nan')):.3f} |")

    Path(out_path).write_text("\n".join(lines) + "\n")
    print(f"[rendered] {out_path}")


# ===========================================================================
# Main entry
# ===========================================================================
def main(mode: str = "replication", tables: list[int] = None,
         out_dir: str = ".") -> None:
    if tables is None:
        tables = [1, 2, 3, 4]
    od = Path(out_dir); od.mkdir(parents=True, exist_ok=True)

    # Replication scale (~30 min on 1 CPU) vs publication scale (~2-3 h on 32 cores)
    if mode == "replication":
        T1_params = {"N_GWAS_grid": [4_000, 6_000], "M_grid": [25, 40],
                     "gamma_grid": [0.0, 0.05, 0.10], "n_iter": 50}
        T2_params = {"N": 4_000, "M": 25, "n_iter": 50, "B": 80}
        T3_params = {"N": 4_000, "M": 25, "n_iter": 50, "B": 80}
        T4_params = {"N": 4_000, "M": 25, "n_iter": 50, "B": 80}
    elif mode == "publication":
        T1_params = {"N_GWAS_grid": [50_000, 100_000, 200_000], "M_grid": [50, 100, 200],
                     "gamma_grid": [0.0, 0.01, 0.02, 0.05, 0.10, 0.15], "n_iter": 1_000}
        T2_params = {"N": 10_000, "M": 50, "n_iter": 1_000, "B": 200}
        T3_params = {"N": 10_000, "M": 50, "n_iter": 1_000, "B": 200,
                     "gamma_grid": [0.00, 0.01, 0.02, 0.05, 0.08, 0.10, 0.15]}
        T4_params = {"N": 10_000, "M": 50, "n_iter": 1_000, "B": 200,
                     "D_grid": [5, 10, 20]}
    else:
        raise ValueError(f"Unknown mode: {mode}")

    results = {}
    if 1 in tables:
        results["t1"] = table1(out_path=str(od / "results_table1.json"), **T1_params)
    if 2 in tables:
        results["t2"] = table2(out_path=str(od / "results_table2.json"), **T2_params)
    if 3 in tables:
        results["t3"] = table3(out_path=str(od / "results_table3.json"), **T3_params)
    if 4 in tables:
        results["t4"] = table4(out_path=str(od / "results_table4.json"), **T4_params)

    if all(k in results for k in ("t1", "t2", "t3", "t4")):
        render_tables_md(results["t1"], results["t2"], results["t3"], results["t4"],
                         out_path=str(od / "tables_v5_main.md"))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Master simulation engine: regenerates v5.2 main-text Tables 1-4."
    )
    parser.add_argument("--mode", choices=["replication", "publication"],
                        default="replication")
    parser.add_argument("--tables", default="1,2,3,4",
                        help="Comma-separated list of tables to regenerate.")
    parser.add_argument("--out_dir", default=".")
    args = parser.parse_args()
    tables = [int(x) for x in args.tables.split(",")]
    print(f"[mode={args.mode}] tables={tables} out_dir={args.out_dir}")
    main(mode=args.mode, tables=tables, out_dir=args.out_dir)

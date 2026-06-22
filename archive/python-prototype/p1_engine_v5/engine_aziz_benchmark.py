"""
engine_aziz_benchmark.py — Head-to-head simulation: cCWR vs per-component pooled MR.

Addresses the corresponding-author directive (Aziz) and cross-exam Q1 of the
v5.1 adversarial review: quantify cCWR's advantage over the natural baseline
that an applied epidemiologist would use today.

4 scenarios:
  A_null:                   alpha_X = (0, 0, 0)        — pure null, calibration check
  B_valid_IV:               alpha_X = (-0.4, -0.4, -0.4) — concordant causal effect
  C_pleiotropy:             alpha_X = (-0.4, -0.4, -0.4), gamma_1 = 0.05 — InSIDE violated
  D_hierarchy_discordant:   alpha_X = (+0.4,  0.0, -0.4) — sign reversal across components

For each scenario × n_iter Monte Carlo replicates:
  1. Simulate cohort.
  2. Fit cCWR DS-CWR via multiplier bootstrap (B_bootstrap).
  3. Fit per-component IVW, MR-Egger, weighted median, MR-PRESSO.
  4. Pool via inverse-variance, Bonferroni, Fisher.
  5. Record:
       point estimates (cCWR's log-DSCWR vs each pooled estimator)
       SEs / CIs
       reject decisions (alpha = 0.05)
       coverage of the true marginal causal direction

Output: results_aziz_benchmark.json — one entry per scenario, each containing
the empirical type-I error / power / coverage of every estimator.

This becomes Table 5 of v5.1 §6.4 (new section).
"""
from __future__ import annotations
import json, time
import numpy as np
from pathlib import Path
from dataclasses import asdict, replace

from .config import P1Config
from .dgp import simulate_cohort
from .kernel import precompute_kernel_matrix
from .multiplier_bootstrap import run_multiplier_bootstrap
from .benchmark_mr import run_per_component_benchmark
from scipy.stats import norm


def run_one_iter(cfg: P1Config, gamma_1: float, iter_idx: int,
                 D_strata: int = 10, B_bootstrap: int = 200) -> dict:
    """
    One Monte Carlo iteration: simulate cohort, run cCWR + per-component MR.

    Returns dict with all estimator results for this iteration.
    """
    seed = cfg.seed * 10_000 + iter_idx
    rng = np.random.default_rng(seed)
    cfg_iter = replace(cfg, gamma_direct=(gamma_1, 0.0, 0.0))
    coh = simulate_cohort(cfg_iter, rng)

    # cCWR via multiplier bootstrap
    H = precompute_kernel_matrix(coh["T"], coh["D"])
    G = coh["G"];  X = coh["X"]
    GtG = (G * G).sum(axis=0)
    bX  = (G.T @ X) / np.maximum(GtG, 1e-12)
    seX = np.sqrt(np.var(X, ddof=1) / np.maximum(GtG, 1e-12))
    Z   = rng.normal(size=(cfg.N_outcome, 3))
    boot = run_multiplier_bootstrap(H, G, X, Z, bX, seX,
                                    D_strata=D_strata, B=B_bootstrap, seed=seed)
    delta_GLS = float(boot["delta_GLS"])
    se_GLS    = float(boot["se_delta_GLS"])
    z_cwr     = delta_GLS / max(se_GLS, 1e-12)
    p_cwr     = 2.0 * (1.0 - norm.cdf(abs(z_cwr)))

    # Per-component pooled MR
    benchmark = run_per_component_benchmark(G, X, coh["T"], coh["D"], rng)

    out = {
        "cCWR": {"log_DSCWR": delta_GLS, "se": se_GLS, "p": p_cwr,
                 "DSCWR": float(np.exp(delta_GLS))},
        **benchmark,
    }
    return out


def run_scenario(name: str, cfg: P1Config, gamma_1: float,
                 n_iter: int, B_bootstrap: int = 200,
                 D_strata: int = 10) -> dict:
    """Aggregate n_iter MC iterations for one scenario."""
    print(f"[{time.strftime('%H:%M:%S')}] Scenario {name}: alpha_X={cfg.alpha_X}, gamma_1={gamma_1}, n_iter={n_iter}")
    results = []
    t0 = time.time()
    n_failures = 0
    for k in range(n_iter):
        try:
            r = run_one_iter(cfg, gamma_1, k, D_strata, B_bootstrap)
            results.append(r)
        except Exception as e:
            n_failures += 1
            if n_failures <= 3:
                print(f"   iter {k} failed: {type(e).__name__}: {e}")
    runtime_s = time.time() - t0

    # Aggregate
    def _summary(estimator: str, key: str = "p", thresh: float = 0.05) -> dict:
        ps = []
        thetas = []
        ses = []
        for r in results:
            if estimator not in r: continue
            entry = r[estimator]
            if key in entry and entry[key] is not None:
                ps.append(entry[key])
            for vk in ("log_DSCWR", "theta"):
                if vk in entry:
                    thetas.append(entry[vk]); break
            if "se" in entry:
                ses.append(entry["se"])
        ps     = np.array(ps,     dtype=float)
        thetas = np.array(thetas, dtype=float)
        ses    = np.array(ses,    dtype=float)
        # Trimmed-mean SE (drop top 10%) — robust against heavy-tailed bootstrap
        ses_finite = ses[np.isfinite(ses)]
        if ses_finite.size > 5:
            cutoff = np.quantile(ses_finite, 0.90)
            trimmed_ses = ses_finite[ses_finite <= cutoff]
        else:
            trimmed_ses = ses_finite
        return {
            "n_valid":          int(np.sum(np.isfinite(ps))),
            "rejection_rate":   float(np.mean(np.isfinite(ps) & (ps < thresh))),
            "mean_estimate":    float(np.mean(thetas)) if thetas.size else float("nan"),
            "median_estimate":  float(np.median(thetas)) if thetas.size else float("nan"),
            "mean_se":          float(np.mean(ses)) if ses.size else float("nan"),
            "median_se":        float(np.median(ses_finite)) if ses_finite.size else float("nan"),
            "trimmed_mean_se":  float(np.mean(trimmed_ses)) if trimmed_ses.size else float("nan"),
        }

    # The list of estimators to report
    estimator_keys = ["cCWR",
                      "IVW_pool_IV", "IVW_pool_Bonf", "IVW_pool_Fisher",
                      "Egger_pool_IV", "Egger_pool_Bonf", "Egger_pool_Fisher",
                      "WM_pool_IV", "WM_pool_Bonf", "WM_pool_Fisher",
                      "PRESSO_p1_global"]
    summary = {}
    for est in estimator_keys:
        if est == "PRESSO_p1_global":
            # use the p_global field
            ps = np.array([r["PRESSO_p1_global"]["p_global"]
                           for r in results if "PRESSO_p1_global" in r], dtype=float)
            summary[est] = {
                "n_valid":        int(np.sum(np.isfinite(ps))),
                "rejection_rate": float(np.mean(np.isfinite(ps) & (ps < 0.05))),
                "mean_estimate":  float("nan"),  # PRESSO global is a test, not a point estimate
                "median_estimate": float("nan"),
                "mean_se":        float("nan"),
            }
        elif est in ("IVW_pool_Fisher", "Egger_pool_Fisher", "WM_pool_Fisher"):
            # Fisher pool produces only a p-value (no point estimate)
            ps = np.array([r[est]["p"] for r in results if est in r], dtype=float)
            summary[est] = {
                "n_valid":        int(np.sum(np.isfinite(ps))),
                "rejection_rate": float(np.mean(np.isfinite(ps) & (ps < 0.05))),
                "mean_estimate":  float("nan"),
                "median_estimate": float("nan"),
                "mean_se":        float("nan"),
            }
        else:
            summary[est] = _summary(est, key="p")

    return {
        "config":     {"alpha_X": cfg.alpha_X, "alpha_S": cfg.alpha_S,
                       "gamma_1": gamma_1, "theta_F": cfg.theta_F,
                       "N_outcome": cfg.N_outcome, "M_snps": cfg.M_snps,
                       "D_strata": D_strata, "B_bootstrap": B_bootstrap,
                       "n_iter": n_iter},
        "n_failures": n_failures,
        "runtime_s":  round(runtime_s, 2),
        "estimators": summary,
    }


def main(out_path: str = "results_aziz_benchmark.json",
         n_iter: int = 100, N: int = 6_000, M: int = 40,
         B_bootstrap: int = 200, D_strata: int = 10,
         seed: int = 20_260_506) -> dict:
    """
    Run the 4-scenario head-to-head benchmark.
    """
    base = P1Config(N_outcome=N, M_snps=M, theta_F=0.8, alpha_S=0.4,
                    n_iter=n_iter, B_bootstrap=B_bootstrap, seed=seed)
    scenarios = [
        ("A_null",                 replace(base, alpha_X=(0.0, 0.0, 0.0)),       0.0),
        ("B_valid_IV",             replace(base, alpha_X=(-0.4, -0.4, -0.4)),    0.0),
        ("C_pleiotropy",           replace(base, alpha_X=(-0.4, -0.4, -0.4)),    0.05),
        ("D_hierarchy_discordant", replace(base, alpha_X=(+0.4,  0.0, -0.4)),    0.0),
    ]
    out = {}
    for name, cfg, g1 in scenarios:
        out[name] = run_scenario(name, cfg, g1, n_iter, B_bootstrap, D_strata)
        s = out[name]["estimators"]
        print(f"   {name:30s}  cCWR rej={s['cCWR']['rejection_rate']:.3f}  "
              f"IVW_IV rej={s['IVW_pool_IV']['rejection_rate']:.3f}  "
              f"IVW_Bonf rej={s['IVW_pool_Bonf']['rejection_rate']:.3f}  "
              f"WM_IV rej={s['WM_pool_IV']['rejection_rate']:.3f}")

    Path(out_path).write_text(json.dumps(out, indent=2, default=str))
    print(f"[saved] {out_path}")
    return out


if __name__ == "__main__":
    import argparse
    from .params import get_params, add_mode_arg
    parser = argparse.ArgumentParser(
        description="Aziz benchmark: head-to-head cCWR vs per-component pooled MR."
    )
    add_mode_arg(parser)
    parser.add_argument("--out", default="results_aziz_benchmark.json")
    parser.add_argument("--n_iter", type=int, default=None)
    parser.add_argument("--N", type=int, default=None)
    parser.add_argument("--M", type=int, default=None)
    parser.add_argument("--B", type=int, default=None,
                        help="Override B_bootstrap.")
    parser.add_argument("--D", type=int, default=10,
                        help="D strata; reduce to 5 at small N.")
    args = parser.parse_args()
    pp = get_params(args.mode)
    n_iter_eff   = args.n_iter if args.n_iter is not None else (50 if args.mode == "replication" else 1_000)
    N_eff        = args.N      if args.N      is not None else pp["N_outcome"]
    M_eff        = args.M      if args.M      is not None else pp["M_snps"]
    B_eff        = args.B      if args.B      is not None else pp["C_C2_B"]
    print(f"[mode={args.mode}] N={N_eff}, M={M_eff}, n_iter={n_iter_eff}, B={B_eff}, D={args.D}")
    main(out_path=args.out, n_iter=n_iter_eff, N=N_eff, M=M_eff,
         B_bootstrap=B_eff, D_strata=args.D, seed=pp["seed"])

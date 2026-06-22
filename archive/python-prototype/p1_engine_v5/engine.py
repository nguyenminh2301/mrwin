"""
engine.py — Run Package A R2 (InSIDE robustness for SDPD).

For each frailty level theta_F ∈ {0, 0.4, 0.8}, simulate n_iter cohorts
under gamma_1 = 0 (the SDPD null) and:

  1. Extract per-SNP beta_X_hat by OLS of X on each SNP.
  2. Extract per-SNP beta_T1_hat by Cox (multiplicative variant) and
     by Aalen (additive variant) for priority-1 mortality.
  3. Fit random-effects MR-Egger of beta_T1_hat on beta_X_hat.
  4. Test H_0: theta_0 = 0 at alpha = 0.05 (two-sided z).
  5. Record rejection (= empirical Type-I error if InSIDE holds; > 0.05
     under shared frailty signals InSIDE violation).

Writes results.json with one record per (scenario, scale).
"""
from __future__ import annotations
import json
import time
import numpy as np
from dataclasses import asdict
from pathlib import Path
from scipy.stats import norm

from .config import P1Config, insider_scenarios
from .dgp import simulate_cohort
from .inference import cox_per_snp, aalen_per_snp, mr_egger_re


def _gwas_X(G: np.ndarray, X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Per-SNP OLS of X on G_m. G is column-standardised."""
    N, M = G.shape
    GtG_diag = (G * G).sum(axis=0)
    beta = (G.T @ X) / np.maximum(GtG_diag, 1e-12)
    resid_var = float(np.var(X, ddof=1))
    se = np.sqrt(resid_var / np.maximum(GtG_diag, 1e-12))
    return beta, se


def run_one_iteration(cfg: P1Config, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    coh = simulate_cohort(cfg, rng)

    G  = coh["G"]
    X  = coh["X"]
    T1 = coh["T"][:, 0]
    D1 = coh["D"][:, 0]

    # Optional ascertainment: keep only individuals whose observed time exceeds
    # the median observed time in the cohort (a stylised prevalent-sample
    # selection). Under shared frailty this selection induces InSIDE violation.
    if cfg.ascertain_on_survival:
        threshold = float(np.median(T1))
        keep = T1 >= threshold
        G  = G[keep]
        X  = X[keep]
        T1 = T1[keep]
        D1 = D1[keep]

    # --- Per-SNP beta_X (instrument-strength side of MR-Egger)
    bX, seX = _gwas_X(G, X)

    # --- Cox-PH variant for priority-1 mortality
    bY_cox, seY_cox = cox_per_snp(G, T1, D1)
    # --- Aalen additive variant for priority-1 mortality
    bY_aal, seY_aal = aalen_per_snp(G, T1, D1)

    eg_cox = mr_egger_re(bX, bY_cox, seY_cox)
    eg_aal = mr_egger_re(bX, bY_aal, seY_aal)

    return {
        "intercept_cox":     eg_cox["intercept"],
        "intercept_se_cox":  eg_cox["intercept_se"],
        "slope_cox":         eg_cox["slope"],
        "slope_se_cox":      eg_cox["slope_se"],
        "Q_cox":             eg_cox["Q"],
        "intercept_aalen":   eg_aal["intercept"],
        "intercept_se_aalen":eg_aal["intercept_se"],
        "slope_aalen":       eg_aal["slope"],
        "slope_se_aalen":    eg_aal["slope_se"],
        "Q_aalen":           eg_aal["Q"],
    }


def run_scenario(name: str, cfg: P1Config) -> dict:
    """Run n_iter iterations under cfg and aggregate."""
    iter_records = []
    z_crit = norm.ppf(1 - 0.025)

    t_start = time.time()
    for k in range(cfg.n_iter):
        # Per-iteration seed: deterministic, audited.
        seed_k = cfg.seed * 10_000 + k
        try:
            rec = run_one_iteration(cfg, seed_k)
        except Exception as exc:
            iter_records.append({"iter": k, "failed": True, "error": str(exc)})
            continue
        iter_records.append({"iter": k, "failed": False, **rec})
    runtime = time.time() - t_start

    # Aggregate: empirical Type-I error of intercept = 0 under each scale.
    successful = [r for r in iter_records if not r["failed"]]
    n_ok = len(successful)
    n_fail = len(iter_records) - n_ok

    def _agg(scale: str) -> dict:
        ic    = np.array([r[f"intercept_{scale}"]    for r in successful])
        ic_se = np.array([r[f"intercept_se_{scale}"] for r in successful])
        z     = ic / np.maximum(ic_se, 1e-12)
        reject = (np.abs(z) > z_crit).astype(float)
        return {
            "n_successful":         n_ok,
            "mean_intercept":       float(ic.mean()),
            "sd_intercept":         float(ic.std(ddof=1)) if n_ok > 1 else float("nan"),
            "mean_intercept_se":    float(ic_se.mean()),
            "type_I_error":         float(reject.mean()),
            "type_I_error_mc_se":   float(np.sqrt(reject.mean() * (1 - reject.mean()) / max(n_ok, 1))),
        }

    return {
        "scenario":  name,
        "config":    asdict(cfg),
        "runtime_s": round(runtime, 2),
        "n_iter":    cfg.n_iter,
        "n_failures": n_fail,
        "cox":       _agg("cox"),
        "aalen":     _agg("aalen"),
    }


def _run_scenario_worker(args):
    name, cfg = args
    return name, run_scenario(name, cfg)


def main(out_path: str = "results_packageA_R2.json", n_proc: int = 1,
         N_outcome: int | None = None, M_snps: int | None = None,
         n_iter: int | None = None, seed: int | None = None) -> dict:
    scenarios = insider_scenarios(N_outcome=N_outcome, M_snps=M_snps,
                                  n_iter=n_iter, seed=seed)
    out = {}

    if n_proc <= 1:
        for name, cfg in scenarios.items():
            print(f"[{time.strftime('%H:%M:%S')}] Running scenario: {name} (theta_F={cfg.theta_F}, ascertain={cfg.ascertain_on_survival}, n_iter={cfg.n_iter})")
            result = run_scenario(name, cfg)
            out[name] = result
            print(
                f"   {name:40s}  Cox Type-I={result['cox']['type_I_error']:.3f}  "
                f"Aalen Type-I={result['aalen']['type_I_error']:.3f}  "
                f"runtime={result['runtime_s']:.1f}s  failures={result['n_failures']}"
            )
    else:
        from multiprocessing import Pool
        items = list(scenarios.items())
        print(f"[{time.strftime('%H:%M:%S')}] Launching {len(items)} scenarios across {n_proc} processes...")
        with Pool(n_proc) as pool:
            for name, result in pool.imap_unordered(_run_scenario_worker, items):
                out[name] = result
                print(
                    f"   {name:40s}  Cox Type-I={result['cox']['type_I_error']:.3f}  "
                    f"Aalen Type-I={result['aalen']['type_I_error']:.3f}  "
                    f"runtime={result['runtime_s']:.1f}s"
                )

    Path(out_path).write_text(json.dumps(out, indent=2))
    print(f"[done] Wrote {out_path}")
    return out


if __name__ == "__main__":
    import argparse
    from .params import get_params, add_mode_arg
    parser = argparse.ArgumentParser(
        description="Package A R2: InSIDE robustness for SDPD MR-Egger intercept."
    )
    add_mode_arg(parser)
    parser.add_argument("--out", default="results_packageA_R2.json")
    parser.add_argument("--n_iter", type=int, default=None,
                        help="Manual override of n_iter (otherwise comes from --mode).")
    parser.add_argument("--N", type=int, default=None,
                        help="Manual override of N_outcome (otherwise comes from --mode).")
    parser.add_argument("--n_proc", type=int, default=1,
                        help="Number of parallel processes for scenarios.")
    args = parser.parse_args()
    p = get_params(args.mode)
    N_outcome_eff = args.N      if args.N      is not None else p["N_outcome"]
    n_iter_eff    = args.n_iter if args.n_iter is not None else p["A_R2_n_iter"]
    print(f"[mode={args.mode}] N={N_outcome_eff}, M={p['M_snps']}, n_iter={n_iter_eff}, n_proc={args.n_proc}")
    main(out_path=args.out, n_proc=args.n_proc,
         N_outcome=N_outcome_eff, M_snps=p["M_snps"],
         n_iter=n_iter_eff, seed=p["seed"])

"""
engine_packageB.py — Run R3 (scale reconciliation for SDPD power) and
R4 (bootstrap fourth-moment diagnostic for the bivariate Delta method).

R3 design
---------
We replicate v5 Table 1 across two scales of the gamma_1 injection plus a
gap column:

  (a) Multiplicative injection: gamma_1 enters via lp_1 += gamma_1 * S_true
      in the Weibull-PH DGP (DGP scale = log-hazard). Fit Cox-MR-Egger.
  (b) Additive injection: gamma_1 enters as a perturbation on the additive
      Aalen rate structure. Fit Aalen-MR-Egger.
  (c) Cross variant (operational SDPD as actually run by v5):
      Multiplicative injection on the Weibull-PH DGP, fit Aalen-MR-Egger.
      The Type-I+power gap between (a) and (c) quantifies the
      multiplicative-additive scale gap.

For each scale we sweep gamma_1 across {0.00, 0.02, 0.05, 0.10, 0.15} and
report empirical power (proportion of MR-Egger intercept p < 0.05).

R4 design
---------
For each of three scenarios reflecting the v5 Scenario A/B/D taxonomy
(B = valid IV; D = weak IV; A = null), simulate one cohort, run the full
multiplier bootstrap with B = 500, and record per-contrast kurtosis and
skewness of ln theta*^(b, d, d-1) and Delta X*^(b, d, d-1).
"""
from __future__ import annotations
import json
import time
import numpy as np
from dataclasses import asdict, replace
from pathlib import Path
from scipy.stats import norm

from .config import P1Config
from .dgp import simulate_cohort
from .inference import cox_per_snp, aalen_per_snp, mr_egger_re
from .kernel import precompute_kernel_matrix
from .multiplier_bootstrap import run_multiplier_bootstrap


# ---------------------------------------------------------------------------
# R3: Scale reconciliation
# ---------------------------------------------------------------------------
def _gwas_X(G: np.ndarray, X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    GtG = (G * G).sum(axis=0)
    beta = (G.T @ X) / np.maximum(GtG, 1e-12)
    se   = np.sqrt(np.var(X, ddof=1) / np.maximum(GtG, 1e-12))
    return beta, se


def run_R3_one_scenario(
    name: str, cfg: P1Config, gamma_1_values: list[float], n_iter: int
) -> dict:
    """
    Run R3 for one (theta_F, ascertainment) configuration over a grid of
    gamma_1 values. Returns per-gamma_1 empirical Type-I/power on both
    Cox- and Aalen-MR-Egger intercept tests.
    """
    z_crit = norm.ppf(1 - 0.025)
    out = {"scenario": name, "config_base": asdict(cfg),
           "gamma_grid": list(gamma_1_values), "n_iter": n_iter, "rows": []}
    t_start = time.time()
    for gamma in gamma_1_values:
        cfg_g = replace(cfg, gamma_direct=(gamma, 0.0, 0.0))
        rej_cox = 0; rej_aal = 0
        ok = 0
        for k in range(n_iter):
            seed_k = cfg_g.seed * 10_000 + k + int(gamma * 100_000)
            rng = np.random.default_rng(seed_k)
            try:
                coh = simulate_cohort(cfg_g, rng)
            except Exception:
                continue
            G = coh["G"]; X = coh["X"]
            T1 = coh["T"][:, 0]; D1 = coh["D"][:, 0]

            bX, _    = _gwas_X(G, X)
            bY_cox, sY_cox = cox_per_snp(G, T1, D1)
            bY_aal, sY_aal = aalen_per_snp(G, T1, D1)
            eg_cox = mr_egger_re(bX, bY_cox, sY_cox)
            eg_aal = mr_egger_re(bX, bY_aal, sY_aal)
            ok += 1
            if abs(eg_cox["intercept"] / max(eg_cox["intercept_se"], 1e-12)) > z_crit:
                rej_cox += 1
            if abs(eg_aal["intercept"] / max(eg_aal["intercept_se"], 1e-12)) > z_crit:
                rej_aal += 1
        rate_cox = rej_cox / max(ok, 1)
        rate_aal = rej_aal / max(ok, 1)
        out["rows"].append({
            "gamma_1":    float(gamma),
            "n_ok":       int(ok),
            "power_cox":  float(rate_cox),
            "power_aal":  float(rate_aal),
            "scale_gap":  float(rate_cox - rate_aal),
            "mc_se_cox":  float(np.sqrt(rate_cox * (1 - rate_cox) / max(ok, 1))),
            "mc_se_aal":  float(np.sqrt(rate_aal * (1 - rate_aal) / max(ok, 1))),
        })
        print(f"   [{name}] gamma_1={gamma:.3f}  Cox-power={rate_cox:.3f}  Aalen-power={rate_aal:.3f}  gap={rate_cox-rate_aal:+.3f}")
    out["runtime_s"] = round(time.time() - t_start, 1)
    return out


def run_R3(out_path: str, n_iter: int = 300, N: int = 6_000, M: int = 40,
           gamma_grid: list[float] | None = None) -> dict:
    """
    Run R3 in the v5 baseline regime (theta_F = 0.8, unascertained,
    matching Section 6.1's DGP and Section 7's data flow).
    """
    if gamma_grid is None:
        gamma_grid = [0.00, 0.02, 0.05, 0.10, 0.15]
    base = P1Config(N_outcome=N, M_snps=M, theta_F=0.8, n_iter=n_iter)
    print(f"[{time.strftime('%H:%M:%S')}] R3 sweep: gamma_1 in {gamma_grid}, n_iter={n_iter}, N={base.N_outcome}, M={base.M_snps}")
    out = run_R3_one_scenario("R3_v5_baseline", base, gamma_grid, n_iter)
    Path(out_path).write_text(json.dumps(out, indent=2))
    print(f"[done] Wrote {out_path}")
    return out


# ---------------------------------------------------------------------------
# R4: Bootstrap fourth-moment diagnostic
# ---------------------------------------------------------------------------
def run_R4_one_scenario(
    name: str, cfg: P1Config, B: int = 500, n_pcs: int = 3
) -> dict:
    """
    Single cohort, one full bootstrap; record kurtosis/skewness diagnostics.
    """
    rng = np.random.default_rng(cfg.seed)
    print(f"[{time.strftime('%H:%M:%S')}] R4 [{name}] simulate cohort N={cfg.N_outcome}, M={cfg.M_snps}, alpha_S={cfg.alpha_S}, alpha_X={cfg.alpha_X}, gamma_1={cfg.gamma_direct[0]}")
    t = time.time()
    coh = simulate_cohort(cfg, rng)
    print(f"   cohort sim: {time.time()-t:.1f}s")

    t = time.time()
    H = precompute_kernel_matrix(coh["T"], coh["D"])
    print(f"   kernel precompute: {time.time()-t:.1f}s")

    # Estimate beta_hat from the cohort itself (one-sample) by per-SNP OLS of X on G.
    G = coh["G"]; X = coh["X"]
    bX, seX = _gwas_X(G, X)
    Z = rng.normal(size=(cfg.N_outcome, n_pcs))   # surrogate PCs

    t = time.time()
    res = run_multiplier_bootstrap(
        H, G, X, Z, bX, seX, D_strata=10, B=B, seed=cfg.seed,
        use_iptw=False,    # IPTW disabled; PCs are noise here
    )
    print(f"   bootstrap (B={B}): {time.time()-t:.1f}s, max kurt(ln theta*) = {res['max_kurt_log_theta']:.2f}, max kurt(Delta X*) = {res['max_kurt_delta_X']:.2f}")

    return {
        "scenario":     name,
        "config":       asdict(cfg),
        "B":            B,
        **{k: v for k, v in res.items() if k not in {"Sigma_ISG", "Sigma_LW"}},
    }


def run_R4(out_path: str, B: int = 500, N: int = 6_000, M: int = 40, seed: int = 20_260_506) -> dict:
    """
    Three R4 scenarios:
      - A: null (alpha_X = 0)
      - B: valid IV (alpha_X = -0.4)
      - D: weak IV (alpha_S = 0.10, alpha_X = -0.4)
    """
    scenarios = {
        "R4_A_null":         P1Config(N_outcome=N, M_snps=M, theta_F=0.8,
                                      alpha_S=0.4,  alpha_X=(0.0, 0.0, 0.0),
                                      gamma_direct=(0.0, 0.0, 0.0), seed=seed),
        "R4_B_valid_IV":     P1Config(N_outcome=N, M_snps=M, theta_F=0.8,
                                      alpha_S=0.4,  alpha_X=(-0.4, -0.4, -0.4),
                                      gamma_direct=(0.0, 0.0, 0.0), seed=seed),
        "R4_D_weak_IV":      P1Config(N_outcome=N, M_snps=M, theta_F=0.8,
                                      alpha_S=0.10, alpha_X=(-0.4, -0.4, -0.4),
                                      gamma_direct=(0.0, 0.0, 0.0), seed=seed),
    }
    out = {}
    for name, cfg in scenarios.items():
        out[name] = run_R4_one_scenario(name, cfg, B=B)
    Path(out_path).write_text(json.dumps(out, indent=2, default=str))
    print(f"[done] Wrote {out_path}")
    return out


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import argparse
    from .params import get_params, add_mode_arg
    p = argparse.ArgumentParser(
        description="Package B: R3 (scale reconciliation) + R4 (bootstrap moments)."
    )
    add_mode_arg(p)
    p.add_argument("--task", choices=["R3", "R4", "both"], default="both")
    p.add_argument("--R3_n_iter", type=int, default=None,
                   help="Manual override (otherwise from --mode).")
    p.add_argument("--R4_B", type=int, default=None,
                   help="Manual override (otherwise from --mode).")
    p.add_argument("--out_dir", default=".")
    args = p.parse_args()
    pp = get_params(args.mode)
    R3_n_iter = args.R3_n_iter if args.R3_n_iter is not None else pp["B_R3_n_iter"]
    R4_B      = args.R4_B      if args.R4_B      is not None else pp["B_R4_B"]
    print(f"[mode={args.mode}] N={pp['N_outcome']}, M={pp['M_snps']}, R3_n_iter={R3_n_iter}, R4_B={R4_B}")
    if args.task in ("R3", "both"):
        run_R3(out_path=f"{args.out_dir}/results_packageB_R3.json",
               n_iter=R3_n_iter, N=pp["N_outcome"], M=pp["M_snps"],
               gamma_grid=pp["B_R3_gamma_grid"])
    if args.task in ("R4", "both"):
        run_R4(out_path=f"{args.out_dir}/results_packageB_R4.json",
               B=R4_B, N=pp["N_outcome"], M=pp["M_snps"], seed=pp["seed"])

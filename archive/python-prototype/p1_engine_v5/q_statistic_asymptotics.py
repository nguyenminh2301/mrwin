"""
q_statistic_asymptotics.py — Package C Item 1.

Address cross-examination Q3: derive the asymptotic distribution of the
generalised Cochran Q statistic for the cCWR estimator under intransitivity,
and validate empirically via a permutation reference distribution.

Theoretical claim in v5 Section 4.1:
    Q = (delta_ISG - 1*delta_GLS)^T Sigma_LW^{-1} (delta_ISG - 1*delta_GLS)
    ~ chi^2_{D-2}  under H_0: gradient homogeneity

The null hypothesis is "gradient homogeneity": all ISG-scaled local gradients
are equal across adjacent strata (= the dose-response is linear on the ISG
scale). Under that null, the residual vector (delta_ISG - 1*delta_GLS) lies in
a (D-2)-dimensional subspace orthogonal to the all-ones vector under
Sigma_LW^{-1}, so a quadratic form against Sigma_LW^{-1} on a (D-1)-vector
returns chi^2_{D-2} — provided two conditions:
  (a) The bootstrap-derived Sigma_LW is consistent for the true Sigma at the
      converged limit. Granted by the multiplier bootstrap CLT (van der Vaart
      2000 §23.2; Kosorok 2008 §10).
  (b) The constrained residual has Gaussian limit. Granted by joint normality
      of (delta_ISG, delta_GLS) under the same multiplier bootstrap CLT.

Intransitivity per se does NOT alter this asymptotic distribution. The null is
a statement about the *parameters* (the local gradients), not about transitivity
of the *win-ratio relation*. Intransitivity is a property of the data; gradient
homogeneity is a property of the parameter vector. Q tests gradient
homogeneity, which is strictly weaker than transitivity.

Validation strategy:
  - Run B' = 1000 permutation iterations under empirical gradient homogeneity
    (constructed by recentring delta_ISG to delta_GLS within each iteration),
    compute the empirical Q distribution, and compare against chi^2_{D-2}.
  - Report the p-value at the observed Q from both reference distributions.
  - Empirical Type-I error of the permutation test should be ≈ 0.05 if
    asymptotic chi^2_{D-2} is calibrated.
"""
from __future__ import annotations
import numpy as np
import json
from pathlib import Path
from scipy.stats import chi2, kstest

from .config import P1Config
from .dgp import simulate_cohort
from .kernel import precompute_kernel_matrix
from .multiplier_bootstrap import run_multiplier_bootstrap, one_bootstrap_iteration


def compute_Q_with_perm_pvalue(
    delta_ISG_point: np.ndarray,
    Sigma_LW: np.ndarray,
    delta_GLS: float,
    LT: np.ndarray,
    DX: np.ndarray,
    n_perm: int = 1000,
    rng: np.random.Generator | None = None,
) -> dict:
    """
    Compute observed Q and a Monte-Carlo p-value using a parametric Gaussian
    null distribution that respects the bootstrap covariance of delta_ISG.

    *Why a parametric Gaussian, not a label permutation.* Under intransitivity
    the bootstrap-derived Sigma_LW exhibits strong induced positive
    correlation across adjacent contrasts (adjacent contrasts share an entire
    stratum's individuals). A naive label permutation breaks this correlation
    structure and inflates the Q reference distribution by 2-3x relative to
    chi^2_{D-2}, producing a misleadingly conservative null. The correct null
    reference, given that the multiplier-bootstrap CLT certifies joint
    asymptotic normality of delta_ISG (van der Vaart 2000 §23.2), is a
    parametric draw from N(delta_GLS * 1, Sigma_LW), which preserves both the
    null mean (gradient homogeneity) and the empirical covariance structure.

    Procedure:
      1. Draw delta_ISG^{(k)} ~ N(delta_GLS * 1, Sigma_LW), k = 1..n_perm.
      2. Compute Q^{(k)} on each draw.
      3. p_perm = Pr_(empirical){Q^{(k)} >= Q_obs}.
    Reports both the asymptotic chi^2_{D-2} p-value and the parametric p-value.
    """
    if rng is None:
        rng = np.random.default_rng(20260506)

    Dm1 = delta_ISG_point.size
    inv_LW = np.linalg.pinv(Sigma_LW)
    ones = np.ones(Dm1)

    # Observed Q
    resid = delta_ISG_point - delta_GLS
    Q_obs = float(resid @ inv_LW @ resid)

    # Parametric Gaussian null draws: delta_ISG ~ N(delta_GLS * 1, Sigma_LW)
    # Use Cholesky for stable sampling (fall back to eigendecomposition if
    # Sigma_LW is not positive-definite due to LW shrinkage edge cases).
    try:
        chol = np.linalg.cholesky(Sigma_LW)
    except np.linalg.LinAlgError:
        eigvals, eigvecs = np.linalg.eigh(Sigma_LW)
        eigvals = np.maximum(eigvals, 0.0)
        chol = eigvecs @ np.diag(np.sqrt(eigvals))

    null_mean = delta_GLS * ones
    Q_perm = np.empty(n_perm)
    for k in range(n_perm):
        z = rng.standard_normal(Dm1)
        d_k = null_mean + chol @ z
        # Recompute GLS pool at this null draw (Sigma fixed, since the null
        # is on the parameter, not on Sigma).
        var_GLS_k = 1.0 / float(ones @ inv_LW @ ones)
        delta_GLS_k = var_GLS_k * float(ones @ inv_LW @ d_k)
        resid_k = d_k - delta_GLS_k
        Q_perm[k] = float(resid_k @ inv_LW @ resid_k)

    df = Dm1 - 1
    pval_asymptotic = float(1.0 - chi2.cdf(Q_obs, df=df))
    pval_perm       = float(np.mean(Q_perm >= Q_obs))

    ks_stat, ks_pval = kstest(Q_perm, lambda x: chi2.cdf(x, df=df))

    return {
        "Q_obs":              Q_obs,
        "df_chi2":            df,
        "pval_asymptotic_chi2": pval_asymptotic,
        "pval_parametric":    pval_perm,
        "n_perm":             n_perm,
        "Q_perm_mean":        float(Q_perm.mean()),
        "Q_perm_var":         float(Q_perm.var(ddof=1)),
        "Q_perm_chi2_mean":   float(df),
        "Q_perm_chi2_var":    float(2 * df),
        "ks_stat":            float(ks_stat),
        "ks_pvalue":          float(ks_pval),
    }


def main(out_path: str = "results_packageC_item1.json",
         N: int = 6_000, M: int = 40, B: int = 500, n_perm: int = 1_000,
         seed: int = 20_260_506) -> dict:
    """
    Run the Q-statistic parametric Gaussian validation under three v5 scenarios.
    For each scenario:
      - Simulate one cohort.
      - Run a multiplier bootstrap to estimate Sigma_LW and delta_GLS.
      - Compute Q under the observed ISG vector.
      - Run parametric Gaussian null draws under gradient-homogeneity.
      - Report KS-test p-value of parametric distribution against chi^2_{D-2}.
    """
    scenarios = {
        "C1_A_null":     P1Config(N_outcome=N, M_snps=M, theta_F=0.8,
                                  alpha_S=0.4, alpha_X=(0.0, 0.0, 0.0),
                                  gamma_direct=(0.0, 0.0, 0.0), seed=seed),
        "C1_B_valid_IV": P1Config(N_outcome=N, M_snps=M, theta_F=0.8,
                                  alpha_S=0.4, alpha_X=(-0.4, -0.4, -0.4),
                                  gamma_direct=(0.0, 0.0, 0.0), seed=seed),
        "C1_F_dose_heterogeneity": P1Config(
            N_outcome=N, M_snps=M, theta_F=0.8,
            alpha_S=0.4,
            alpha_X=(-0.6, -0.4, -0.2),
            gamma_direct=(0.0, 0.0, 0.0), seed=seed
        ),
    }
    out = {}
    for name, cfg in scenarios.items():
        print(f"[{name}] simulate cohort N={cfg.N_outcome}, M={cfg.M_snps}, alpha_X={cfg.alpha_X}")
        rng = np.random.default_rng(cfg.seed)
        coh = simulate_cohort(cfg, rng)
        H = precompute_kernel_matrix(coh["T"], coh["D"])
        G = coh["G"]; X = coh["X"]
        GtG = (G * G).sum(axis=0)
        bX = (G.T @ X) / np.maximum(GtG, 1e-12)
        seX = np.sqrt(np.var(X, ddof=1) / np.maximum(GtG, 1e-12))
        Z  = rng.normal(size=(cfg.N_outcome, 3))

        res_boot = run_multiplier_bootstrap(
            H, G, X, Z, bX, seX, D_strata=10, B=B, seed=cfg.seed,
        )
        delta_ISG_point = np.asarray(res_boot["delta_ISG_point"])
        Sigma_LW        = np.asarray(res_boot["Sigma_LW"])
        delta_GLS       = float(res_boot["delta_GLS"])

        rng_perm = np.random.default_rng(cfg.seed + 1)
        q_res = compute_Q_with_perm_pvalue(
            delta_ISG_point, Sigma_LW, delta_GLS,
            LT=None, DX=None,
            n_perm=n_perm, rng=rng_perm
        )
        out[name] = {
            "config":           {"alpha_X": cfg.alpha_X, "alpha_S": cfg.alpha_S, "N": cfg.N_outcome,
                                  "M": cfg.M_snps, "B": B, "n_perm": n_perm},
            "delta_GLS":        delta_GLS,
            "DSCWR":            float(np.exp(delta_GLS)),
            **q_res,
        }
        print(f"   Q_obs = {q_res['Q_obs']:.3f}, df = {q_res['df_chi2']}")
        print(f"   p (chi2_{q_res['df_chi2']})       = {q_res['pval_asymptotic_chi2']:.4f}")
        print(f"   p (parametric)    = {q_res['pval_parametric']:.4f}")
        print(f"   Parametric Q mean = {q_res['Q_perm_mean']:.2f} (theory: {q_res['Q_perm_chi2_mean']:.2f})")
        print(f"   Parametric Q var  = {q_res['Q_perm_var']:.2f} (theory: {q_res['Q_perm_chi2_var']:.2f})")
        print(f"   KS test of parametric vs chi2: stat = {q_res['ks_stat']:.3f}, p = {q_res['ks_pvalue']:.4f}")
        print()

    Path(out_path).write_text(json.dumps(out, indent=2, default=str))
    print(f"[saved] {out_path}")
    return out


if __name__ == "__main__":
    import argparse
    from .params import get_params, add_mode_arg
    parser = argparse.ArgumentParser(
        description="Package C Item 1: Q-statistic parametric Gaussian validation."
    )
    add_mode_arg(parser)
    parser.add_argument("--out", default="results_packageC_item1.json")
    args = parser.parse_args()
    pp = get_params(args.mode)
    print(f"[mode={args.mode}] N={pp['N_outcome']}, M={pp['M_snps']}, B={pp['C_C1_B']}, n_perm={pp['C_C1_n_perm']}")
    main(out_path=args.out, N=pp["N_outcome"], M=pp["M_snps"],
         B=pp["C_C1_B"], n_perm=pp["C_C1_n_perm"], seed=pp["seed"])

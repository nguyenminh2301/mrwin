"""
benchmark_mr.py — Aziz benchmark: per-component pooled MR vs cCWR.

Addresses cross-exam Q1: quantify cCWR's empirical advantage over the natural
baseline an applied epidemiologist would use today: separately fit IVW /
MR-Egger / weighted-median / MR-PRESSO on each component outcome, then pool.

The 3 component outcomes:
  k = 1: time-to-priority-1 (mortality)
  k = 2: time-to-priority-2 (HF hospitalisation)
  k = 3: time-to-priority-3 (renal decline)

For each component k and each per-SNP exposure summary (b_X, se_X) we compute
the per-SNP outcome summary (b_Y_k, se_Y_k) on the Aalen additive scale (per
the v5 §5.2 mandate). Then we run four MR estimators per component:

  1. IVW (inverse-variance-weighted)
  2. MR-Egger (with intercept)
  3. Weighted median (50% breakdown)
  4. MR-PRESSO (outlier-robust; we implement the global test variant)

Then we pool across components by three pooling rules:
  - Bonferroni: report min(p_1, p_2, p_3) * 3, with point estimate from the
    component with the smallest p-value
  - Inverse-variance: precision-weighted average of the three component
    estimates
  - Fisher's combined p-value: -2 * sum(log(p_k)) ~ chi^2_{2K}

Each pooled estimator gives a (point estimate, SE/CI, p-value) triple. We
compare each against the cCWR DS-CWR on:
  - Type-I error at gamma_1 = 0 (Scenario A, B null configurations)
  - Power at gamma_1 > 0 (Scenario C: pleiotropy)
  - Power against true causal effect (Scenario B valid IV)
  - Hierarchy discordance: alpha_X = (+0.4, 0.0, -0.4) — mortality up, eGFR
    decline down. This is the most informative scenario for the comparison
    because the hierarchy-aware cCWR should agree with the *priority-1
    component* MR (mortality) while pooled rules will partially cancel.
"""
from __future__ import annotations
import numpy as np
from typing import Tuple, Dict
from scipy.stats import chi2, norm

from .inference import aalen_per_snp


# ---------------------------------------------------------------------------
# Per-component outcome summaries
# ---------------------------------------------------------------------------
def per_component_summaries(
    G: np.ndarray,
    T_components: np.ndarray,   # (N, K)
    D_components: np.ndarray,   # (N, K)
) -> Dict[int, Tuple[np.ndarray, np.ndarray]]:
    """
    For each priority k in 1..K, fit Aalen additive per-SNP summaries.
    Returns {k: (beta_Y_k, se_Y_k)} with k ∈ {1, 2, 3} = priorities.
    """
    K = T_components.shape[1]
    out = {}
    for k in range(K):
        b, se = aalen_per_snp(G, T_components[:, k], D_components[:, k])
        out[k + 1] = (b, se)   # 1-indexed priority
    return out


# ---------------------------------------------------------------------------
# IVW (inverse-variance-weighted)
# ---------------------------------------------------------------------------
def ivw(b_X: np.ndarray, b_Y: np.ndarray, se_Y: np.ndarray) -> Tuple[float, float, float]:
    """
    Standard fixed-effect IVW on summary data.
    Returns (theta_hat, se, two-sided p-value).
    Reference: Burgess, Butterworth, Thompson (2013).
    """
    w = b_X ** 2 / np.maximum(se_Y ** 2, 1e-20)
    theta = np.sum(w * b_Y / b_X) / np.sum(w)
    se    = 1.0 / np.sqrt(np.sum(w))
    z     = theta / max(se, 1e-12)
    p     = 2.0 * (1.0 - norm.cdf(abs(z)))
    return float(theta), float(se), float(p)


# ---------------------------------------------------------------------------
# MR-Egger (with intercept)
# ---------------------------------------------------------------------------
def mr_egger(b_X: np.ndarray, b_Y: np.ndarray, se_Y: np.ndarray) -> Tuple[float, float, float, float]:
    """
    MR-Egger fixed-effect regression: b_Y ~ alpha + beta * b_X, weighted by 1/se_Y^2.
    Returns (slope_hat, slope_se, slope_p, intercept_p).
    """
    w = 1.0 / np.maximum(se_Y ** 2, 1e-20)
    sw = np.sum(w)
    swx = np.sum(w * b_X);   swy = np.sum(w * b_Y)
    swxx = np.sum(w * b_X**2); swxy = np.sum(w * b_X * b_Y)
    den = sw * swxx - swx**2
    if den < 1e-20:
        return 0.0, np.inf, 1.0, 1.0
    slope = (sw * swxy - swx * swy) / den
    intercept = (swy - slope * swx) / sw
    # Residual variance + standard SE formulae for weighted least squares
    fit = intercept + slope * b_X
    rv = np.sum(w * (b_Y - fit) ** 2) / max(b_X.size - 2, 1)
    slope_se = np.sqrt(rv * sw / max(den, 1e-20))
    intercept_se = np.sqrt(rv * swxx / max(den, 1e-20))
    slope_p = 2.0 * (1.0 - norm.cdf(abs(slope / max(slope_se, 1e-12))))
    intercept_p = 2.0 * (1.0 - norm.cdf(abs(intercept / max(intercept_se, 1e-12))))
    return float(slope), float(slope_se), float(slope_p), float(intercept_p)


# ---------------------------------------------------------------------------
# Weighted median
# ---------------------------------------------------------------------------
def weighted_median_mr(b_X: np.ndarray, b_Y: np.ndarray, se_Y: np.ndarray,
                       n_boot: int = 200, rng: np.random.Generator | None = None) -> Tuple[float, float, float]:
    """
    Weighted median MR estimator (Bowden 2016).
    The point estimate is the weighted median of theta_m = b_Y_m / b_X_m, with
    weights w_m = b_X_m^2 / se_Y_m^2. The SE is from a parametric bootstrap.
    Returns (theta_hat, se, two-sided p).
    """
    if rng is None:
        rng = np.random.default_rng(0)
    M = b_X.size
    theta_per_snp = b_Y / b_X
    w = b_X ** 2 / np.maximum(se_Y ** 2, 1e-20)
    w = w / w.sum()
    order = np.argsort(theta_per_snp)
    cum = np.cumsum(w[order])
    median_idx = np.searchsorted(cum, 0.5)
    median_idx = min(median_idx, M - 1)
    theta_hat = float(theta_per_snp[order[median_idx]])

    # Parametric bootstrap for SE
    se_X = np.full_like(b_X, 1e-3)   # default; caller can pass actual se_X if known
    boots = np.empty(n_boot)
    for k in range(n_boot):
        bX_b = b_X + rng.normal(scale=se_X)
        bY_b = b_Y + rng.normal(scale=se_Y)
        theta_b = bY_b / np.where(np.abs(bX_b) < 1e-6, 1e-6, bX_b)
        w_b = bX_b ** 2 / np.maximum(se_Y ** 2, 1e-20)
        w_b /= w_b.sum()
        ord_b = np.argsort(theta_b)
        cum_b = np.cumsum(w_b[ord_b])
        idx_b = min(np.searchsorted(cum_b, 0.5), M - 1)
        boots[k] = theta_b[ord_b[idx_b]]
    se = float(np.std(boots, ddof=1))
    z  = theta_hat / max(se, 1e-12)
    p  = 2.0 * (1.0 - norm.cdf(abs(z)))
    return theta_hat, se, float(p)


# ---------------------------------------------------------------------------
# MR-PRESSO (global test only — sufficient for benchmarking)
# ---------------------------------------------------------------------------
def mr_presso_global(b_X: np.ndarray, b_Y: np.ndarray, se_Y: np.ndarray,
                      n_perm: int = 1000, rng: np.random.Generator | None = None) -> Tuple[float, float, float, float]:
    """
    MR-PRESSO global outlier test (Verbanck 2018).

    1. Fit IVW. Compute residual sum of squares RSS_obs.
    2. Permute pairings of (b_X, b_Y) under H_0 of no horizontal pleiotropy.
    3. Compute null distribution of RSS_perm.
    4. Global p-value = Pr(RSS_perm >= RSS_obs).

    Returns (ivw_theta, ivw_se, ivw_p, presso_global_p).
    Note: this tests for the *presence* of pleiotropy. If significant, the
    correct downstream action is outlier removal + refit, but for benchmarking
    we report the global test as a calibrated diagnostic.
    """
    if rng is None:
        rng = np.random.default_rng(0)
    theta_ivw, se_ivw, p_ivw = ivw(b_X, b_Y, se_Y)
    fitted = theta_ivw * b_X
    resid = (b_Y - fitted) / np.maximum(se_Y, 1e-20)
    rss_obs = float(np.sum(resid ** 2))

    rss_perm = np.empty(n_perm)
    for k in range(n_perm):
        perm = rng.permutation(b_X.size)
        bY_p = b_Y[perm]
        theta_p, _, _ = ivw(b_X, bY_p, se_Y)
        resid_p = (bY_p - theta_p * b_X) / np.maximum(se_Y, 1e-20)
        rss_perm[k] = np.sum(resid_p ** 2)
    p_global = float(np.mean(rss_perm >= rss_obs))
    return theta_ivw, se_ivw, p_ivw, p_global


# ---------------------------------------------------------------------------
# Pool a list of (theta, se, p) triples by various rules
# ---------------------------------------------------------------------------
def pool_inverse_variance(comp_results: list) -> Tuple[float, float, float]:
    """
    Inverse-variance pooled estimate across components.
    Each comp_result = (theta_k, se_k, p_k).
    """
    thetas = np.array([r[0] for r in comp_results])
    ses    = np.array([r[1] for r in comp_results])
    ws     = 1.0 / np.maximum(ses ** 2, 1e-20)
    pooled = float(np.sum(ws * thetas) / np.sum(ws))
    se_p   = float(1.0 / np.sqrt(np.sum(ws)))
    z      = pooled / max(se_p, 1e-12)
    p      = 2.0 * (1.0 - norm.cdf(abs(z)))
    return pooled, se_p, p


def pool_fisher(comp_results: list) -> float:
    """Fisher's combined-p method. Returns combined p-value (no point estimate)."""
    ps = np.array([r[2] for r in comp_results])
    ps = np.clip(ps, 1e-300, 1 - 1e-12)
    chi2_stat = -2.0 * np.sum(np.log(ps))
    df = 2 * len(ps)
    return float(1.0 - chi2.cdf(chi2_stat, df=df))


def pool_bonferroni(comp_results: list) -> Tuple[float, float, float]:
    """
    Bonferroni-adjusted: take the most-significant component, multiply its
    p-value by K. Returns the best component's (theta, se, K * p).
    """
    K = len(comp_results)
    ps = [r[2] for r in comp_results]
    best = int(np.argmin(ps))
    th, se, p = comp_results[best]
    return th, se, min(1.0, K * p)


# ---------------------------------------------------------------------------
# Top-level: run all benchmarks on one cohort
# ---------------------------------------------------------------------------
def run_per_component_benchmark(
    G: np.ndarray, X: np.ndarray,
    T_comp: np.ndarray, D_comp: np.ndarray,
    rng: np.random.Generator,
) -> dict:
    """
    Returns a dict with one entry per (estimator, pooling_rule) for each of
    the 3 components. Provides:
      'IVW_pool_IV'      = IVW per component, then inverse-variance pool
      'IVW_pool_Bonf'    = IVW per component, then Bonferroni
      'IVW_pool_Fisher'  = IVW per component, then Fisher's combined p
      'Egger_pool_IV'    = MR-Egger per component, then inverse-variance pool
      'WM_pool_IV'       = weighted-median per component, then inverse-variance pool
      'PRESSO_p1'        = MR-PRESSO global test on priority-1 component only
      'cCWR_p1_only'     = (placeholder for cCWR; computed externally)
    """
    # Per-SNP exposure summary
    GtG = (G * G).sum(axis=0)
    b_X = (G.T @ X) / np.maximum(GtG, 1e-12)
    se_X = np.sqrt(np.var(X, ddof=1) / np.maximum(GtG, 1e-12))

    # Per-SNP per-component outcome summaries
    comp_summaries = per_component_summaries(G, T_comp, D_comp)

    # Per-component MR results
    component_ivw = {}
    component_eg  = {}
    component_wm  = {}
    component_presso = {}
    for k in (1, 2, 3):
        bY, sY = comp_summaries[k]
        component_ivw[k]    = ivw(b_X, bY, sY)
        eg_slope, eg_se, eg_p, eg_int_p = mr_egger(b_X, bY, sY)
        component_eg[k]     = (eg_slope, eg_se, eg_p)   # ignore intercept_p for pooling
        component_wm[k]     = weighted_median_mr(b_X, bY, sY, n_boot=100, rng=rng)
        component_presso[k] = mr_presso_global(b_X, bY, sY, n_perm=200, rng=rng)

    # Pool
    out = {}
    for est_name, comp_dict in [("IVW", component_ivw),
                                ("Egger", component_eg),
                                ("WM", component_wm)]:
        comp_list = [comp_dict[k] for k in (1, 2, 3)]
        th_iv, se_iv, p_iv     = pool_inverse_variance(comp_list)
        th_bf, se_bf, p_bf     = pool_bonferroni(comp_list)
        p_fisher               = pool_fisher(comp_list)
        out[f"{est_name}_pool_IV"]     = {"theta": th_iv, "se": se_iv, "p": p_iv}
        out[f"{est_name}_pool_Bonf"]   = {"theta": th_bf, "se": se_bf, "p": p_bf}
        out[f"{est_name}_pool_Fisher"] = {"theta": np.nan, "se": np.nan, "p": p_fisher}

    # MR-PRESSO global on priority-1 only (treated as the primary endpoint)
    th_pr, se_pr, p_pr_ivw, p_pr_global = component_presso[1]
    out["PRESSO_p1_global"] = {"theta": th_pr, "se": se_pr, "p_ivw": p_pr_ivw,
                               "p_global": p_pr_global}

    return out

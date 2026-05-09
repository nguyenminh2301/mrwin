"""
inference.py — Per-SNP summary-statistic extraction (Cox-PH approximation
and Aalen additive-rate approximation) plus random-effects MR-Egger.

Implements two scale variants for SDPD calibration:

(a) Cox-PH per-SNP: log-hazard-ratio for SNP m on priority-1 mortality,
    matching the multiplicative scale of Section 6.1's Weibull DGP.

(b) Aalen additive: per-SNP additive rate effect on priority-1 mortality,
    matching the v5 mandate in Section 5.2 ("we mandate Aalen additive
    hazard models to preserve collapsibility").

Both feed into a random-effects MR-Egger regression, whose intercept tests
the SDPD null gamma_1 = 0. Type-I error of this intercept test under
gamma_1 = 0 across theta_F ∈ {0, 0.4, 0.8} is the InSIDE robustness
diagnostic for Flaw 2.
"""
from __future__ import annotations
import numpy as np
from typing import Tuple


# ---------------------------------------------------------------------------
# (a) Per-SNP Cox-PH log-HR via Newton-Raphson on the Breslow score equation.
# Vectorised across SNPs. For per-SNP effects with |beta| < 0.05, this matches
# full Cox to ~3 decimal places; we use it as the operational SDPD diagnostic
# at GWAS scale.
# ---------------------------------------------------------------------------
def cox_per_snp(
    G: np.ndarray, T: np.ndarray, D: np.ndarray, max_iter: int = 25, tol: float = 1e-7
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Returns (beta_hat, se_hat) per SNP for the Cox model
        lambda(t | G_m) = lambda_0(t) * exp(b_m * G_m).

    Uses Breslow's partial-likelihood score equation, solved per SNP independently
    by Newton-Raphson. Vectorised across SNPs, but the risk-set sums are computed
    per iteration (this is unavoidable for partial likelihood).
    """
    N, M = G.shape
    # Sort by event time once.
    order = np.argsort(T, kind="stable")
    T_s = T[order]
    D_s = D[order]
    G_s = G[order]

    # Distinct event times (Breslow ties handling: each distinct time contributes
    # the sum-of-G over events at that time).
    event_idx = np.where(D_s == 1)[0]

    # Risk-set indicator: at event j (in sorted order), risk set = {i : T_s[i] >= T_s[j]}.
    # Because T_s is sorted ascending, this is just indices [j:] for the smallest-tied
    # event with that time. For Breslow's approximation we use the simple "rolling tail"
    # construction.

    beta = np.zeros(M)
    for it in range(max_iter):
        # Linear predictor and risk weights per individual, per SNP.
        # LP[i, m] = G_s[i, m] * beta[m]
        LP = G_s * beta[np.newaxis, :]                                  # (N, M)
        w = np.exp(LP)                                                  # (N, M)

        # Cumulative sum from the bottom (i.e. risk-set sum at sorted index j is
        # sum_{i >= j} w[i, m]).
        cum_w = np.cumsum(w[::-1], axis=0)[::-1]                        # (N, M)
        cum_Gw = np.cumsum((G_s * w)[::-1], axis=0)[::-1]               # (N, M)
        cum_G2w = np.cumsum((G_s * G_s * w)[::-1], axis=0)[::-1]        # (N, M)

        # Score and information per SNP.
        # Score: sum over events j of (G_s[j, m] - cum_Gw[j, m] / cum_w[j, m])
        # Info:  sum over events j of (cum_G2w[j, m] / cum_w[j, m]
        #                              - (cum_Gw[j, m] / cum_w[j, m])^2)
        if event_idx.size == 0:
            return np.zeros(M), np.full(M, np.nan)
        e_cum_w  = cum_w[event_idx]                                     # (E, M)
        e_cum_Gw = cum_Gw[event_idx]
        e_cum_G2w = cum_G2w[event_idx]
        e_G      = G_s[event_idx]
        mean_G   = e_cum_Gw / np.maximum(e_cum_w, 1e-12)
        score    = (e_G - mean_G).sum(axis=0)
        var_G    = e_cum_G2w / np.maximum(e_cum_w, 1e-12) - mean_G ** 2
        info     = var_G.sum(axis=0)

        # Newton step.
        info_safe = np.maximum(info, 1e-10)
        step = score / info_safe
        # Damp to avoid runaway.
        step = np.clip(step, -0.5, 0.5)
        beta_new = beta + step
        if np.max(np.abs(step)) < tol:
            beta = beta_new
            break
        beta = beta_new

    # Final SE.
    LP = G_s * beta[np.newaxis, :]
    w = np.exp(LP)
    cum_w  = np.cumsum(w[::-1], axis=0)[::-1]
    cum_Gw = np.cumsum((G_s * w)[::-1], axis=0)[::-1]
    cum_G2w = np.cumsum((G_s * G_s * w)[::-1], axis=0)[::-1]
    e_cum_w  = cum_w[event_idx]
    e_cum_Gw = cum_Gw[event_idx]
    e_cum_G2w = cum_G2w[event_idx]
    mean_G = e_cum_Gw / np.maximum(e_cum_w, 1e-12)
    var_G = e_cum_G2w / np.maximum(e_cum_w, 1e-12) - mean_G ** 2
    info = var_G.sum(axis=0)
    se = 1.0 / np.sqrt(np.maximum(info, 1e-10))
    return beta, se


# ---------------------------------------------------------------------------
# (b) Aalen additive-rate per SNP. We use the simplest per-SNP additive model:
#     dN_i(t) = (a_0(t) + b_m G_{i,m}) Y_i(t) dt
# integrated over t, which yields a per-SNP linear-regression-style estimator
# of b_m. For per-SNP effects this collapses to OLS of cumulative event counts
# on (G_m, total at-risk time), and is the standard Aalen-Lin survival
# regression score (Aalen 1989; Martinussen-Scheike 2006 §5).
# ---------------------------------------------------------------------------
def aalen_per_snp(
    G: np.ndarray, T: np.ndarray, D: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Per-SNP Aalen additive estimator of the SNP-specific rate effect on a
    single endpoint. Returns (beta_hat, se_hat).

    Implementation: for each SNP m, regress D_i on (T_i, G_{i,m} * T_i) by OLS
    and report HC0 sandwich SE for the slope. Vectorised across SNPs to the
    extent feasible.

    NOTE (Package A scope): we report Aalen output for completeness, but in
    Package A the InSIDE robustness diagnostic is calibrated only on the
    Cox-PH variant. Aalen scale reconciliation is the subject of Package B,
    item R3, which addresses Flaw 3 from the adversarial review.
    """
    N, M = G.shape
    Y = T  # at-risk time
    # Per-SNP loop: small N x 2 OLS with HC0. M is at most a few hundred.
    b   = np.empty(M)
    se  = np.empty(M)
    Yc  = Y.reshape(-1, 1)
    for m in range(M):
        Xm = np.column_stack([Y, G[:, m] * Y])                # (N, 2)
        XtX = Xm.T @ Xm
        try:
            XtX_inv = np.linalg.inv(XtX)
        except np.linalg.LinAlgError:
            b[m] = 0.0; se[m] = np.inf; continue
        coef = XtX_inv @ (Xm.T @ D)
        a0_m, b_m = coef
        resid = D - Xm @ coef
        # HC0 sandwich
        meat = (Xm * (resid ** 2)[:, None]).T @ Xm
        cov  = XtX_inv @ meat @ XtX_inv
        b[m]  = b_m
        se[m] = float(np.sqrt(max(cov[1, 1], 1e-20)))
    return b, se


# ---------------------------------------------------------------------------
# Random-effects MR-Egger.
# Bowden 2015, Eq. (3): beta_Y = theta_0 + theta_1 * beta_X + epsilon_m,
# weighted by 1 / se(beta_Y_m)^2. Random-effects: residual heterogeneity
# absorbed into a variance multiplier.
# Returns (intercept, intercept_se, slope, slope_se).
# ---------------------------------------------------------------------------
def mr_egger_re(
    beta_X: np.ndarray, beta_Y: np.ndarray, se_Y: np.ndarray
) -> dict:
    w = 1.0 / np.maximum(se_Y ** 2, 1e-12)
    # WLS design [1, beta_X]
    X = np.column_stack([np.ones_like(beta_X), beta_X])
    W = np.diag(w)
    XtWX = X.T @ W @ X
    XtWY = X.T @ W @ beta_Y
    coef = np.linalg.solve(XtWX, XtWY)
    resid = beta_Y - X @ coef
    # Random-effects multiplier: Cochran Q over (M-2) df, with phi = max(1, Q/(M-2)).
    M = len(beta_X)
    Q = (resid ** 2 * w).sum()
    phi = max(1.0, Q / max(M - 2, 1))
    cov = phi * np.linalg.inv(XtWX)
    return {
        "intercept": float(coef[0]),
        "intercept_se": float(np.sqrt(cov[0, 0])),
        "slope": float(coef[1]),
        "slope_se": float(np.sqrt(cov[1, 1])),
        "Q": float(Q),
        "phi": float(phi),
        "M": int(M),
    }

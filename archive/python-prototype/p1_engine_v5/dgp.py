"""
dgp.py — Illness-death multi-state DGP with shared gamma frailty.

State diagram (matches Section 6.1 of CausalWinStatistic_v5.md):

    State 0 (healthy) ──► State 1 (HF hosp, priority 2)
       │                       │
       │                       ▼
       ├───────────► State 3 (death, priority 1, absorbing)
       │                       ▲
       ▼                       │
    State 2 (renal, priority 3)─┘

All transition intensities multiplied by a shared frailty W ~ Gamma(1/theta_F, theta_F),
which couples mortality with both intermediate states. Under theta_F > 0,
ascertainment on survival to GWAS analysis induces correlation between
SNP–exposure and SNP–mortality summary statistics — the InSIDE violation
of Flaw 2 in the v5 review.

Notation matches Section 6.1:
    lambda_{rs}(t | X, U, S, W) = W * lambda_{0,rs} * nu * t^{nu-1}
                                  * exp(alpha_{rs} X + nu_{rs} U + gamma_{rs} S)

For simplicity and consistency with the v5 simulation specification,
we use the same alpha, nu, gamma vectors per priority across all transitions
that involve that priority (e.g. transitions 0->3 and 1->3, 2->3 share the
priority-1 parameters).
"""
from __future__ import annotations
import numpy as np
from .config import P1Config


def _weibull_inv(rng, rate, shape, lp, N):
    """Sample T ~ Weibull(rate * exp(lp), shape) by inverse CDF."""
    u = rng.uniform(0.0, 1.0, N)
    return (-np.log(u) / (rate * np.exp(lp))) ** (1.0 / shape)


def simulate_cohort(cfg: P1Config, rng: np.random.Generator) -> dict:
    """
    Simulate a single outcome cohort under the illness-death DGP with shared
    frailty. Returns observed per-priority event times and indicators that
    feed both the cCWR estimator (not used in Package A) and the per-SNP
    summary-statistic extraction used by SDPD.

    Returned arrays:
        G        : (N, M)    standardised genotype matrix
        true_betas: (M,)     per-SNP true effects on X (target of GWAS)
        S_true   : (N,)      latent PRS = G @ true_betas
        X        : (N,)      continuous exposure
        U        : (N,)      unmeasured confounder ~ N(0,1)
        W        : (N,)      shared frailty (or ones if theta_F = 0)
        T        : (N, 3)    observed times per priority (1=death, 2=HF, 3=renal)
        D        : (N, 3)    observed event indicators per priority
        C        : (N,)      administrative censoring time
    """
    N = cfg.N_outcome
    M = cfg.M_snps
    K = 3

    # Genotypes
    mafs = rng.uniform(cfg.maf_low, cfg.maf_high, M)
    G = np.column_stack([rng.binomial(2, p, N) for p in mafs]).astype(np.float64)
    G_std = (G - G.mean(axis=0)) / (G.std(axis=0) + 1e-12)

    # True per-SNP effects on X (the GWAS target)
    true_betas = rng.normal(0.0, cfg.sigma_beta, M)
    S_true = G_std @ true_betas

    # Confounder, frailty, exposure
    U = rng.normal(0.0, 1.0, N)
    if cfg.theta_F > 0.0:
        # Gamma(shape=1/theta_F, scale=theta_F): mean 1, var theta_F.
        W = rng.gamma(shape=1.0 / cfg.theta_F, scale=cfg.theta_F, size=N)
    else:
        W = np.ones(N)
    X = cfg.alpha_S * S_true + cfg.alpha_U * U + rng.normal(0.0, 1.0, N)

    # Latent cause-specific times for the 0->{1,2,3} transitions out of the
    # initial (healthy) state. The shared frailty W enters multiplicatively
    # via log W in the linear predictor.
    log_W = np.log(W)
    lp = []
    for k in range(K):
        lp_k = (
            cfg.alpha_X[k] * X
            + cfg.nu_U[k] * U
            + cfg.gamma_direct[k] * S_true
            + log_W
        )
        lp.append(lp_k)

    # Map: priority 1 = death (index 0 in alpha_X tuple if we follow the v5 order
    # priority = (death, HF, renal). The DGP's baseline_haz tuple in the same
    # order: baseline_haz[0]=death, [1]=HF, [2]=renal. We adopt this convention
    # consistently here.
    t_death = _weibull_inv(rng, cfg.baseline_haz[0], cfg.shape_weibull, lp[0], N)
    t_hf    = _weibull_inv(rng, cfg.baseline_haz[1], cfg.shape_weibull, lp[1], N)
    t_renal = _weibull_inv(rng, cfg.baseline_haz[2], cfg.shape_weibull, lp[2], N)

    # Administrative censoring
    C = rng.exponential(1.0 / cfg.censoring_rate, N)
    C = np.minimum(C, cfg.max_follow_up)

    # Priority-1 (death): observed = min(t_death, C)
    T1 = np.minimum(t_death, C)
    D1 = (t_death <= C).astype(int)

    # Priority-2 (HF): can only occur before death OR before censoring
    T2 = np.minimum(np.minimum(t_hf, C), t_death)
    D2 = ((t_hf <= C) & (t_hf <= t_death)).astype(int)

    # Priority-3 (renal): same logic
    T3 = np.minimum(np.minimum(t_renal, C), t_death)
    D3 = ((t_renal <= C) & (t_renal <= t_death)).astype(int)

    return {
        "G": G_std,
        "mafs": mafs,
        "true_betas": true_betas,
        "S_true": S_true,
        "U": U,
        "W": W,
        "X": X,
        "T": np.column_stack([T1, T2, T3]),
        "D": np.column_stack([D1, D2, D3]),
        "C": C,
    }

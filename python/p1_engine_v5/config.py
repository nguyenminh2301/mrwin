"""
config.py — P1 v5 simulation engine, Package A configuration.

Single source of truth for every tunable. asdict(config) is written
verbatim into the JSON output for audit purposes.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Tuple


@dataclass
class P1Config:
    # ---- DGP: illness-death with shared frailty (Section 6.1 of v5) ----
    N_outcome: int = 6_000
    M_snps: int = 40
    alpha_S: float = 0.4              # instrument strength on X
    alpha_U: float = 0.8              # confounder strength on X
    # Per-priority structural log-hazard ratios on X.
    # Priority order: (priority-1 = death/absorbing,
    #                  priority-2 = HF hospitalisation,
    #                  priority-3 = renal decline)
    alpha_X: Tuple[float, float, float] = (-0.4, -0.4, -0.4)
    # Hierarchy-contaminating pleiotropy on each priority.
    # gamma[0] = direct mortality pleiotropy (gamma_1 in v5 notation).
    gamma_direct: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    nu_U: Tuple[float, float, float] = (0.5, 0.5, 0.5)
    baseline_haz: Tuple[float, float, float] = (0.02, 0.05, 0.10)
    shape_weibull: float = 1.2
    censoring_rate: float = 0.05
    max_follow_up: float = 10.0

    # ---- Shared frailty: W ~ Gamma(1/theta_F, theta_F), E[W]=1, Var[W]=theta_F.
    # theta_F = 0 disables frailty (degenerate W = 1).
    theta_F: float = 0.8

    # ---- Ascertainment: if True, restrict the analysis cohort to individuals
    # surviving past the median observed mortality time. Mirrors a prevalent-
    # sample GWAS where cases are oversampled relative to early decedents.
    ascertain_on_survival: bool = False

    # ---- MC and bootstrap iteration counts ----
    n_iter: int = 400
    B_bootstrap: int = 200            # not used in Package A (no MB needed)
    seed: int = 20260506

    # ---- SNP genetics ----
    maf_low: float = 0.10
    maf_high: float = 0.40
    sigma_beta: float = 0.05          # SD of true per-SNP effect on X


# ---- Three calibration scenarios for Package A R2 (InSIDE robustness) ----
def insider_scenarios(N_outcome: int | None = None,
                      M_snps: int | None = None,
                      n_iter: int | None = None,
                      seed: int | None = None) -> dict:
    """
    Returns six scenarios spanning theta_F in {0, 0.4, 0.8} crossed with
    "ascertained" in {False, True}. Optional overrides for N, M, n_iter, seed
    are applied to the base config; otherwise the dataclass defaults are used.

    Ascertainment: the post-hoc review (Flaw 2 of v5 review) attributes the
    InSIDE violation to "ascertainment on survival to GWAS analysis". Under
    the manuscript's exact one-sample, no-ascertainment setup this leakage
    does not materialise; we demonstrate this by running both variants:

      * ascertained = False  → manuscript's actual setup. Expected: nominal
        Type-I error across all theta_F. Confirms SDPD validity in v5's
        precise data-flow.
      * ascertained = True   → restrict to individuals surviving past the
        median observed event time before fitting MR-Egger. Expected:
        Type-I error inflation at theta_F > 0. Maps the boundary at
        which InSIDE fails — relevant when the GWAS cohort is itself a
        prevalent (long-survival) sample.

    Both sets together resolve Flaw 2 quantitatively.
    """
    base_kwargs = {}
    if N_outcome is not None: base_kwargs["N_outcome"] = N_outcome
    if M_snps    is not None: base_kwargs["M_snps"]    = M_snps
    if n_iter    is not None: base_kwargs["n_iter"]    = n_iter
    if seed      is not None: base_kwargs["seed"]      = seed
    base = P1Config(**base_kwargs)
    return {
        "F0_no_frailty_unascertained":   _replace(base, theta_F=0.0, ascertain_on_survival=False),
        "F1_moderate_unascertained":     _replace(base, theta_F=0.4, ascertain_on_survival=False),
        "F2_v5_unascertained":           _replace(base, theta_F=0.8, ascertain_on_survival=False),
        "F0_no_frailty_ascertained":     _replace(base, theta_F=0.0, ascertain_on_survival=True),
        "F1_moderate_ascertained":       _replace(base, theta_F=0.4, ascertain_on_survival=True),
        "F2_v5_ascertained":             _replace(base, theta_F=0.8, ascertain_on_survival=True),
    }


def _replace(cfg: P1Config, **kw) -> P1Config:
    from dataclasses import replace
    return replace(cfg, **kw)

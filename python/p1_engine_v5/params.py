"""
params.py — Centralised replication-vs-publication parameter switchboard.

Two parameter sets:

  REPLICATION       — exact parameters Claude executed in the May 2026
                      revision session on a 1-CPU container. Reproduces
                      every number reported in:
                        - results_packageA_R1.json (analytical, no MC)
                        - results_packageA_R2.json
                        - results_packageB_R3.json
                        - results_packageB_R4.json
                        - results_packageC_item1.json
                        - results_packageC_item2.json
                      Total wall-clock for the full sweep on 1 CPU: ~5–7 min.
                      For confirming that the engine is intact and bug-free
                      before scaling up.

  PUBLICATION       — parameters recommended for the v5.1 → IJE submission.
                      Corresponds to the v5 manuscript's canonical scale
                      (N=10,000 per cohort, M=50 SNPs, 1,000 MC iterations,
                      B=1,000 bootstrap iterations).
                      Total wall-clock estimate on a 32-core cluster:
                      ~2–3 hours for the full sweep.
                      For the actual numbers reported in the submitted
                      manuscript and supplements.

USAGE
=====

From any engine module:

    from p1_engine_v5.params import get_params, MODE_REPLICATION, MODE_PUBLICATION
    p = get_params(mode="publication")     # or "replication"
    cfg = P1Config(N_outcome=p["N_outcome"], M_snps=p["M_snps"], n_iter=p["n_iter"], ...)

From the CLI (every engine entry-point supports --mode):

    python -m p1_engine_v5.engine --mode replication
    python -m p1_engine_v5.engine_packageB --task R3 --mode publication
    python -m p1_engine_v5.q_statistic_asymptotics --mode publication

NOTE ON MODE INTERACTIONS
=========================

- Replication-mode runs use Claude's audited seeds and produce *bitwise
  identical* results to results_packageA_R2.json etc. on any platform with
  numpy >= 2.0 (verify by checksum after running).
- Publication-mode runs share the same seeding scheme (cfg.seed * 10_000 + k)
  but at scaled-up N, M, n_iter, B; results are NOT bitwise identical to
  replication, but Monte-Carlo conclusions (Type-I error within 0.005 of
  reported) should hold.
- Both modes use the same DGP (dgp.py), the same SDPD (inference.py), the
  same multiplier bootstrap (multiplier_bootstrap.py). The difference is
  parameter values only.
"""

MODE_REPLICATION = "replication"
MODE_PUBLICATION = "publication"


# =============================================================================
# REPLICATION  (what Claude actually ran on the 1-CPU container, May 2026)
# =============================================================================
REPLICATION = {
    # ---- DGP / cohort
    "N_outcome":       6_000,        # Cohort size
    "M_snps":          40,           # Number of independent SNPs
    "alpha_S":         0.4,          # Instrument strength
    "alpha_U":         0.8,          # Confounder strength
    "alpha_X_protective": (-0.4, -0.4, -0.4),
    "theta_F":         0.8,          # Frailty variance (matches Section 6.1)
    "shape_weibull":   1.2,
    "censoring_rate":  0.05,
    "max_follow_up":   10.0,
    "sigma_beta":      0.05,
    # ---- Package A R2 (InSIDE robustness, 6 scenarios)
    "A_R2_n_iter":     400,          # iterations per scenario
    # ---- Package B R3 (scale reconciliation across gamma_1 grid)
    "B_R3_n_iter":     300,
    "B_R3_gamma_grid": [0.00, 0.02, 0.05, 0.10, 0.15],
    # ---- Package B R4 (bootstrap moment diagnostics, 3 scenarios)
    "B_R4_B":          500,          # bootstrap iterations
    # ---- Package C Item 1 (Q-statistic parametric Gaussian validation)
    "C_C1_B":          500,          # multiplier bootstrap iterations
    "C_C1_n_perm":     1000,         # parametric Gaussian draws
    # ---- Package C Item 2 (wall-clock benchmark)
    "C_C2_grid_N":     [4_000, 6_000, 10_000],
    "C_C2_M":          50,
    "C_C2_B":          200,
    "C_C2_sparse_test_N": 10_000,    # for sparse-pair extrapolation
    # ---- Common
    "seed":            20_260_506,
    "D_strata":        10,
}


# =============================================================================
# PUBLICATION  (recommended for the v5.1 IJE submission)
# =============================================================================
PUBLICATION = {
    # ---- DGP / cohort  (matches v5 main text Section 6.1 declared scale)
    "N_outcome":       10_000,
    "M_snps":          50,
    "alpha_S":         0.4,
    "alpha_U":         0.8,
    "alpha_X_protective": (-0.4, -0.4, -0.4),
    "theta_F":         0.8,
    "shape_weibull":   1.2,
    "censoring_rate":  0.05,
    "max_follow_up":   10.0,
    "sigma_beta":      0.05,
    # ---- Package A R2: 1,000 iterations × 6 scenarios
    # MC SE at p=0.05, n=1000 ≈ 0.0069 (vs 0.011 at n=400).
    # Detects 2pp Type-I inflation cleanly.
    "A_R2_n_iter":     1_000,
    # ---- Package B R3: 1,000 iterations × 5 gamma_1 values
    # Same MC SE bound; resolves Cox-Aalen power gap to 1pp.
    "B_R3_n_iter":     1_000,
    "B_R3_gamma_grid": [0.00, 0.02, 0.05, 0.10, 0.15],
    # ---- Package B R4: 1,000 bootstrap iterations
    # Kurtosis SE at B=1000 ~ sqrt(24/B) = 0.155, vs 0.219 at B=500;
    # tightens the kurtosis estimate enough to call "Gaussian" with confidence.
    "B_R4_B":          1_000,
    # ---- Package C Item 1: 1,000 bootstrap × 10,000 parametric Gaussian draws
    # KS-test resolution ≈ sqrt(2/n) = 0.014; detects 0.05 KS deviation reliably.
    "C_C1_B":          1_000,
    "C_C1_n_perm":     10_000,
    # ---- Package C Item 2: extended grid + cluster benchmark target
    "C_C2_grid_N":     [10_000, 50_000, 100_000],
    "C_C2_M":          50,
    "C_C2_B":          1_000,
    "C_C2_sparse_test_N": 100_000,   # closer to UK Biobank for tighter projection
    # ---- Common
    "seed":            20_260_506,
    "D_strata":        10,
}


def get_params(mode: str) -> dict:
    """Return parameter dictionary for the requested mode."""
    if mode == MODE_REPLICATION:
        return dict(REPLICATION)
    if mode == MODE_PUBLICATION:
        return dict(PUBLICATION)
    raise ValueError(f"Unknown mode: {mode!r}. Use 'replication' or 'publication'.")


def add_mode_arg(parser):
    """Helper: add --mode flag to an argparse parser."""
    parser.add_argument(
        "--mode", choices=["replication", "publication"], default="replication",
        help="Parameter set: 'replication' (Claude's audited 1-CPU run, ~5 min) "
             "or 'publication' (canonical IJE-submission scale, ~2–3 h on 32-core cluster)."
    )

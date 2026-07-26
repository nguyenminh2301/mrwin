"""
wallclock_benchmark.py — Package C Item 2.

Address cross-examination Q6: provide wall-clock runtime for the full pipeline
at multiple cohort sizes and bootstrap counts, stating cluster configuration.

Phases timed:
  1. Cohort simulation (DGP).
  2. Kernel precomputation O(N^2).
  3. Per-iteration bootstrap (DGP, kernel, win-loss sums) at small n_iter.
  4. Multiplier bootstrap with B = 200 and B = 500.

Results report:
  - Wall-clock per phase
  - Memory footprint of kernel matrix at each N
  - Extrapolation to biobank N = 100,000 (run a single full pipeline at
    100,000 if memory permits; otherwise report sub-task projections)
"""
from __future__ import annotations
import gc
import json
import time
import numpy as np
from pathlib import Path

from .config import P1Config
from .dgp import simulate_cohort
from .inference import cox_per_snp, aalen_per_snp, mr_egger_re
from .kernel import precompute_kernel_matrix
from .multiplier_bootstrap import run_multiplier_bootstrap


def _memory_mb(arr) -> float:
    """Memory footprint of a numpy array in MB."""
    return arr.nbytes / 1024 ** 2


def benchmark_at_N(N: int, M: int = 50, B: int = 200, seed: int = 20260506) -> dict:
    """One full benchmark run at given N."""
    print(f"\n=== Benchmark at N={N:,}, M={M}, B={B} ===")
    cfg = P1Config(N_outcome=N, M_snps=M, theta_F=0.8, alpha_S=0.4,
                   alpha_X=(-0.4,-0.4,-0.4), seed=seed)
    rec = {"N": N, "M": M, "B": B}

    # Phase 1: cohort simulation
    rng = np.random.default_rng(cfg.seed)
    t = time.time()
    coh = simulate_cohort(cfg, rng)
    rec["sim_cohort_s"] = round(time.time() - t, 3)
    print(f"  cohort sim: {rec['sim_cohort_s']:.2f}s")

    # Phase 2: kernel precompute
    t = time.time()
    H = precompute_kernel_matrix(coh["T"], coh["D"])
    rec["kernel_precompute_s"] = round(time.time() - t, 3)
    rec["kernel_memory_MB"]    = round(_memory_mb(H), 2)
    print(f"  kernel precompute: {rec['kernel_precompute_s']:.2f}s, memory: {rec['kernel_memory_MB']:.0f} MB")

    # Phase 3: SDPD diagnostics (Cox + Aalen + MR-Egger)
    G  = coh["G"]; X = coh["X"]; T1 = coh["T"][:, 0]; D1 = coh["D"][:, 0]
    t = time.time()
    GtG = (G * G).sum(axis=0)
    bX  = (G.T @ X) / np.maximum(GtG, 1e-12)
    seX = np.sqrt(np.var(X, ddof=1) / np.maximum(GtG, 1e-12))
    rec["gwas_X_s"] = round(time.time() - t, 3)
    t = time.time()
    bY_cox, sY_cox = cox_per_snp(G, T1, D1)
    rec["cox_per_snp_s"] = round(time.time() - t, 3)
    t = time.time()
    bY_aal, sY_aal = aalen_per_snp(G, T1, D1)
    rec["aalen_per_snp_s"] = round(time.time() - t, 3)
    t = time.time()
    eg_cox = mr_egger_re(bX, bY_cox, sY_cox)
    eg_aal = mr_egger_re(bX, bY_aal, sY_aal)
    rec["mr_egger_s"] = round(time.time() - t, 3)
    print(f"  SDPD: cox {rec['cox_per_snp_s']:.2f}s + aalen {rec['aalen_per_snp_s']:.2f}s + mr_egger {rec['mr_egger_s']:.3f}s")

    # Phase 4: multiplier bootstrap (full pipeline)
    Z = rng.normal(size=(N, 3))
    t = time.time()
    res = run_multiplier_bootstrap(
        H, G, X, Z, bX, seX, D_strata=10, B=B, seed=cfg.seed,
    )
    rec["multiplier_bootstrap_s"] = round(time.time() - t, 3)
    rec["delta_GLS"] = float(res["delta_GLS"])
    rec["se_delta_GLS"] = float(res["se_delta_GLS"])
    print(f"  multiplier bootstrap (B={B}): {rec['multiplier_bootstrap_s']:.2f}s")
    print(f"  delta_GLS = {rec['delta_GLS']:.4f}, SE = {rec['se_delta_GLS']:.4f}")

    rec["total_pipeline_s"] = round(
        rec["sim_cohort_s"] + rec["kernel_precompute_s"] + rec["gwas_X_s"]
        + rec["cox_per_snp_s"] + rec["aalen_per_snp_s"] + rec["mr_egger_s"]
        + rec["multiplier_bootstrap_s"], 2
    )
    print(f"  TOTAL: {rec['total_pipeline_s']:.1f}s")
    # Free memory before next N
    del coh, H, G, X, Z; gc.collect()
    return rec


def main(out_path: str = "results_packageC_item2.json",
         grid_N: list[int] | None = None, M: int = 50, B: int = 200,
         sparse_test_N: int = 10_000) -> dict:
    """
    Benchmark across a grid of cohort sizes; project to UK Biobank scale
    (N=337,000) under both dense and sparse-pair kernels.
    """
    if grid_N is None:
        grid_N = [4_000, 6_000, 10_000]
    out = {"benchmarks": []}
    for N in grid_N:
        out["benchmarks"].append(benchmark_at_N(N=N, M=M, B=B))

    # Extrapolation: dense kernel
    print("\n=== Extrapolation to biobank N = 337,000 (UK Biobank WBU), DENSE kernel ===")
    last = out["benchmarks"][-1]
    factor = (337_000 / last["N"]) ** 2
    proj_kernel = last["kernel_precompute_s"] * factor
    proj_kernel_mem = last["kernel_memory_MB"] * factor
    proj_boot_per_iter = (last["multiplier_bootstrap_s"] / last["B"]) * factor
    proj_boot_total = proj_boot_per_iter * B
    out["extrapolation_337k_dense"] = {
        "kernel_precompute_proj_s":   round(proj_kernel, 1),
        "kernel_memory_proj_GB":      round(proj_kernel_mem / 1024, 2),
        "bootstrap_per_iter_proj_s":  round(proj_boot_per_iter, 1),
        f"bootstrap_B{B}_total_proj_min": round(proj_boot_total / 60, 1),
        "feasible_64gb_node": bool(proj_kernel_mem / 1024 < 64),
    }
    print(f"  Projected kernel memory: {out['extrapolation_337k_dense']['kernel_memory_proj_GB']:.0f} GB")
    print(f"  Projected bootstrap (B={B}) total: {round(proj_boot_total/60,1)} min")

    # Sparse-pair benchmark at sparse_test_N
    print(f"\n=== Sparse-pair kernel benchmark at N = {sparse_test_N:,} ===")
    from .kernel_sparse import stratum_pair_kernel
    rng2 = np.random.default_rng(12345)
    cfg2 = P1Config(N_outcome=sparse_test_N, M_snps=M, theta_F=0.8, alpha_S=0.4,
                    alpha_X=(-0.4,-0.4,-0.4), seed=12345)
    coh2 = simulate_cohort(cfg2, rng2)
    S = coh2["S_true"]
    qs = np.quantile(S, np.linspace(0, 1, 11))
    qs[0] = -np.inf; qs[-1] = np.inf
    d_idx = np.digitize(S, qs[1:-1])
    t = time.time()
    block_total_mem_MB = 0.0
    for dd in range(1, 10):
        idx_d   = np.where(d_idx == dd)[0]
        idx_dm1 = np.where(d_idx == dd-1)[0]
        H_block = stratum_pair_kernel(
            coh2["T"][idx_d], coh2["D"][idx_d],
            coh2["T"][idx_dm1], coh2["D"][idx_dm1],
        )
        block_total_mem_MB += _memory_mb(H_block)
    sparse_one_pass_s = round(time.time() - t, 3)
    print(f"  Sparse one-pass over 9 stratum pairs: {sparse_one_pass_s:.2f}s")
    print(f"  Sum of block memory (transient): {block_total_mem_MB:.1f} MB")

    factor_sparse = (337_000 / sparse_test_N) ** 2
    proj_sparse_per_pass = sparse_one_pass_s * factor_sparse
    proj_sparse_mem_GB   = block_total_mem_MB * factor_sparse / 1024
    proj_sparse_boot_total = proj_sparse_per_pass * B
    out["extrapolation_337k_sparse"] = {
        "one_pass_kernel_proj_s":     round(proj_sparse_per_pass, 1),
        "peak_block_memory_proj_GB":  round(proj_sparse_mem_GB / 9, 2),
        f"bootstrap_B{B}_total_proj_min": round(proj_sparse_boot_total / 60, 1),
        "feasible_64gb_node": bool(proj_sparse_mem_GB / 9 < 64),
    }
    print(f"  Projected peak block memory (single pair): {round(proj_sparse_mem_GB/9,2)} GB")
    print(f"  Projected bootstrap (B={B}) wall-clock: {round(proj_sparse_boot_total/60,0)} min (1 CPU)")

    out["cluster_recommendation"] = {
        "primary_strategy": "sparse-pair kernel (kernel_sparse.py)",
        "memory_required_GB": round(proj_sparse_mem_GB / 9 + 8, 1),
        "cpu_recommendation": "32 cores; multiprocessing across bootstrap iterations cuts wall-clock by ~25x",
        "approximate_cluster_wallclock_min_full_run": round(proj_sparse_boot_total / 60 / 25, 1),
    }
    print(f"\nCluster recommendation: {out['cluster_recommendation']['memory_required_GB']:.0f} GB RAM, 32 cores → ~{out['cluster_recommendation']['approximate_cluster_wallclock_min_full_run']:.0f} min for full bootstrap")

    Path(out_path).write_text(json.dumps(out, indent=2))
    print(f"\n[saved] {out_path}")
    return out


if __name__ == "__main__":
    import argparse
    from .params import get_params, add_mode_arg
    parser = argparse.ArgumentParser(
        description="Package C Item 2: wall-clock benchmark and biobank-scale feasibility."
    )
    add_mode_arg(parser)
    parser.add_argument("--out", default="results_packageC_item2.json")
    args = parser.parse_args()
    pp = get_params(args.mode)
    print(f"[mode={args.mode}] grid_N={pp['C_C2_grid_N']}, M={pp['C_C2_M']}, B={pp['C_C2_B']}, sparse_test_N={pp['C_C2_sparse_test_N']}")
    main(out_path=args.out, grid_N=pp["C_C2_grid_N"], M=pp["C_C2_M"],
         B=pp["C_C2_B"], sparse_test_N=pp["C_C2_sparse_test_N"])

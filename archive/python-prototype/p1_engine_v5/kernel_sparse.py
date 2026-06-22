"""
kernel_sparse.py — Memory-efficient kernel: precompute only adjacent-stratum
pairs.

For the cCWR estimator with D deciles, only pairs (i, j) with i in stratum d
and j in stratum d-1 are ever evaluated. Storing the full N×N kernel matrix
wastes a factor of D in memory. At biobank scale (N = 337,000, D = 10) this
is the difference between 106 GB (infeasible) and 10.6 GB (feasible on a
high-memory node).

Trade-off versus kernel.py:
  + 1/D-fold memory reduction.
  - Kernel must be recomputed when stratification changes (each bootstrap
    iteration). For the multiplier bootstrap at B = 200 this multiplies
    bootstrap runtime by D / 1 = factor of D. NOT a free lunch.
  - HOWEVER: the per-iteration kernel block is only (N/D) × (N/D), so the
    pairwise comparison cost per iteration is (N/D)^2 rather than N^2.
    Total work = B × (D-1) × (N/D)^2 = B (D-1)/D^2 × N^2.
    Compared to the dense kernel approach which does B × (D-1) × (N/D)^2
    PLUS one N^2 precompute, the asymptotic complexity is identical for
    moderate B; the dense approach has a setup cost amortised over B,
    while the sparse approach has no setup but slightly higher per-iter cost.

For Package C Item 2, the dense kernel.py remains the default for moderate-N
runs; this sparse_kernel.py module enables the biobank-scale feasibility
claim by handling the case where the dense matrix exceeds memory.
"""
from __future__ import annotations
import numpy as np


def stratum_pair_kernel(
    T_d:    np.ndarray, D_d:    np.ndarray,
    T_dm1:  np.ndarray, D_dm1:  np.ndarray,
) -> np.ndarray:
    """
    Compute the (n_d, n_{d-1}) int8 kernel block for one stratum pair only.

    Args:
        T_d, D_d   : (n_d, K) outcomes for stratum d.
        T_dm1, D_dm1 : (n_{d-1}, K) outcomes for stratum d-1.

    Returns:
        H_block    : (n_d, n_{d-1}) int8 matrix with values {-1, 0, +1}.
    """
    n_i, K = T_d.shape
    n_j    = T_dm1.shape[0]
    H = np.zeros((n_i, n_j), dtype=np.int8)
    decided = np.zeros((n_i, n_j), dtype=bool)

    for k in range(K):
        Tk_i = T_d[:, k][:,   None]   # (n_i, 1)
        Tk_j = T_dm1[:, k][None, :]   # (1, n_j)
        Dk_j = D_dm1[:, k][None, :]
        Dk_i = D_d[:, k][:,   None]
        i_wins  = (Dk_j == 1) & (Tk_i > Tk_j)
        i_loses = (Dk_i == 1) & (Tk_j > Tk_i)
        new_win  = i_wins  & ~decided
        new_loss = i_loses & ~decided
        H[new_win] = 1
        H[new_loss] = -1
        decided |= (new_win | new_loss)
    return H


def stratum_pair_logcwr(
    T_d:    np.ndarray, D_d:    np.ndarray,
    T_dm1:  np.ndarray, D_dm1:  np.ndarray,
    weights_d:   np.ndarray | None = None,
    weights_dm1: np.ndarray | None = None,
) -> float:
    """
    Compute weighted log(theta^{(d, d-1)}) without materialising the full
    N×N kernel. Bootstrap iterations call this directly.
    """
    H_block = stratum_pair_kernel(T_d, D_d, T_dm1, D_dm1)
    if weights_d is None and weights_dm1 is None:
        sum_wins   = float(np.sum(H_block == 1))
        sum_losses = float(np.sum(H_block == -1))
    else:
        wi = weights_d[:, None]
        wj = weights_dm1[None, :]
        W  = wi * wj
        sum_wins   = float(np.sum(W * (H_block == 1)))
        sum_losses = float(np.sum(W * (H_block == -1)))
    return float(np.log(max(sum_wins, 1e-12) / max(sum_losses, 1e-12)))

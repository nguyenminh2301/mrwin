"""
kernel.py — Vectorised hierarchical comparison kernel.

For the v5 manuscript's K=3 priority order (priority-1=death, priority-2=HF,
priority-3=renal) the kernel h(O_i, O_j) is computed by sequential adjudication:
the first informative priority terminates the comparison.

Since the kernel depends only on the observed event times/indicators and not on
stratum assignment, we precompute an (N, N) matrix once per cohort and reuse it
across all bootstrap iterations. Memory at N=6000: int8(6000, 6000) ≈ 36 MB.

Output convention:
    H[i, j] = +1  if individual i wins the pairwise comparison against j
    H[i, j] = -1  if i loses to j
    H[i, j] =  0  if the pair ties at all priorities

Antisymmetric: H[i, j] = -H[j, i] strictly.
"""
from __future__ import annotations
import numpy as np


def precompute_kernel_matrix(T: np.ndarray, D: np.ndarray, block_size: int = 4_000) -> np.ndarray:
    """
    Compute the full (N, N) int8 hierarchical-kernel matrix in row-blocks.

    Memory model: the dense H is int8 of size N^2 bytes. Transient intermediates
    are bool arrays of size block_size * N. For N=30,000 with block_size=4,000,
    transient memory is ~120 MB instead of ~900 MB.

    Args:
        T : (N, K) observed times per priority. Order: priority-1 = death,
            priority-2 = HF, priority-3 = renal.
        D : (N, K) event indicators per priority.

    Returns:
        H : (N, N) int8 matrix with values in {-1, 0, +1}. Diagonal is zero.
    """
    N, K = T.shape
    H = np.zeros((N, N), dtype=np.int8)

    # We'll process row-blocks [a:b] x [0:N] at each iteration.
    # We need the "decided" bool to persist across priority levels for each block.
    for a in range(0, N, block_size):
        b = min(a + block_size, N)
        decided = np.zeros((b - a, N), dtype=bool)
        # Diagonal entries within this block are "decided" (i==j)
        for i in range(a, b):
            decided[i - a, i] = True

        for k in range(K):
            Tk = T[:, k]; Dk = D[:, k]
            Tk_i = Tk[a:b, None]                # (block, 1)
            Tk_j = Tk[None, :]                  # (1, N)
            Dk_i = Dk[a:b, None]
            Dk_j = Dk[None, :]
            i_wins  = (Dk_j == 1) & (Tk_i > Tk_j)
            i_loses = (Dk_i == 1) & (Tk_j > Tk_i)
            new_win  = i_wins  & ~decided
            new_loss = i_loses & ~decided
            H_block = H[a:b]
            H_block[new_win]  = 1
            H_block[new_loss] = -1
            H[a:b] = H_block
            decided |= (new_win | new_loss)
            # Free intermediate bool arrays explicitly before next priority level
            del i_wins, i_loses, new_win, new_loss

    np.fill_diagonal(H, 0)
    return H


def stratum_win_loss(
    H: np.ndarray, idx_d: np.ndarray, idx_dm1: np.ndarray,
    weights: np.ndarray | None = None
) -> tuple[float, float, float]:
    """
    Compute weighted win/loss/tie sums for the stratum pair (d, d-1).

    Args:
        H        : precomputed (N, N) kernel matrix.
        idx_d    : indices of individuals in stratum d.
        idx_dm1  : indices of individuals in stratum d-1.
        weights  : optional (N,) per-individual weight vector (xi*omega).

    Returns:
        (sum_wins, sum_losses, sum_total)
        where for the (d, d-1) contrast,
          sum_wins = sum_{i in d, j in d-1} w_i w_j * 1{H[i,j]=+1}
          sum_losses = sum w_i w_j * 1{H[i,j]=-1}
          sum_total  = sum w_i w_j  (over all pairs, including ties)
    """
    H_block = H[np.ix_(idx_d, idx_dm1)]
    if weights is None:
        n_i = idx_d.size
        n_j = idx_dm1.size
        sum_wins   = float(np.sum(H_block == 1))
        sum_losses = float(np.sum(H_block == -1))
        sum_total  = float(n_i * n_j)
    else:
        wi = weights[idx_d][:, None]
        wj = weights[idx_dm1][None, :]
        W  = wi * wj
        sum_wins   = float(np.sum(W * (H_block == 1)))
        sum_losses = float(np.sum(W * (H_block == -1)))
        sum_total  = float(np.sum(W))

    return sum_wins, sum_losses, sum_total


def stratum_log_cwr(
    H: np.ndarray, idx_d: np.ndarray, idx_dm1: np.ndarray,
    weights: np.ndarray | None = None
) -> float:
    """
    Compute ln(theta^{(d, d-1)}) = ln(tau_1 / tau_{-1}) for one contrast.

    Returns:
        ln(theta) at the input stratum/weights, with a small floor on
        numerator and denominator to avoid log(0).
    """
    sw, sl, st = stratum_win_loss(H, idx_d, idx_dm1, weights)
    sw = max(sw, 1e-12)
    sl = max(sl, 1e-12)
    return float(np.log(sw / sl))

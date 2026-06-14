"""
kernel_fast.py — WP13 / Phase II S1 reference: subquadratic single-endpoint
(K=1) hierarchical win/loss.

The dense backend (kernel.py / kernel_sparse.py) compares every cross-pair, so a
stratum pair (d, d-1) costs O(n_d * n_{d-1}) and the bootstrap repeats this every
iteration: B * Theta(N^2 / D). For a single endpoint the win/loss counts
factorise into weighted dominance counts that a sorted sweep evaluates in
O(N log N) — independent of D — with bit-for-bit identical results.

Reduction (single endpoint; high set H = stratum d, low set L = stratum d-1):

    W   = sum_{j in L, status_j=1} w_j * ( sum_{i in H : t_i > t_j} w_i )
    Lo  = sum_{i in H, status_i=1} w_i * ( sum_{j in L : t_j > t_i} w_j )
    Tot = ( sum_{i in H} w_i ) * ( sum_{j in L} w_j )

with STRICT time inequalities (equal times tie), matching the dense kernel rule:

    high i wins  over low j  <=>  status_low_j = 1 and t_high_i > t_low_j
    high i loses to   low j  <=>  status_high_i = 1 and t_low_j  > t_high_i

This module is intentionally dependency-free (no numpy import) so the parity
proof runs anywhere. Inputs may be Python lists or numpy arrays (any indexable).
"""
from __future__ import annotations
from typing import Sequence, Optional


# ---------------------------------------------------------------------------
# Naive dense reference (the correctness oracle, mirrors R kernel.py K=1)
# ---------------------------------------------------------------------------
def dense_pair_win_loss_1d(
    t_high: Sequence[float], s_high: Sequence[int], w_high: Sequence[float],
    t_low: Sequence[float], s_low: Sequence[int], w_low: Sequence[float],
) -> tuple[float, float, float]:
    """O(n_high * n_low) brute-force win/loss for one stratum pair (K=1)."""
    W = 0.0
    Lo = 0.0
    for i in range(len(t_high)):
        ti, si, wi = t_high[i], s_high[i], w_high[i]
        for j in range(len(t_low)):
            tj, sj, wj = t_low[j], s_low[j], w_low[j]
            if sj == 1 and ti > tj:
                W += wi * wj
            elif si == 1 and tj > ti:
                Lo += wi * wj
    tot = sum(w_high) * sum(w_low)
    return W, Lo, tot


# ---------------------------------------------------------------------------
# Fast single-pair sweep: O(n log n)
# ---------------------------------------------------------------------------
def fast_pair_win_loss_1d(
    t_high: Sequence[float], s_high: Sequence[int],
    t_low: Sequence[float], s_low: Sequence[int],
    w_high: Optional[Sequence[float]] = None,
    w_low: Optional[Sequence[float]] = None,
) -> tuple[float, float, float]:
    """Weighted win/loss for one stratum pair via a descending-time sweep."""
    n_high = len(t_high)
    n_low = len(t_low)
    if w_high is None and w_low is None:
        w_high = [1.0] * n_high
        w_low = [1.0] * n_low
    elif w_high is None or w_low is None:
        raise ValueError("weights must be supplied for both sides together")

    # entries: (time, side, status, weight); side 0 = high, 1 = low
    entries = [(t_high[i], 0, s_high[i], w_high[i]) for i in range(n_high)]
    entries += [(t_low[j], 1, s_low[j], w_low[j]) for j in range(n_low)]
    entries.sort(key=lambda e: e[0], reverse=True)  # descending time

    acc_high = 0.0
    acc_low = 0.0
    W = 0.0
    Lo = 0.0
    n = len(entries)
    i = 0
    while i < n:
        # one maximal group of equal time
        j = i
        t0 = entries[i][0]
        while j < n and entries[j][0] == t0:
            j += 1
        # score against STRICTLY greater times (accumulators from earlier groups)
        for k in range(i, j):
            _, side, status, w = entries[k]
            if side == 1 and status == 1:
                W += w * acc_high
            elif side == 0 and status == 1:
                Lo += w * acc_low
        # fold this tied group into the accumulators
        for k in range(i, j):
            _, side, _, w = entries[k]
            if side == 0:
                acc_high += w
            else:
                acc_low += w
        i = j

    tot = sum(w_high) * sum(w_low)
    return W, Lo, tot


# ---------------------------------------------------------------------------
# Fast all-adjacent-contrasts single pass: O(N log N) for every (d, d-1)
# ---------------------------------------------------------------------------
def sweep_all_adjacent_1d(
    time: Sequence[float], status: Sequence[int], stratum: Sequence[int],
    weight: Optional[Sequence[float]] = None, n_strata: Optional[int] = None,
    time_order: Optional[Sequence[int]] = None,
) -> dict[int, tuple[float, float, float]]:
    """
    Compute wins/losses/total for every adjacent contrast (d, d-1) in ONE pass.

    `time_order` (descending-time index order) is reusable across bootstrap
    iterations because the times are fixed; pass it to skip the per-iteration
    sort. Returns {d: (wins, losses, total)} keyed by the high stratum index d.
    """
    n = len(time)
    if weight is None:
        weight = [1.0] * n
    if n_strata is None:
        n_strata = max(stratum) if n else 0
    if time_order is None:
        time_order = sorted(range(n), key=lambda x: time[x], reverse=True)

    size = n_strata + 2  # allow indices 0 .. n_strata+1
    acc = [0.0] * size
    totw = [0.0] * size
    W = [0.0] * size
    Lo = [0.0] * size
    for x in range(n):
        totw[stratum[x]] += weight[x]

    i = 0
    m = len(time_order)
    while i < m:
        j = i
        t0 = time[time_order[i]]
        while j < m and time[time_order[j]] == t0:
            j += 1
        # score the equal-time group against strictly-greater-time accumulators
        for k in range(i, j):
            x = time_order[k]
            if status[x] == 1:
                s = stratum[x]
                W[s + 1] += weight[x] * acc[s + 1]   # x is LOW side of (s+1, s)
                Lo[s] += weight[x] * acc[s - 1]       # x is HIGH side of (s, s-1)
        # fold the group into per-stratum accumulators
        for k in range(i, j):
            x = time_order[k]
            acc[stratum[x]] += weight[x]
        i = j

    out: dict[int, tuple[float, float, float]] = {}
    for d in range(2, n_strata + 1):
        if totw[d] > 0 and totw[d - 1] > 0:
            out[d] = (W[d], Lo[d], totw[d] * totw[d - 1])
    return out


def dense_adjacent_win_loss(
    time: Sequence[float], status: Sequence[int], stratum: Sequence[int],
    weight: Optional[Sequence[float]] = None, n_strata: Optional[int] = None,
) -> dict[int, tuple[float, float, float]]:
    """O(N^2 / D) brute-force adjacent contrasts — the correctness oracle."""
    n = len(time)
    if weight is None:
        weight = [1.0] * n
    if n_strata is None:
        n_strata = max(stratum) if n else 0
    out: dict[int, tuple[float, float, float]] = {}
    for d in range(2, n_strata + 1):
        hi = [x for x in range(n) if stratum[x] == d]
        lo = [x for x in range(n) if stratum[x] == d - 1]
        if not hi or not lo:
            continue
        out[d] = dense_pair_win_loss_1d(
            [time[x] for x in hi], [status[x] for x in hi], [weight[x] for x in hi],
            [time[x] for x in lo], [status[x] for x in lo], [weight[x] for x in lo],
        )
    return out

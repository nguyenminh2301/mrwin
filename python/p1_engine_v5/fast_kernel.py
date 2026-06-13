"""
fast_kernel.py — Sub-quadratic reference implementation of the right-censored
hierarchical win/loss totals.

This is the executable reference (oracle) for the algorithm specified in
`inst/spec/fast-hierarchical-win-algorithm.md`. It computes EXACTLY the same
two-group weighted win/loss totals as the O(n_i n_j) kernel in
``kernel.py`` / ``kernel_sparse.py`` (and ``R/kernel.R``), but in
``O(3^K * N log^{K-1} N)`` time and ``O(N)`` memory.

Kernel semantics (i = high group, j = low group), priorities k = 0..K-1,
first deciding priority wins:

    WIN_k(i,j):  D_low[j,k] == 1 and T_high[i,k] >  T_low[j,k]   (i wins)
    LOSS_k(i,j): D_high[i,k] == 1 and T_low[j,k] >  T_high[i,k]  (i loses)
    otherwise undecided -> move to priority k+1.

Method (see spec Sections 5-6):

1. Telescope the first-decided-win indicator:
       Wins = sum_k ( prod_{m<k} (1 - WIN_m - LOSS_m) ) * WIN_k
2. Expand the prefix product. Each earlier priority m is ABSENT (+1),
   WIN (-1), or LOSS (-1). This yields (3^K - 1)/2 signed terms.
3. Each term is a weighted, strict, two-group multidimensional dominance
   count (after one-sided event-indicator filtering and per-coordinate sign
   flips that turn every comparison into the uniform form x_i > y_j).
4. Dominance counts are solved with a Fenwick tree (d<=2) and value-pivot
   CDQ divide-and-conquer (d>=3), all tie-safe via coordinate compression.

Exactness: for counts (unweighted), all intermediate sums are integers and
float64 is exact up to 2^53; at biobank scale the total pair count
(~1e11) is far below that bound, so the unweighted path returns bit-exact
integer totals.

Public API:
    fast_pair_win_loss(T_high, D_high, T_low, D_low, w_high=None, w_low=None)
    fast_pair_logcwr(...)            -> float
    fast_adjacent_win_loss(T, D, strata, weights=None)  -> list of dicts
"""
from __future__ import annotations

import itertools
import numpy as np

__all__ = [
    "fast_pair_win_loss",
    "fast_pair_logcwr",
    "fast_adjacent_win_loss",
]

# Prefix-priority states in the decomposition.
_ABSENT, _WIN, _LOSS = 0, 1, 2


# --------------------------------------------------------------------------
# Dominance engine
# --------------------------------------------------------------------------
class _Fenwick:
    """Fenwick / binary indexed tree over a compressed coordinate axis."""

    __slots__ = ("n", "t")

    def __init__(self, n: int):
        self.n = n
        self.t = np.zeros(n + 1, dtype=np.float64)

    def add(self, i: int, v: float) -> None:
        i += 1
        while i <= self.n:
            self.t[i] += v
            i += i & (-i)

    def prefix(self, i: int) -> float:
        """Sum over compressed ranks [0..i]."""
        i += 1
        s = 0.0
        while i > 0:
            s += self.t[i]
            i -= i & (-i)
        return s


def _compress(values: np.ndarray) -> dict:
    return {v: r for r, v in enumerate(np.unique(values))}


def _dom1(X, wx, Y, wy):
    """sum wx_i wy_j over X[i,0] > Y[j,0] (strict)."""
    order = np.argsort(X[:, 0], kind="mergesort")
    xs = X[order, 0]
    cw = np.concatenate([[0.0], np.cumsum(wx[order])])
    total = cw[-1]
    # count weight with xs <= y  -> cw[searchsorted(..., 'right')]
    idx = np.searchsorted(xs, Y[:, 0], side="right")
    return float(np.sum(wy * (total - cw[idx])))


def _dom2(X, wx, Y, wy):
    """Strict 2-D dominance: X[i] > Y[j] in both coords."""
    if len(X) == 0 or len(Y) == 0:
        return 0.0
    comp = _compress(np.concatenate([X[:, 1], Y[:, 1]]))
    fw = _Fenwick(len(comp))
    # Sweep dim 0 descending; within a tied dim-0 block, run all B-queries
    # before inserting tied A-points so equal dim-0 values never count.
    events = [(X[i, 0], 1, comp[X[i, 1]], float(wx[i])) for i in range(len(X))]
    events += [(Y[j, 0], 0, comp[Y[j, 1]], float(wy[j])) for j in range(len(Y))]
    events.sort(key=lambda e: (-e[0], e[1]))
    ans = 0.0
    inserted = 0.0
    for _, typ, rank, w in events:
        if typ == 0:  # B-query: inserted A with dim1 strictly greater
            ans += w * (inserted - fw.prefix(rank))
        else:  # A-insert
            fw.add(rank, w)
            inserted += w
    return ans


def _dom_cdq(X, wx, Y, wy):
    """General d-dim strict dominance via value-pivot CDQ on dim 0."""
    d = X.shape[1]
    if len(X) == 0 or len(Y) == 0:
        return 0.0
    if d == 1:
        return _dom1(X, wx, Y, wy)
    if d == 2:
        return _dom2(X, wx, Y, wy)
    vals = np.concatenate([X[:, 0], Y[:, 0]])
    distinct = np.unique(vals)
    if distinct.size == 1:
        return 0.0  # no strict domination on dim 0 is possible
    pivot = distinct[distinct.size // 2]
    if pivot == distinct[-1]:
        pivot = distinct[distinct.size // 2 - 1]
    xu = X[:, 0] > pivot  # upper A (strictly larger dim 0)
    yu = Y[:, 0] > pivot
    # cross: upper-A dominates lower-B on dim 0 (strict); residual dims 1..d-1
    ans = _dom_cdq(X[xu][:, 1:], wx[xu], Y[~yu][:, 1:], wy[~yu])
    ans += _dom_cdq(X[xu], wx[xu], Y[yu], wy[yu])
    ans += _dom_cdq(X[~xu], wx[~xu], Y[~yu], wy[~yu])
    return ans


def _dominance_sum(X, wx, Y, wy):
    if X.ndim != 2:
        raise ValueError("X must be 2-D")
    if X.shape[1] == 0:  # no conditions: full weighted cross product
        return float(wx.sum() * wy.sum())
    return _dom_cdq(X, wx, Y, wy)


# --------------------------------------------------------------------------
# Decomposition into signed dominance terms
# --------------------------------------------------------------------------
def _term_value(TH, DH, wH, TL, DL, wL, decide_k, decide_kind, prefix):
    """One decomposition term: deciding relation `decide_kind` at coord
    `decide_k`, with earlier coords fixed by `prefix` (coord -> WIN/LOSS).
    Returns the signed weighted dominance count.
    """
    conds = [(decide_k, decide_kind)] + list(prefix.items())
    imask = np.ones(len(TH), dtype=bool)
    jmask = np.ones(len(TL), dtype=bool)
    xcols, ycols = [], []
    for c, kind in conds:
        if kind == _WIN:  # D_low[j,c]=1 and T_high[i,c] > T_low[j,c]
            jmask &= DL[:, c] == 1
            xcols.append(TH[:, c])
            ycols.append(TL[:, c])
        else:  # LOSS: D_high[i,c]=1 and T_low[j,c] > T_high[i,c] -> negate
            imask &= DH[:, c] == 1
            xcols.append(-TH[:, c])
            ycols.append(-TL[:, c])

    X = np.column_stack(xcols)[imask]
    Y = np.column_stack(ycols)[jmask]
    val = _dominance_sum(X, wH[imask], Y, wL[jmask])
    return ((-1) ** len(prefix)) * val


def _decompose_total(TH, DH, wH, TL, DL, wL, decide_kind):
    K = TH.shape[1]
    total = 0.0
    for k in range(K):  # deciding priority
        for states in itertools.product((_ABSENT, _WIN, _LOSS), repeat=k):
            prefix = {m: s for m, s in enumerate(states) if s != _ABSENT}
            total += _term_value(TH, DH, wH, TL, DL, wL, k, decide_kind, prefix)
    return total


# --------------------------------------------------------------------------
# Public API
# --------------------------------------------------------------------------
def _prep(T, D):
    T = np.asarray(T, dtype=np.float64)
    D = np.asarray(D, dtype=np.int64)
    if T.ndim == 1:
        T = T[:, None]
        D = D[:, None]
    if T.shape != D.shape:
        raise ValueError("T and D must have the same shape")
    if not np.all(np.isfinite(T)):
        raise ValueError("T must be finite")
    if not np.isin(D, (0, 1)).all():
        raise ValueError("D must be binary 0/1")
    return T, D


def fast_pair_win_loss(T_high, D_high, T_low, D_low, w_high=None, w_low=None):
    """Weighted win/loss/total for one high-vs-low stratum pair.

    Returns ``(wins, losses, total)`` identical to
    ``kernel.stratum_win_loss`` on the corresponding dense block. ``wins`` and
    ``losses`` are returned as Python ``int`` when both weight vectors are
    ``None`` (bit-exact counts), else as ``float``.
    """
    TH, DH = _prep(T_high, D_high)
    TL, DL = _prep(T_low, D_low)
    if TH.shape[1] != TL.shape[1]:
        raise ValueError("high and low groups need the same number of priorities")
    nH, nL = TH.shape[0], TL.shape[0]

    unweighted = w_high is None and w_low is None
    wH = np.ones(nH) if w_high is None else np.asarray(w_high, dtype=np.float64)
    wL = np.ones(nL) if w_low is None else np.asarray(w_low, dtype=np.float64)
    if wH.shape != (nH,) or wL.shape != (nL,):
        raise ValueError("weight lengths must match the group sizes")
    if np.any(wH < 0) or np.any(wL < 0) or not np.all(np.isfinite(wH)) \
            or not np.all(np.isfinite(wL)):
        raise ValueError("weights must be finite and non-negative")

    wins = _decompose_total(TH, DH, wH, TL, DL, wL, _WIN)
    losses = _decompose_total(TH, DH, wH, TL, DL, wL, _LOSS)
    total = float(wH.sum() * wL.sum())

    if unweighted:
        return int(round(wins)), int(round(losses)), nH * nL
    return wins, losses, total


def fast_pair_logcwr(T_high, D_high, T_low, D_low,
                     w_high=None, w_low=None, floor=1e-12):
    """log(theta) = log(wins/losses) for one stratum pair, floored."""
    wins, losses, _ = fast_pair_win_loss(
        T_high, D_high, T_low, D_low, w_high, w_low
    )
    return float(np.log(max(wins, floor) / max(losses, floor)))


def fast_adjacent_win_loss(T, D, strata, weights=None):
    """Win/loss/total for every observed adjacent stratum pair (d, d-1).

    Mirrors ``R/kernel.R::mrwin_sparse_adjacent_win_loss``. Returns a list of
    dicts with keys ``high, low, wins, losses, total``.
    """
    T, D = _prep(T, D)
    strata = np.asarray(strata, dtype=np.int64)
    if strata.shape[0] != T.shape[0]:
        raise ValueError("strata length must match the number of rows")
    if weights is not None:
        weights = np.asarray(weights, dtype=np.float64)

    levels = np.unique(strata)
    out = []
    for d in levels[1:]:
        low = d - 1
        if low not in levels:
            continue
        hi = np.where(strata == d)[0]
        lo = np.where(strata == low)[0]
        wH = None if weights is None else weights[hi]
        wL = None if weights is None else weights[lo]
        wins, losses, total = fast_pair_win_loss(
            T[hi], D[hi], T[lo], D[lo], wH, wL
        )
        out.append({"high": int(d), "low": int(low),
                    "wins": wins, "losses": losses, "total": total})
    return out

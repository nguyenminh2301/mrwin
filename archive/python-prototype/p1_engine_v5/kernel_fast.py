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


# ===========================================================================
# C2: hierarchical K>=2 fast win/loss.
# ===========================================================================
def dense_pair_win_loss_kd(t_high, s_high, t_low, s_low, w_high=None, w_low=None):
    """Brute-force K-priority win/loss oracle (mirrors the dense kernel)."""
    nh = len(t_high)
    nl = len(t_low)
    wh = [1.0] * nh if w_high is None else list(w_high)
    wl = [1.0] * nl if w_low is None else list(w_low)
    K = len(t_high[0]) if nh else 0
    W = 0.0
    Lo = 0.0
    for i in range(nh):
        for j in range(nl):
            res = 0
            for k in range(K):
                ti = t_high[i][k]; tj = t_low[j][k]
                di = s_high[i][k]; dj = s_low[j][k]
                if dj == 1 and ti > tj:
                    res = 1
                    break
                if di == 1 and tj > ti:
                    res = -1
                    break
            if res == 1:
                W += wh[i] * wl[j]
            elif res == -1:
                Lo += wh[i] * wl[j]
    return W, Lo, sum(wh) * sum(wl)


class _BIT:
    """Fenwick tree for weighted prefix sums over 1..n."""
    __slots__ = ("n", "t")

    def __init__(self, n):
        self.n = n
        self.t = [0.0] * (n + 1)

    def add(self, i, v):
        while i <= self.n:
            self.t[i] += v
            i += i & -i

    def prefix(self, i):
        s = 0.0
        while i > 0:
            s += self.t[i]
            i -= i & -i
        return s


def _dom2d_count(ax, ay, aw, bx, by, bw, x_le, y_a_gt_b):
    """
    Sum_{i,j} aw_i * bw_j over A x B with:
      x:  ax_i <= bx_j   (x_le=True)   or   ax_i >= bx_j   (x_le=False)
      y:  ay_i >  by_j   (y_a_gt_b=True) or  ay_i <  by_j   (y_a_gt_b=False)
    Strict y, inclusive x. O((|A|+|B|) log(|A|+|B|)).
    """
    na = len(ax)
    nb = len(bx)
    if na == 0 or nb == 0:
        return 0.0
    # compress y over the union of ay and by
    ys = sorted(set(ay) | set(by))
    rank = {v: i + 1 for i, v in enumerate(ys)}  # 1-based
    bit = _BIT(len(ys))

    # event order over x: process A-inserts and B-queries by x. For x_le we sweep
    # ascending and insert A with ax<=bx before querying (A before B at equal x);
    # for x>= we sweep descending and insert A with ax>=bx before querying.
    a_items = [(ax[i], ay[i], aw[i]) for i in range(na)]
    b_items = [(bx[j], by[j], bw[j]) for j in range(nb)]
    if x_le:
        a_items.sort(key=lambda e: e[0])
        b_items.sort(key=lambda e: e[0])
    else:
        a_items.sort(key=lambda e: e[0], reverse=True)
        b_items.sort(key=lambda e: e[0], reverse=True)

    total_inserted = 0.0
    res = 0.0
    ai = 0
    na_ = len(a_items)
    for bx_j, by_j, bw_j in b_items:
        if x_le:
            # insert all A with ax <= bx_j
            while ai < na_ and a_items[ai][0] <= bx_j:
                _, ay_i, aw_i = a_items[ai]
                bit.add(rank[ay_i], aw_i)
                total_inserted += aw_i
                ai += 1
        else:
            # insert all A with ax >= bx_j (descending sweep)
            while ai < na_ and a_items[ai][0] >= bx_j:
                _, ay_i, aw_i = a_items[ai]
                bit.add(rank[ay_i], aw_i)
                total_inserted += aw_i
                ai += 1
        r = rank[by_j]
        if y_a_gt_b:
            # ay > by_j : total - (sum of ay <= by_j)
            s = total_inserted - bit.prefix(r)
        else:
            # ay < by_j : sum of ay <= (by_j - 1 rank)
            s = bit.prefix(r - 1)
        res += bw_j * s
    return res


def _col(rows, k):
    return [r[k] for r in rows]


def _sub(seq, idx):
    return [seq[i] for i in idx]


def fast_pair_win_loss_2d(t_high, s_high, t_low, s_low, w_high=None, w_low=None):
    """
    O(N log N) win/loss for K=2 priorities. t_*/s_* are sequences of length-2
    rows. Decomposes into level-1 separations (1D sweep) plus, for pairs tied at
    level 1, level-2 separations split over the four (status1_high, status1_low)
    regimes; the two mixed regimes use a 2D dominance count.
    """
    nh = len(t_high)
    nl = len(t_low)
    wh = [1.0] * nh if w_high is None else list(w_high)
    wl = [1.0] * nl if w_low is None else list(w_low)

    th1 = _col(t_high, 0); th2 = _col(t_high, 1)
    sh1 = _col(s_high, 0); sh2 = _col(s_high, 1)
    tl1 = _col(t_low, 0); tl2 = _col(t_low, 1)
    sl1 = _col(s_low, 0); sl2 = _col(s_low, 1)

    # Level-1 separations (the K=1 problem on column 1).
    W, Lo, tot = fast_pair_win_loss_1d(th1, sh1, tl1, sl1, wh, wl)

    # --- regime (high d1=0, low d1=0): always tie at 1 -> level-2 1D on subsets
    Hi = [i for i in range(nh) if sh1[i] == 0]
    Lj = [j for j in range(nl) if sl1[j] == 0]
    if Hi and Lj:
        w2, l2, _ = fast_pair_win_loss_1d(
            _sub(th2, Hi), _sub(sh2, Hi), _sub(tl2, Lj), _sub(sl2, Lj),
            _sub(wh, Hi), _sub(wl, Lj))
        W += w2; Lo += l2

    # --- regime (high d1=1, low d1=1): tie at 1 iff t1 equal -> group by t1
    Hi = [i for i in range(nh) if sh1[i] == 1]
    Lj = [j for j in range(nl) if sl1[j] == 1]
    if Hi and Lj:
        from collections import defaultdict
        gh = defaultdict(list); gl = defaultdict(list)
        for i in Hi:
            gh[th1[i]].append(i)
        for j in Lj:
            gl[tl1[j]].append(j)
        for v, his in gh.items():
            ljs = gl.get(v)
            if not ljs:
                continue
            w2, l2, _ = fast_pair_win_loss_1d(
                _sub(th2, his), _sub(sh2, his), _sub(tl2, ljs), _sub(sl2, ljs),
                _sub(wh, his), _sub(wl, ljs))
            W += w2; Lo += l2

    # --- regime (high d1=0, low d1=1): tie at 1 iff th1 <= tl1
    H = [i for i in range(nh) if sh1[i] == 0]
    # W: low needs d2=1, win at 2 (th2 > tl2)
    Bw = [j for j in range(nl) if sl1[j] == 1 and sl2[j] == 1]
    W += _dom2d_count(_sub(th1, H), _sub(th2, H), _sub(wh, H),
                      _sub(tl1, Bw), _sub(tl2, Bw), _sub(wl, Bw),
                      x_le=True, y_a_gt_b=True)
    # Lo: high needs d2=1, loss at 2 (tl2 > th2)
    Hl = [i for i in range(nh) if sh1[i] == 0 and sh2[i] == 1]
    Bl = [j for j in range(nl) if sl1[j] == 1]
    Lo += _dom2d_count(_sub(th1, Hl), _sub(th2, Hl), _sub(wh, Hl),
                       _sub(tl1, Bl), _sub(tl2, Bl), _sub(wl, Bl),
                       x_le=True, y_a_gt_b=False)

    # --- regime (high d1=1, low d1=0): tie at 1 iff th1 >= tl1
    H = [i for i in range(nh) if sh1[i] == 1]
    Bw = [j for j in range(nl) if sl1[j] == 0 and sl2[j] == 1]
    W += _dom2d_count(_sub(th1, H), _sub(th2, H), _sub(wh, H),
                      _sub(tl1, Bw), _sub(tl2, Bw), _sub(wl, Bw),
                      x_le=False, y_a_gt_b=True)
    Hl = [i for i in range(nh) if sh1[i] == 1 and sh2[i] == 1]
    Bl = [j for j in range(nl) if sl1[j] == 0]
    Lo += _dom2d_count(_sub(th1, Hl), _sub(th2, Hl), _sub(wh, Hl),
                       _sub(tl1, Bl), _sub(tl2, Bl), _sub(wl, Bl),
                       x_le=False, y_a_gt_b=False)

    return W, Lo, tot


# ---- K=3 fast path -------------------------------------------------------
class _BIT2D:
    """2D Fenwick tree for weighted point-add / prefix-sum over n2 x n3."""
    __slots__ = ("n2", "n3", "t")

    def __init__(self, n2, n3):
        self.n2 = n2
        self.n3 = n3
        self.t = [[0.0] * (n3 + 1) for _ in range(n2 + 1)]

    def add(self, i2, i3, w):
        a = i2
        while a <= self.n2:
            row = self.t[a]
            b = i3
            while b <= self.n3:
                row[b] += w
                b += b & -b
            a += a & -a

    def prefix(self, i2, i3):
        s = 0.0
        a = i2
        while a > 0:
            row = self.t[a]
            b = i3
            while b > 0:
                s += row[b]
                b -= b & -b
            a -= a & -a
        return s


def _dom1d_strict(ax, aw, bx, bw, gt):
    """Sum_{i,j} aw_i bw_j with a.x > b.x (gt) or a.x < b.x (!gt)."""
    import bisect
    if not ax or not bx:
        return 0.0
    order = sorted(range(len(ax)), key=lambda i: ax[i])
    xs = [ax[i] for i in order]
    cum = [0.0]
    for i in order:
        cum.append(cum[-1] + aw[i])
    total = cum[-1]
    res = 0.0
    for j in range(len(bx)):
        if gt:  # a.x > b.x : total - (a.x <= b.x)
            s = total - cum[bisect.bisect_right(xs, bx[j])]
        else:   # a.x < b.x
            s = cum[bisect.bisect_left(xs, bx[j])]
        res += bw[j] * s
    return res


def _dom3d_count(a1, a2, a3, aw, b1, b2, b3, bw, r1, r2, s3):
    """
    Sum_{i,j} aw_i bw_j with:
      dim1: a1_i <= b1_j (r1='le') or a1_i >= b1_j (r1='ge')
      dim2: a2_i <= b2_j (r2='le') or a2_i >= b2_j (r2='ge')
      dim3: a3_i >  b3_j (s3='gt') or a3_i <  b3_j (s3='lt')   [strict]
    CDQ divide-and-conquer on dim1, reducing each cross-merge to the validated
    2D dominance count on (dim2, dim3). O((|A|+|B|) log^2 N) time, O(N) memory.
    """
    na = len(a1)
    nb = len(b1)
    if na == 0 or nb == 0:
        return 0.0
    # events: (key1, is_b, x2, x3, w); A (is_b=0) before B (is_b=1) on dim1 ties
    ev = [(a1[i], 0, a2[i], a3[i], aw[i]) for i in range(na)]
    ev += [(b1[j], 1, b2[j], b3[j], bw[j]) for j in range(nb)]
    if r1 == 'le':
        ev.sort(key=lambda e: (e[0], e[1]))
    else:
        ev.sort(key=lambda e: (-e[0], e[1]))

    x_le = (r2 == 'le')
    y_gt = (s3 == 'gt')
    res = [0.0]

    def cdq(lo, hi):
        if hi - lo <= 1:
            return
        mid = (lo + hi) // 2
        cdq(lo, mid)
        cdq(mid, hi)
        ax2 = []; ax3 = []; aww = []
        for k in range(lo, mid):
            e = ev[k]
            if e[1] == 0:
                ax2.append(e[2]); ax3.append(e[3]); aww.append(e[4])
        if not ax2:
            return
        bx2 = []; bx3 = []; bww = []
        for k in range(mid, hi):
            e = ev[k]
            if e[1] == 1:
                bx2.append(e[2]); bx3.append(e[3]); bww.append(e[4])
        if not bx2:
            return
        res[0] += _dom2d_count(ax2, ax3, aww, bx2, bx3, bww,
                               x_le=x_le, y_a_gt_b=y_gt)

    cdq(0, na + nb)
    return res[0]


def _dom3d_count_dense(a1, a2, a3, aw, b1, b2, b3, bw, r1, r2, s3):
    """Reference O(N^2)-memory 2D-BIT implementation; correctness oracle only."""
    na = len(a1)
    nb = len(b1)
    if na == 0 or nb == 0:
        return 0.0
    ys2 = sorted(set(a2) | set(b2))
    ys3 = sorted(set(a3) | set(b3))
    rk2 = {v: i + 1 for i, v in enumerate(ys2)}
    rk3 = {v: i + 1 for i, v in enumerate(ys3)}
    n2 = len(ys2); n3 = len(ys3)
    bit = _BIT2D(n2, n3)

    a_order = sorted(range(na), key=lambda i: a1[i], reverse=(r1 == 'ge'))
    b_order = sorted(range(nb), key=lambda j: b1[j], reverse=(r1 == 'ge'))

    def rect(i2lo, i2hi, i3lo, i3hi):
        if i2lo > i2hi or i3lo > i3hi:
            return 0.0
        return (bit.prefix(i2hi, i3hi) - bit.prefix(i2lo - 1, i3hi)
                - bit.prefix(i2hi, i3lo - 1) + bit.prefix(i2lo - 1, i3lo - 1))

    res = 0.0
    ai = 0
    for j in b_order:
        b1j = b1[j]
        if r1 == 'le':
            while ai < na and a1[a_order[ai]] <= b1j:
                i = a_order[ai]; bit.add(rk2[a2[i]], rk3[a3[i]], aw[i]); ai += 1
        else:
            while ai < na and a1[a_order[ai]] >= b1j:
                i = a_order[ai]; bit.add(rk2[a2[i]], rk3[a3[i]], aw[i]); ai += 1
        r2b = rk2[b2[j]]; r3b = rk3[b3[j]]
        # dim2 range
        if r2 == 'le':
            i2lo, i2hi = 1, r2b
        else:
            i2lo, i2hi = r2b, n2
        # dim3 strict range
        if s3 == 'gt':
            i3lo, i3hi = r3b + 1, n3
        else:
            i3lo, i3hi = 1, r3b - 1
        res += bw[j] * rect(i2lo, i2hi, i3lo, i3hi)
    return res


def _tie_kind(di, dj):
    """Constraint on (t_high, t_low) for the pair to TIE at a level."""
    if di == 0 and dj == 0:
        return 'any'
    if di == 0 and dj == 1:
        return 'le'   # tie iff t_high <= t_low
    if di == 1 and dj == 0:
        return 'ge'   # tie iff t_high >= t_low
    return 'eq'        # both events: tie iff equal


def _level3_count(ax1, ax2, ax3, aw, bx1, bx2, bx3, bw, c1, c2, s3):
    """
    Sum over A x B of aw bw with tie-constraint c1 on dim1, c2 on dim2
    (each in {'any','le','ge','eq'}) and strict s3 on dim3 ('gt'/'lt').
    Reduces 'any' (drop) and 'eq' (group-by) dims, then dispatches to the
    1D/2D/3D counter.
    """
    if not ax1 or not bx1:
        return 0.0
    if c1 == 'eq':
        from collections import defaultdict
        ga = defaultdict(list); gb = defaultdict(list)
        for i in range(len(ax1)):
            ga[ax1[i]].append(i)
        for j in range(len(bx1)):
            gb[bx1[j]].append(j)
        tot = 0.0
        for v in ga.keys() & gb.keys():
            ia = ga[v]; jb = gb[v]
            tot += _level3_count(
                [ax1[i] for i in ia], [ax2[i] for i in ia], [ax3[i] for i in ia], [aw[i] for i in ia],
                [bx1[j] for j in jb], [bx2[j] for j in jb], [bx3[j] for j in jb], [bw[j] for j in jb],
                'any', c2, s3)
        return tot
    if c2 == 'eq':
        from collections import defaultdict
        ga = defaultdict(list); gb = defaultdict(list)
        for i in range(len(ax2)):
            ga[ax2[i]].append(i)
        for j in range(len(bx2)):
            gb[bx2[j]].append(j)
        tot = 0.0
        for v in ga.keys() & gb.keys():
            ia = ga[v]; jb = gb[v]
            tot += _level3_count(
                [ax1[i] for i in ia], [ax2[i] for i in ia], [ax3[i] for i in ia], [aw[i] for i in ia],
                [bx1[j] for j in jb], [bx2[j] for j in jb], [bx3[j] for j in jb], [bw[j] for j in jb],
                c1, 'any', s3)
        return tot
    # c1, c2 in {'any','le','ge'}
    active = [(d, c) for d, c in ((1, c1), (2, c2)) if c in ('le', 'ge')]
    if len(active) == 0:
        return _dom1d_strict(ax3, aw, bx3, bw, gt=(s3 == 'gt'))
    if len(active) == 1:
        d, c = active[0]
        ax = ax1 if d == 1 else ax2
        bx = bx1 if d == 1 else bx2
        return _dom2d_count(ax, ax3, aw, bx, bx3, bw,
                            x_le=(c == 'le'), y_a_gt_b=(s3 == 'gt'))
    return _dom3d_count(ax1, ax2, ax3, aw, bx1, bx2, bx3, bw, c1, c2, s3)


def fast_pair_win_loss_3d(t_high, s_high, t_low, s_low, w_high=None, w_low=None):
    """O(N log^2 N) win/loss for K=3 priorities. t_*/s_* are length-3 rows."""
    nh = len(t_high); nl = len(t_low)
    wh = [1.0] * nh if w_high is None else list(w_high)
    wl = [1.0] * nl if w_low is None else list(w_low)

    # terms k=1 and k=2 via the K=2 fast path on columns (1,2)
    th12 = [(r[0], r[1]) for r in t_high]; sh12 = [(r[0], r[1]) for r in s_high]
    tl12 = [(r[0], r[1]) for r in t_low]; sl12 = [(r[0], r[1]) for r in s_low]
    W, Lo, tot = fast_pair_win_loss_2d(th12, sh12, tl12, sl12, wh, wl)

    th1 = _col(t_high, 0); th2 = _col(t_high, 1); th3 = _col(t_high, 2)
    sh1 = _col(s_high, 0); sh2 = _col(s_high, 1); sh3 = _col(s_high, 2)
    tl1 = _col(t_low, 0); tl2 = _col(t_low, 1); tl3 = _col(t_low, 2)
    sl1 = _col(s_low, 0); sl2 = _col(s_low, 1); sl3 = _col(s_low, 2)

    # term k=3: tie at levels 1 and 2, separated at level 3, over 16 regimes
    for di1 in (0, 1):
        for di2 in (0, 1):
            Hbase = [i for i in range(nh) if sh1[i] == di1 and sh2[i] == di2]
            if not Hbase:
                continue
            for dj1 in (0, 1):
                for dj2 in (0, 1):
                    c1 = _tie_kind(di1, dj1)
                    c2 = _tie_kind(di2, dj2)
                    Lbase = [j for j in range(nl) if sl1[j] == dj1 and sl2[j] == dj2]
                    if not Lbase:
                        continue
                    # WIN at 3: low has event (sl3==1), th3 > tl3
                    Hw = Hbase
                    Lw = [j for j in Lbase if sl3[j] == 1]
                    if Hw and Lw:
                        W += _level3_count(
                            _sub(th1, Hw), _sub(th2, Hw), _sub(th3, Hw), _sub(wh, Hw),
                            _sub(tl1, Lw), _sub(tl2, Lw), _sub(tl3, Lw), _sub(wl, Lw),
                            c1, c2, 'gt')
                    # LOSS at 3: high has event (sh3==1), tl3 > th3 (th3 < tl3)
                    Hl = [i for i in Hbase if sh3[i] == 1]
                    Ll = Lbase
                    if Hl and Ll:
                        Lo += _level3_count(
                            _sub(th1, Hl), _sub(th2, Hl), _sub(th3, Hl), _sub(wh, Hl),
                            _sub(tl1, Ll), _sub(tl2, Ll), _sub(tl3, Ll), _sub(wl, Ll),
                            c1, c2, 'lt')
    return W, Lo, tot


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

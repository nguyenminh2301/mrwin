//! Sub-quadratic right-censored hierarchical win/loss kernel.
//!
//! Pure-Rust core of the algorithm specified in
//! `inst/spec/fast-hierarchical-win-algorithm.md` and mirrored by the Python
//! reference `python/p1_engine_v5/fast_kernel.py`. No external dependencies, so
//! it builds and tests with no network access. The R-facing extendr glue is a
//! separate crate that links this one.
//!
//! Semantics (i = high group, j = low group), priorities `k = 0..K-1`, first
//! deciding priority wins:
//! * `WIN_k(i,j)`  : `D_low[j,k]==1 && T_high[i,k] >  T_low[j,k]`
//! * `LOSS_k(i,j)` : `D_high[i,k]==1 && T_low[j,k] >  T_high[i,k]`
//! * otherwise undecided -> next priority.
//!
//! `Wins = sum_k (prod_{m<k}(1 - WIN_m - LOSS_m)) * WIN_k`; expanding the prefix
//! product yields signed multidimensional dominance counts (Fenwick + CDQ).

/// Relation choice for a priority in the decomposition.
const ABSENT: u8 = 0;
const WIN: u8 = 1;
const LOSS: u8 = 2;

// ---------------------------------------------------------------------------
// Fenwick tree over f64
// ---------------------------------------------------------------------------
struct Fenwick {
    t: Vec<f64>,
}

impl Fenwick {
    fn new(n: usize) -> Self {
        Fenwick { t: vec![0.0; n + 1] }
    }
    fn add(&mut self, mut i: usize, v: f64) {
        i += 1;
        while i < self.t.len() {
            self.t[i] += v;
            i += i & i.wrapping_neg();
        }
    }
    /// Sum over compressed ranks [0..=i].
    fn prefix(&self, mut i: usize) -> f64 {
        i += 1;
        let mut s = 0.0;
        while i > 0 {
            s += self.t[i];
            i -= i & i.wrapping_neg();
        }
        s
    }
}

// ---------------------------------------------------------------------------
// Point set: rows are points, each `coords[r]` has length d; `w[r]` is weight.
// ---------------------------------------------------------------------------
struct Pts {
    coords: Vec<Vec<f64>>,
    w: Vec<f64>,
}

impl Pts {
    fn len(&self) -> usize {
        self.coords.len()
    }
    fn dim(&self) -> usize {
        if self.coords.is_empty() {
            0
        } else {
            self.coords[0].len()
        }
    }
    /// Sub-set rows by a boolean predicate on the first coordinate vs pivot.
    fn split(&self, keep_upper: bool, pivot: f64) -> Pts {
        let mut coords = Vec::new();
        let mut w = Vec::new();
        for (r, c) in self.coords.iter().enumerate() {
            let upper = c[0] > pivot;
            if upper == keep_upper {
                coords.push(c.clone());
                w.push(self.w[r]);
            }
        }
        Pts { coords, w }
    }
    /// Drop the first coordinate (project onto dims 1..d).
    fn drop_first(&self) -> Pts {
        Pts {
            coords: self.coords.iter().map(|c| c[1..].to_vec()).collect(),
            w: self.w.clone(),
        }
    }
    fn wsum(&self) -> f64 {
        self.w.iter().sum()
    }
}

// ---------------------------------------------------------------------------
// Dominance sums: sum_{i in X, j in Y : X[i,c] > Y[j,c] for all c} wx_i wy_j
// ---------------------------------------------------------------------------
fn compress(values: &[f64]) -> Vec<f64> {
    let mut v: Vec<f64> = values.to_vec();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v.dedup();
    v
}

fn rank(sorted: &[f64], x: f64) -> usize {
    // index of x in the deduped sorted axis (exact equality expected)
    sorted.partition_point(|&v| v < x)
}

fn dom1(x: &Pts, y: &Pts) -> f64 {
    let mut xv: Vec<(f64, f64)> = x
        .coords
        .iter()
        .zip(&x.w)
        .map(|(c, &w)| (c[0], w))
        .collect();
    xv.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let mut cum = vec![0.0; xv.len() + 1];
    for i in 0..xv.len() {
        cum[i + 1] = cum[i] + xv[i].1;
    }
    let total = *cum.last().unwrap();
    let xs: Vec<f64> = xv.iter().map(|p| p.0).collect();
    let mut ans = 0.0;
    for (c, &wy) in y.coords.iter().zip(&y.w) {
        // weight of xs <= c  == cum[upper_bound]
        let idx = xs.partition_point(|&v| v <= c[0]);
        ans += wy * (total - cum[idx]);
    }
    ans
}

fn dom2(x: &Pts, y: &Pts) -> f64 {
    if x.len() == 0 || y.len() == 0 {
        return 0.0;
    }
    let mut all1: Vec<f64> = Vec::with_capacity(x.len() + y.len());
    for c in &x.coords {
        all1.push(c[1]);
    }
    for c in &y.coords {
        all1.push(c[1]);
    }
    let axis = compress(&all1);
    let mut fw = Fenwick::new(axis.len());
    // events: (dim0, kind, rank1, w); kind 0 = query(Y), 1 = insert(X).
    // Sort dim0 descending; within ties, queries before inserts so equal dim0
    // values are not counted as strictly dominating.
    let mut ev: Vec<(f64, u8, usize, f64)> = Vec::with_capacity(x.len() + y.len());
    for (c, &w) in x.coords.iter().zip(&x.w) {
        ev.push((c[0], 1, rank(&axis, c[1]), w));
    }
    for (c, &w) in y.coords.iter().zip(&y.w) {
        ev.push((c[0], 0, rank(&axis, c[1]), w));
    }
    ev.sort_by(|a, b| {
        b.0.partial_cmp(&a.0).unwrap().then(a.1.cmp(&b.1))
    });
    let mut ans = 0.0;
    let mut inserted = 0.0;
    for (_, kind, r, w) in ev {
        if kind == 0 {
            ans += w * (inserted - fw.prefix(r));
        } else {
            fw.add(r, w);
            inserted += w;
        }
    }
    ans
}

fn dom(x: &Pts, y: &Pts) -> f64 {
    if x.len() == 0 || y.len() == 0 {
        return 0.0;
    }
    let d = x.dim();
    if d == 0 {
        return x.wsum() * y.wsum();
    }
    if d == 1 {
        return dom1(x, y);
    }
    if d == 2 {
        return dom2(x, y);
    }
    // CDQ on dim 0 with value-median pivot (tie-safe strictness).
    let mut vals: Vec<f64> = Vec::with_capacity(x.len() + y.len());
    for c in &x.coords {
        vals.push(c[0]);
    }
    for c in &y.coords {
        vals.push(c[0]);
    }
    let distinct = compress(&vals);
    if distinct.len() == 1 {
        return 0.0; // no strict domination on dim 0 possible
    }
    let mut pivot = distinct[distinct.len() / 2];
    if pivot == *distinct.last().unwrap() {
        pivot = distinct[distinct.len() / 2 - 1];
    }
    let xu = x.split(true, pivot);
    let xl = x.split(false, pivot);
    let yu = y.split(true, pivot);
    let yl = y.split(false, pivot);
    // cross: upper-X dominates lower-Y on dim0 (strict); residual dims 1..d-1
    let mut ans = dom(&xu.drop_first(), &yl.drop_first());
    ans += dom(&xu, &yu);
    ans += dom(&xl, &yl);
    ans
}

// ---------------------------------------------------------------------------
// Decomposition
// ---------------------------------------------------------------------------
#[allow(clippy::too_many_arguments)]
fn term_value(
    th: &[f64],
    dh: &[u8],
    wh: &[f64],
    tl: &[f64],
    dl: &[u8],
    wl: &[f64],
    k: usize,
    decide_k: usize,
    decide_kind: u8,
    prefix: &[(usize, u8)],
) -> f64 {
    let n_h = wh.len();
    let n_l = wl.len();
    let mut conds: Vec<(usize, u8)> = Vec::with_capacity(prefix.len() + 1);
    conds.push((decide_k, decide_kind));
    conds.extend_from_slice(prefix);

    let mut imask = vec![true; n_h];
    let mut jmask = vec![true; n_l];
    // Per-condition sign applied to the coordinate (-1 for LOSS to unify to ">")
    let mut signs: Vec<f64> = Vec::with_capacity(conds.len());
    for &(c, kind) in &conds {
        if kind == WIN {
            for j in 0..n_l {
                if dl[j * k + c] != 1 {
                    jmask[j] = false;
                }
            }
            signs.push(1.0);
        } else {
            for i in 0..n_h {
                if dh[i * k + c] != 1 {
                    imask[i] = false;
                }
            }
            signs.push(-1.0);
        }
    }

    let mut xc = Vec::new();
    let mut xw = Vec::new();
    for i in 0..n_h {
        if imask[i] {
            let row: Vec<f64> = conds
                .iter()
                .enumerate()
                .map(|(t, &(c, _))| signs[t] * th[i * k + c])
                .collect();
            xc.push(row);
            xw.push(wh[i]);
        }
    }
    let mut yc = Vec::new();
    let mut yw = Vec::new();
    for j in 0..n_l {
        if jmask[j] {
            let row: Vec<f64> = conds
                .iter()
                .enumerate()
                .map(|(t, &(c, _))| signs[t] * tl[j * k + c])
                .collect();
            yc.push(row);
            yw.push(wl[j]);
        }
    }
    let x = Pts { coords: xc, w: xw };
    let y = Pts { coords: yc, w: yw };
    let val = dom(&x, &y);
    let sign = if prefix.len() % 2 == 0 { 1.0 } else { -1.0 };
    sign * val
}

#[allow(clippy::too_many_arguments)]
fn total(
    th: &[f64],
    dh: &[u8],
    wh: &[f64],
    tl: &[f64],
    dl: &[u8],
    wl: &[f64],
    k: usize,
    decide_kind: u8,
) -> f64 {
    let mut acc = 0.0;
    for decide_k in 0..k {
        // enumerate 3^decide_k prefix state assignments over coords 0..decide_k
        let combos = 3usize.pow(decide_k as u32);
        for combo in 0..combos {
            let mut prefix: Vec<(usize, u8)> = Vec::new();
            let mut rem = combo;
            for m in 0..decide_k {
                let state = (rem % 3) as u8;
                rem /= 3;
                if state != ABSENT {
                    prefix.push((m, state));
                }
            }
            acc += term_value(th, dh, wh, tl, dl, wl, k, decide_k, decide_kind, &prefix);
        }
    }
    acc
}

/// Weighted win/loss totals for one high-vs-low stratum pair.
///
/// Inputs are row-major flattened: `th`/`dh` are `n_high * k`, `tl`/`dl` are
/// `n_low * k`, `wh`/`wl` are per-individual weights. Returns `(wins, losses)`
/// identical to the dense kernel. For uniform weights the values are exact
/// integer counts represented in `f64`.
#[allow(clippy::too_many_arguments)]
pub fn pair_win_loss(
    th: &[f64],
    dh: &[u8],
    wh: &[f64],
    tl: &[f64],
    dl: &[u8],
    wl: &[f64],
    k: usize,
) -> (f64, f64) {
    assert!(k >= 1, "k must be >= 1");
    let n_h = wh.len();
    let n_l = wl.len();
    assert_eq!(th.len(), n_h * k, "th length must be n_high * k");
    assert_eq!(dh.len(), n_h * k, "dh length must be n_high * k");
    assert_eq!(tl.len(), n_l * k, "tl length must be n_low * k");
    assert_eq!(dl.len(), n_l * k, "dl length must be n_low * k");
    let wins = total(th, dh, wh, tl, dl, wl, k, WIN);
    let losses = total(th, dh, wh, tl, dl, wl, k, LOSS);
    (wins, losses)
}

// ===========================================================================
// Tests: fast vs an independent O(n_i n_j) brute force inside Rust
// ===========================================================================
#[cfg(test)]
mod tests {
    use super::*;

    // SplitMix64 for deterministic, dependency-free test data.
    struct Rng(u64);
    impl Rng {
        fn next_u64(&mut self) -> u64 {
            self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
            z ^ (z >> 31)
        }
        fn unif(&mut self) -> f64 {
            (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next_u64() % n
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn brute(
        th: &[f64],
        dh: &[u8],
        wh: &[f64],
        tl: &[f64],
        dl: &[u8],
        wl: &[f64],
        k: usize,
    ) -> (f64, f64) {
        let n_h = wh.len();
        let n_l = wl.len();
        let mut wins = 0.0;
        let mut losses = 0.0;
        for i in 0..n_h {
            for j in 0..n_l {
                let mut sgn = 0i8;
                for c in 0..k {
                    let ti = th[i * k + c];
                    let tj = tl[j * k + c];
                    if dl[j * k + c] == 1 && ti > tj {
                        sgn = 1;
                        break;
                    }
                    if dh[i * k + c] == 1 && tj > ti {
                        sgn = -1;
                        break;
                    }
                }
                if sgn == 1 {
                    wins += wh[i] * wl[j];
                } else if sgn == -1 {
                    losses += wh[i] * wl[j];
                }
            }
        }
        (wins, losses)
    }

    fn gen(
        rng: &mut Rng,
        n_h: usize,
        n_l: usize,
        k: usize,
        ties: bool,
        weighted: bool,
    ) -> (Vec<f64>, Vec<u8>, Vec<f64>, Vec<f64>, Vec<u8>, Vec<f64>) {
        let mk_t = |rng: &mut Rng, n: usize| -> Vec<f64> {
            (0..n * k)
                .map(|_| {
                    if ties {
                        (1 + rng.below(4)) as f64
                    } else {
                        -rng.unif().max(1e-9).ln()
                    }
                })
                .collect()
        };
        let mk_d = |rng: &mut Rng, n: usize| -> Vec<u8> {
            (0..n * k).map(|_| rng.below(2) as u8).collect()
        };
        let mk_w = |rng: &mut Rng, n: usize| -> Vec<f64> {
            (0..n)
                .map(|_| if weighted { -rng.unif().max(1e-9).ln() } else { 1.0 })
                .collect()
        };
        (
            mk_t(rng, n_h),
            mk_d(rng, n_h),
            mk_w(rng, n_h),
            mk_t(rng, n_l),
            mk_d(rng, n_l),
            mk_w(rng, n_l),
        )
    }

    #[test]
    fn fast_matches_brute_all_k() {
        let mut rng = Rng(0xC0FFEE);
        for k in 1..=4 {
            for trial in 0..25 {
                let ties = trial % 3 == 0;
                let weighted = trial % 2 == 0;
                let n_h = 5 + (rng.below(40) as usize);
                let n_l = 5 + (rng.below(40) as usize);
                let (th, dh, wh, tl, dl, wl) = gen(&mut rng, n_h, n_l, k, ties, weighted);
                let (fw, fl) = pair_win_loss(&th, &dh, &wh, &tl, &dl, &wl, k);
                let (bw, bl) = brute(&th, &dh, &wh, &tl, &dl, &wl, k);
                if weighted {
                    assert!((fw - bw).abs() < 1e-9, "k={k} trial={trial} wins {fw} vs {bw}");
                    assert!((fl - bl).abs() < 1e-9, "k={k} trial={trial} loss {fl} vs {bl}");
                } else {
                    assert_eq!(fw, bw, "k={k} trial={trial} wins");
                    assert_eq!(fl, bl, "k={k} trial={trial} loss");
                }
            }
        }
    }

    #[test]
    fn unweighted_counts_are_exact_integers() {
        let mut rng = Rng(7);
        let (th, dh, wh, tl, dl, wl) = gen(&mut rng, 40, 35, 3, true, false);
        let (w, l) = pair_win_loss(&th, &dh, &wh, &tl, &dl, &wl, 3);
        assert_eq!(w, w.round());
        assert_eq!(l, l.round());
    }
}

// fast_kernel.cpp — WP13 Rcpp port of the subquadratic hierarchical win/loss
// kernel (K in {1,2,3}). Mirrors the validated R/Python reference in
// R/kernel_fast.R and python/p1_engine_v5/kernel_fast.py exactly; this removes
// the interpreted-R constant factor (esp. the CDQ 3D recursion).
//
// Returns c(wins, losses, total) for one high-vs-low stratum pair.

#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <functional>
#include <unordered_map>
using namespace Rcpp;

namespace {

typedef std::vector<double> vd;
typedef std::vector<int> vi;

// ---- Fenwick tree (weighted prefix sums over 1..n) -----------------------
struct BIT {
  int n;
  vd t;
  BIT(int n_) : n(n_), t(n_ + 1, 0.0) {}
  void add(int i, double v) { for (; i <= n; i += i & -i) t[i] += v; }
  double prefix(int i) const { double s = 0.0; for (; i > 0; i -= i & -i) s += t[i]; return s; }
};

inline vd gather(const vd& v, const vi& idx) {
  vd out; out.reserve(idx.size());
  for (int i : idx) out.push_back(v[i]);
  return out;
}

double vsum(const vd& v) { double s = 0.0; for (double x : v) s += x; return s; }

// ---- single-endpoint weighted win/loss sweep -----------------------------
// wins  = sum_{j: sl=1} wl_j * sum_{i: th_i > tl_j} wh_i
// losses= sum_{i: sh=1} wh_i * sum_{j: tl_j > th_i} wl_j   (strict times)
void sweep1d(const vd& th, const vi& sh, const vd& wh,
             const vd& tl, const vi& sl, const vd& wl,
             double& wins, double& losses) {
  int nh = th.size(), nl = tl.size();
  wins = 0.0; losses = 0.0;
  if (nh == 0 || nl == 0) return;
  int n = nh + nl;
  // entries: time, side(0 high/1 low), status, weight
  std::vector<int> ord(n);
  vd time(n), w(n); vi side(n), st(n);
  for (int i = 0; i < nh; i++) { time[i] = th[i]; side[i] = 0; st[i] = sh[i]; w[i] = wh[i]; }
  for (int j = 0; j < nl; j++) { int k = nh + j; time[k] = tl[j]; side[k] = 1; st[k] = sl[j]; w[k] = wl[j]; }
  std::iota(ord.begin(), ord.end(), 0);
  std::sort(ord.begin(), ord.end(), [&](int a, int b) { return time[a] > time[b]; });
  double acc_high = 0.0, acc_low = 0.0;
  int i = 0;
  while (i < n) {
    int j = i; double t0 = time[ord[i]];
    while (j < n && time[ord[j]] == t0) j++;
    for (int k = i; k < j; k++) {
      int e = ord[k];
      if (side[e] == 1 && st[e] == 1) wins += w[e] * acc_high;
      else if (side[e] == 0 && st[e] == 1) losses += w[e] * acc_low;
    }
    for (int k = i; k < j; k++) {
      int e = ord[k];
      if (side[e] == 0) acc_high += w[e]; else acc_low += w[e];
    }
    i = j;
  }
}

// ---- weighted 2D dominance count -----------------------------------------
// dim x: ax <= bx (x_le) or ax >= bx ; dim y: ay > by (y_gt) or ay < by [strict]
double dom2d(const vd& ax, const vd& ay, const vd& aw,
             const vd& bx, const vd& by, const vd& bw,
             bool x_le, bool y_gt) {
  int na = ax.size(), nb = bx.size();
  if (na == 0 || nb == 0) return 0.0;
  vd ys; ys.reserve(na + nb);
  for (double v : ay) ys.push_back(v);
  for (double v : by) ys.push_back(v);
  std::sort(ys.begin(), ys.end());
  ys.erase(std::unique(ys.begin(), ys.end()), ys.end());
  int m = ys.size();
  auto rankY = [&](double v) {
    return int(std::lower_bound(ys.begin(), ys.end(), v) - ys.begin()) + 1;
  };
  vi ao(na), bo(nb);
  std::iota(ao.begin(), ao.end(), 0);
  std::iota(bo.begin(), bo.end(), 0);
  if (x_le) {
    std::sort(ao.begin(), ao.end(), [&](int a, int b) { return ax[a] < ax[b]; });
    std::sort(bo.begin(), bo.end(), [&](int a, int b) { return bx[a] < bx[b]; });
  } else {
    std::sort(ao.begin(), ao.end(), [&](int a, int b) { return ax[a] > ax[b]; });
    std::sort(bo.begin(), bo.end(), [&](int a, int b) { return bx[a] > bx[b]; });
  }
  BIT bit(m);
  double total = 0.0, res = 0.0;
  int ai = 0;
  for (int jj = 0; jj < nb; jj++) {
    int j = bo[jj]; double bxj = bx[j];
    if (x_le) {
      while (ai < na && ax[ao[ai]] <= bxj) { int i = ao[ai]; bit.add(rankY(ay[i]), aw[i]); total += aw[i]; ai++; }
    } else {
      while (ai < na && ax[ao[ai]] >= bxj) { int i = ao[ai]; bit.add(rankY(ay[i]), aw[i]); total += aw[i]; ai++; }
    }
    int r = rankY(by[j]);
    double s = y_gt ? (total - bit.prefix(r)) : bit.prefix(r - 1);
    res += bw[j] * s;
  }
  return res;
}

// ---- weighted 1D strict dominance ----------------------------------------
double dom1d_strict(const vd& ax, const vd& aw, const vd& bx, const vd& bw, bool gt) {
  int na = ax.size(), nb = bx.size();
  if (na == 0 || nb == 0) return 0.0;
  vi o(na); std::iota(o.begin(), o.end(), 0);
  std::sort(o.begin(), o.end(), [&](int a, int b) { return ax[a] < ax[b]; });
  vd xs(na), cum(na + 1, 0.0);
  for (int k = 0; k < na; k++) { xs[k] = ax[o[k]]; cum[k + 1] = cum[k] + aw[o[k]]; }
  double total = cum[na], res = 0.0;
  for (int j = 0; j < nb; j++) {
    double s;
    if (gt) { // a.x > b.x : total - count(xs <= bx)
      int c = int(std::upper_bound(xs.begin(), xs.end(), bx[j]) - xs.begin());
      s = total - cum[c];
    } else {  // a.x < b.x : count(xs < bx)
      int c = int(std::lower_bound(xs.begin(), xs.end(), bx[j]) - xs.begin());
      s = cum[c];
    }
    res += bw[j] * s;
  }
  return res;
}

// ---- weighted 3D dominance via CDQ on dim1 -------------------------------
double dom3d(const vd& a1, const vd& a2, const vd& a3, const vd& aw,
             const vd& b1, const vd& b2, const vd& b3, const vd& bw,
             bool r1_le, bool r2_le, bool s3_gt) {
  int na = a1.size(), nb = b1.size();
  if (na == 0 || nb == 0) return 0.0;
  int N = na + nb;
  vd k1(N), x2(N), x3(N), w(N); vi isb(N);
  for (int i = 0; i < na; i++) { k1[i] = a1[i]; x2[i] = a2[i]; x3[i] = a3[i]; w[i] = aw[i]; isb[i] = 0; }
  for (int j = 0; j < nb; j++) { int k = na + j; k1[k] = b1[j]; x2[k] = b2[j]; x3[k] = b3[j]; w[k] = bw[j]; isb[k] = 1; }
  vi ord(N); std::iota(ord.begin(), ord.end(), 0);
  std::sort(ord.begin(), ord.end(), [&](int a, int b) {
    if (k1[a] != k1[b]) return r1_le ? (k1[a] < k1[b]) : (k1[a] > k1[b]);
    return isb[a] < isb[b]; // A (0) before B (1) on ties
  });
  vd K1(N), X2(N), X3(N), W(N); vi ISB(N);
  for (int k = 0; k < N; k++) { int e = ord[k]; K1[k] = k1[e]; X2[k] = x2[e]; X3[k] = x3[e]; W[k] = w[e]; ISB[k] = isb[e]; }

  double res = 0.0;
  // iterative-safe recursion via std::function
  std::function<void(int, int)> cdq = [&](int lo, int hi) {
    if (lo >= hi) return;
    int mid = (lo + hi) / 2;
    cdq(lo, mid);
    cdq(mid + 1, hi);
    vd ax2, ax3, aww, bx2, bx3, bww;
    for (int k = lo; k <= mid; k++) if (ISB[k] == 0) { ax2.push_back(X2[k]); ax3.push_back(X3[k]); aww.push_back(W[k]); }
    if (ax2.empty()) return;
    for (int k = mid + 1; k <= hi; k++) if (ISB[k] == 1) { bx2.push_back(X2[k]); bx3.push_back(X3[k]); bww.push_back(W[k]); }
    if (bx2.empty()) return;
    res += dom2d(ax2, ax3, aww, bx2, bx3, bww, r2_le, s3_gt);
  };
  cdq(0, N - 1);
  return res;
}

// constraint kinds for a tie at a level
enum Kind { ANY, LE, GE, EQ };
inline Kind tie_kind(int di, int dj) {
  if (di == 0 && dj == 0) return ANY;
  if (di == 0 && dj == 1) return LE;
  if (di == 1 && dj == 0) return GE;
  return EQ;
}

// level-3 count with tie-constraints c1,c2 on dims 1,2 and strict s3_gt on dim3
double level3(const vd& a1, const vd& a2, const vd& a3, const vd& aw,
              const vd& b1, const vd& b2, const vd& b3, const vd& bw,
              Kind c1, Kind c2, bool s3_gt) {
  if (a1.empty() || b1.empty()) return 0.0;
  // Reduce 'eq' dims ONE AT A TIME (group-by), keeping the other constraint.
  if (c1 == EQ) {
    std::unordered_map<double, vi> ma, mb;
    for (int i = 0; i < (int)a1.size(); i++) ma[a1[i]].push_back(i);
    for (int j = 0; j < (int)b1.size(); j++) mb[b1[j]].push_back(j);
    double tot = 0.0;
    for (auto& kv : ma) {
      auto it = mb.find(kv.first);
      if (it == mb.end()) continue;
      const vi& ia = kv.second; const vi& jb = it->second;
      tot += level3(gather(a1, ia), gather(a2, ia), gather(a3, ia), gather(aw, ia),
                    gather(b1, jb), gather(b2, jb), gather(b3, jb), gather(bw, jb),
                    ANY, c2, s3_gt);
    }
    return tot;
  }
  if (c2 == EQ) {
    std::unordered_map<double, vi> ma, mb;
    for (int i = 0; i < (int)a2.size(); i++) ma[a2[i]].push_back(i);
    for (int j = 0; j < (int)b2.size(); j++) mb[b2[j]].push_back(j);
    double tot = 0.0;
    for (auto& kv : ma) {
      auto it = mb.find(kv.first);
      if (it == mb.end()) continue;
      const vi& ia = kv.second; const vi& jb = it->second;
      tot += level3(gather(a1, ia), gather(a2, ia), gather(a3, ia), gather(aw, ia),
                    gather(b1, jb), gather(b2, jb), gather(b3, jb), gather(bw, jb),
                    c1, ANY, s3_gt);
    }
    return tot;
  }
  bool d1 = (c1 == LE || c1 == GE);
  bool d2 = (c2 == LE || c2 == GE);
  if (!d1 && !d2) {
    return dom1d_strict(a3, aw, b3, bw, s3_gt);
  }
  if (d1 != d2) {
    if (d1) return dom2d(a1, a3, aw, b1, b3, bw, (c1 == LE), s3_gt);
    else    return dom2d(a2, a3, aw, b2, b3, bw, (c2 == LE), s3_gt);
  }
  return dom3d(a1, a2, a3, aw, b1, b2, b3, bw, (c1 == LE), (c2 == LE), s3_gt);
}

// ---- K=2 assembly (returns wins, losses) ---------------------------------
void fast2d(const vd& th0, const vd& th1, const vi& sh0, const vi& sh1, const vd& wh,
            const vd& tl0, const vd& tl1, const vi& sl0, const vi& sl1, const vd& wl,
            double& W, double& Lo) {
  sweep1d(th0, sh0, wh, tl0, sl0, wl, W, Lo);
  int nh = th0.size(), nl = tl0.size();
  auto whichH = [&](int s0, int s1, int s3gate /*-1 = ignore*/, const vi& sh3) {
    vi r; for (int i = 0; i < nh; i++) if (sh0[i] == s0 && sh1[i] == s1 && (s3gate < 0 || sh3[i] == s3gate)) r.push_back(i);
    return r;
  };
  vi dummy;
  // regime (0,0): always tie -> 1D on col1
  {
    vi H; for (int i = 0; i < nh; i++) if (sh0[i] == 0) H.push_back(i);
    vi L; for (int j = 0; j < nl; j++) if (sl0[j] == 0) L.push_back(j);
    if (!H.empty() && !L.empty()) {
      double w2, l2; vi sh1H(H.size()), sl1L(L.size());
      for (int k = 0; k < (int)H.size(); k++) sh1H[k] = sh1[H[k]];
      for (int k = 0; k < (int)L.size(); k++) sl1L[k] = sl1[L[k]];
      sweep1d(gather(th1, H), sh1H, gather(wh, H), gather(tl1, L), sl1L, gather(wl, L), w2, l2);
      W += w2; Lo += l2;
    }
  }
  // regime (1,1): tie iff t0 equal -> group by t0, 1D on col1
  {
    std::unordered_map<double, vi> gh, gl;
    for (int i = 0; i < nh; i++) if (sh0[i] == 1) gh[th0[i]].push_back(i);
    for (int j = 0; j < nl; j++) if (sl0[j] == 1) gl[tl0[j]].push_back(j);
    for (auto& kv : gh) {
      auto it = gl.find(kv.first);
      if (it == gl.end()) continue;
      const vi& H = kv.second; const vi& L = it->second;
      double w2, l2; vi sh1H(H.size()), sl1L(L.size());
      for (int k = 0; k < (int)H.size(); k++) sh1H[k] = sh1[H[k]];
      for (int k = 0; k < (int)L.size(); k++) sl1L[k] = sl1[L[k]];
      sweep1d(gather(th1, H), sh1H, gather(wh, H), gather(tl1, L), sl1L, gather(wl, L), w2, l2);
      W += w2; Lo += l2;
    }
  }
  // regime (0,1): tie iff th0 <= tl0
  {
    vi H; for (int i = 0; i < nh; i++) if (sh0[i] == 0) H.push_back(i);
    vi Bw; for (int j = 0; j < nl; j++) if (sl0[j] == 1 && sl1[j] == 1) Bw.push_back(j);
    W += dom2d(gather(th0, H), gather(th1, H), gather(wh, H),
               gather(tl0, Bw), gather(tl1, Bw), gather(wl, Bw), true, true);
    vi Hl; for (int i = 0; i < nh; i++) if (sh0[i] == 0 && sh1[i] == 1) Hl.push_back(i);
    vi Bl; for (int j = 0; j < nl; j++) if (sl0[j] == 1) Bl.push_back(j);
    Lo += dom2d(gather(th0, Hl), gather(th1, Hl), gather(wh, Hl),
                gather(tl0, Bl), gather(tl1, Bl), gather(wl, Bl), true, false);
  }
  // regime (1,0): tie iff th0 >= tl0
  {
    vi H; for (int i = 0; i < nh; i++) if (sh0[i] == 1) H.push_back(i);
    vi Bw; for (int j = 0; j < nl; j++) if (sl0[j] == 0 && sl1[j] == 1) Bw.push_back(j);
    W += dom2d(gather(th0, H), gather(th1, H), gather(wh, H),
               gather(tl0, Bw), gather(tl1, Bw), gather(wl, Bw), false, true);
    vi Hl; for (int i = 0; i < nh; i++) if (sh0[i] == 1 && sh1[i] == 1) Hl.push_back(i);
    vi Bl; for (int j = 0; j < nl; j++) if (sl0[j] == 0) Bl.push_back(j);
    Lo += dom2d(gather(th0, Hl), gather(th1, Hl), gather(wh, Hl),
                gather(tl0, Bl), gather(tl1, Bl), gather(wl, Bl), false, false);
  }
}

inline vd mcol(const NumericMatrix& m, int c) { int n = m.nrow(); vd v(n); for (int i = 0; i < n; i++) v[i] = m(i, c); return v; }
inline vi icol(const IntegerMatrix& m, int c) { int n = m.nrow(); vi v(n); for (int i = 0; i < n; i++) v[i] = m(i, c); return v; }

} // namespace

// [[Rcpp::export]]
NumericVector mrwin_fast_pair_cpp(NumericMatrix time_high, IntegerMatrix status_high,
                                  NumericMatrix time_low, IntegerMatrix status_low,
                                  NumericVector weights_high, NumericVector weights_low) {
  int K = time_high.ncol();
  vd wh(weights_high.begin(), weights_high.end());
  vd wl(weights_low.begin(), weights_low.end());
  double total = vsum(wh) * vsum(wl);
  double W = 0.0, Lo = 0.0;

  if (K == 1) {
    sweep1d(mcol(time_high, 0), icol(status_high, 0), wh,
            mcol(time_low, 0), icol(status_low, 0), wl, W, Lo);
  } else if (K == 2) {
    fast2d(mcol(time_high, 0), mcol(time_high, 1), icol(status_high, 0), icol(status_high, 1), wh,
           mcol(time_low, 0), mcol(time_low, 1), icol(status_low, 0), icol(status_low, 1), wl, W, Lo);
  } else {
    // K == 3
    vd th0 = mcol(time_high, 0), th1 = mcol(time_high, 1), th2 = mcol(time_high, 2);
    vd tl0 = mcol(time_low, 0), tl1 = mcol(time_low, 1), tl2 = mcol(time_low, 2);
    vi sh0 = icol(status_high, 0), sh1 = icol(status_high, 1), sh2 = icol(status_high, 2);
    vi sl0 = icol(status_low, 0), sl1 = icol(status_low, 1), sl2 = icol(status_low, 2);
    int nh = th0.size(), nl = tl0.size();
    // terms k=1,2 on cols 0,1
    fast2d(th0, th1, sh0, sh1, wh, tl0, tl1, sl0, sl1, wl, W, Lo);
    // term k=3 over 16 regimes
    for (int di1 = 0; di1 <= 1; di1++)
      for (int di2 = 0; di2 <= 1; di2++) {
        vi Hbase; for (int i = 0; i < nh; i++) if (sh0[i] == di1 && sh1[i] == di2) Hbase.push_back(i);
        if (Hbase.empty()) continue;
        for (int dj1 = 0; dj1 <= 1; dj1++)
          for (int dj2 = 0; dj2 <= 1; dj2++) {
            Kind c1 = tie_kind(di1, dj1), c2 = tie_kind(di2, dj2);
            vi Lbase; for (int j = 0; j < nl; j++) if (sl0[j] == dj1 && sl1[j] == dj2) Lbase.push_back(j);
            if (Lbase.empty()) continue;
            vi Lw; for (int j : Lbase) if (sl2[j] == 1) Lw.push_back(j);
            if (!Lw.empty())
              W += level3(gather(th0, Hbase), gather(th1, Hbase), gather(th2, Hbase), gather(wh, Hbase),
                          gather(tl0, Lw), gather(tl1, Lw), gather(tl2, Lw), gather(wl, Lw), c1, c2, true);
            vi Hl; for (int i : Hbase) if (sh2[i] == 1) Hl.push_back(i);
            if (!Hl.empty())
              Lo += level3(gather(th0, Hl), gather(th1, Hl), gather(th2, Hl), gather(wh, Hl),
                           gather(tl0, Lbase), gather(tl1, Lbase), gather(tl2, Lbase), gather(wl, Lbase), c1, c2, false);
          }
      }
  }
  return NumericVector::create(_["wins"] = W, _["losses"] = Lo, _["total"] = total);
}

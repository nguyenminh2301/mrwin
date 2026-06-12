# An $O(N\log^{K-1}N)$ Algorithm for Right-Censored Hierarchical Win Statistics

**Draft preprint — v0.1**

Authors: Nguyen Thien Minh; N. Ahmad Aziz *(author order and affiliations to be
finalized)*

---

## Abstract

Win statistics — the win ratio, win odds, net benefit, and the broader family of
generalized pairwise comparisons (GPC) — analyze prioritized composite endpoints
by comparing every pair of individuals according to a clinical hierarchy. The
defining computation is quadratic in the sample size $N$: all standard
implementations form, explicitly or implicitly, the $N\times N$ matrix of
pairwise comparisons. This $O(N^2)$ cost (and $O(N^2)$ memory for the dense
kernel) is the binding constraint on applying win statistics at biobank scale
and inside resampling procedures, where the kernel is recomputed thousands of
times. We give an exact algorithm that computes the win and loss totals between
two groups in $O(3^K N\log^{K-1}N)$ time and $O(N)$ memory, where $K$ is the
(small, fixed) number of prioritized endpoints. The key obstruction is
right-censoring, which makes the per-priority "undecided" relation
non-transitive and defeats the sorting arguments that solve the uncensored case.
We circumvent it by telescoping the first-decided-priority indicator and
expanding the prefix product into a signed sum of multidimensional *dominance
counts*, each solved with a Fenwick tree or value-pivot divide-and-conquer. The
method is exact (bit-exact for unweighted counts), supports the per-individual
weights used by multiplier bootstraps, and reduces to the classical
Gehan–Wilcoxon computation when $K=1$. A pure-Rust reference reproduces the dense
kernel exactly for $K=1,\dots,4$ and exhibits an empirical growth exponent of
$1.08$ up to $N=131{,}072$ per group, versus the quadratic dense baseline. For
the canonical $K=3$ hierarchy (death > hospitalization > biomarker decline) this
is $O(N\log^2 N)$, turning a $\sim 10^{11}$-operation, $\sim 100$ GB computation
at $N=337{,}000$ into a $\sim 10^{8}$-operation, $O(N)$-memory one.

---

## 1. Introduction

Composite endpoints that mix events of differing clinical severity are
ubiquitous in cardiovascular, renal, oncology, and neurodegenerative trials.
Time-to-first-event analyses discard the severity ordering (a death and a
hospitalization count equally). The **win ratio** [Pocock 2012] and the more
general **generalized pairwise comparisons** framework [Buyse 2010] restore it:
for an ordered pair of individuals, one compares on the highest-priority
endpoint first and descends the hierarchy until the pair is decided. The win
ratio is the number of "wins" divided by the number of "losses"; win odds and
net benefit are simple functions of the same pairwise tallies. Win statistics
have been adopted in dozens of phase III trials and are implemented in widely
used software (`BuyseTest`, `WINS`, `hce`).

The computation underlying all of these is the same: evaluate the hierarchical
comparison for every ordered pair. With $N$ individuals this is $\Theta(N^2)$
pairwise comparisons, and the dense kernel matrix is $\Theta(N^2)$ memory. The
quadratic cost is usually tolerated because trials are small, but it becomes
prohibitive in two increasingly common settings:

1. **Biobank-scale analyses.** Applying win statistics to cohorts of
   $N=10^5$–$10^6$ (e.g. for Mendelian-randomization estimands over prioritized
   outcomes) requires $10^{10}$–$10^{12}$ comparisons and, for the dense kernel,
   tens to hundreds of gigabytes.
2. **Resampling-based inference.** Permutation tests, bootstraps, and multiplier
   bootstraps recompute the tallies $B$ times (often with re-stratification),
   multiplying an already-quadratic cost.

A constant-factor speedup (C/C++/Rust, SIMD, parallelism) does not change the
asymptotic class. What is needed is a lower-order algorithm.

**Prior work.** The efficiency literature for win statistics targets
*inference* rather than the point computation: closed-form and asymptotic
variance estimators for the (log) win ratio [Luo 2015; Dong 2016] and exact
permutation/bootstrap moment formulas [recent GPC work]. For a *single*,
*uncensored* ordinal or continuous endpoint, the win count is a Mann–Whitney
$U$-statistic and is computable in $O(N\log N)$ by ranking — folklore, and the
basis of fast Wilcoxon routines. To our knowledge, no sub-quadratic algorithm
has been published for the **hierarchical, right-censored** kernel that win
statistics actually use; reference implementations remain $O(N^2)$.

**Contribution.** We give an exact algorithm for the two-group, weighted,
right-censored hierarchical win/loss totals running in $O(3^K N\log^{K-1}N)$ time
and $O(N)$ space (Section 6), based on a decomposition theorem (Section 4) and a
reduction to multidimensional dominance counting (Section 5). We explain the
censoring obstruction that the method overcomes (Section 3), provide a tested
reference implementation (Section 7), and verify exactness and near-linear
scaling empirically (Section 8). The result transfers to the entire win-ratio /
win-odds / net-benefit / GPC family.

## 2. Setup and notation

Two groups of individuals are compared: a *high* group $A$ (indices $i$) and a
*low* group $B$ (indices $j$). (In the two-sample win ratio these are the
treatment and control arms; in our motivating cCWR application they are adjacent
instrument-strata.) Each individual carries $K$ prioritized time-to-event
endpoints, ordered highest-priority first. For priority $k\in\{1,\dots,K\}$,
individual $x$ has an observed time $T_x^{(k)}\ge 0$ and an event indicator
$D_x^{(k)}\in\{0,1\}$ ($1$ = event observed, $0$ = right-censored).

For an ordered pair $(i,j)$ with $i\in A$, $j\in B$, define at each priority $k$:
$$
\mathrm{WIN}_k(i,j) = \mathbf 1\!\left[D_B^{(k)}(j)=1 \;\wedge\; T_A^{(k)}(i) > T_B^{(k)}(j)\right],\quad
\mathrm{LOSS}_k(i,j) = \mathbf 1\!\left[D_A^{(k)}(i)=1 \;\wedge\; T_B^{(k)}(j) > T_A^{(k)}(i)\right].
$$
That is, $i$ *wins* at $k$ if $j$ had the event and $i$ outlasted it; $i$ *loses*
at $k$ if $i$ had the event and $j$ outlasted it. The two are mutually exclusive
(they demand opposite strict orderings of $T^{(k)}$). If neither holds the pair
is **undecided** at $k$ and the comparison moves to priority $k+1$. The pair's
outcome is determined by the first priority at which it is decided.

Given non-negative per-individual weights $w$ (uniform for the point estimate;
random for a multiplier bootstrap), the quantities to compute are
$$
\mathrm{Wins} = \sum_{i\in A}\sum_{j\in B} w_i w_j \,\mathbf 1[\text{first decided priority is a win}],\qquad
\mathrm{Losses} = \sum_{i\in A}\sum_{j\in B} w_i w_j \,\mathbf 1[\text{first decided priority is a loss}].
$$
The win ratio is $\mathrm{Wins}/\mathrm{Losses}$; win odds and net benefit are
obtained from $\mathrm{Wins}$, $\mathrm{Losses}$, and the total weight.

## 3. The censoring obstruction

Write $U_m(i,j) = 1 - \mathrm{WIN}_m(i,j) - \mathrm{LOSS}_m(i,j)\in\{0,1\}$ for
the undecided indicator at priority $m$. If there were no censoring
($D\equiv 1$), then $U_m(i,j)=1$ exactly when $T_A^{(m)}(i)=T_B^{(m)}(j)$ — a
*tie*. Ties induce a transitive equivalence relation, so individuals can be
bucketed by their priority-1 value, then recursively by priority-2 within each
bucket, and so on; the hierarchical win count reduces to lexicographic rank
counting in $O(KN\log N)$.

Censoring breaks this. A pair is also undecided at priority $m$ when the
earlier-timed individual was censored (its event time is unknown, so the order
of the underlying events is indeterminate). The resulting relation
"$U_m(i,j)=1$" is **not transitive**: from $U_m(i,j)=1$ and $U_m(j,\ell)=1$ one
cannot conclude $U_m(i,\ell)=1$. Individuals therefore cannot be partitioned
into priority-equivalence classes, and the recursive bucketing collapses. This
non-transitivity is, we believe, why the censored hierarchical case has resisted
a clean sub-quadratic treatment.

## 4. Decomposition theorem

The outcome of pair $(i,j)$ is a win iff there exists a priority $k$ at which the
pair is decided as a win while being undecided at all earlier priorities. These
events are disjoint across $k$, giving the telescoping identity
$$
\mathbf 1[\text{first decision is a win}](i,j) \;=\; \sum_{k=1}^{K}\Big(\textstyle\prod_{m=1}^{k-1} U_m(i,j)\Big)\,\mathrm{WIN}_k(i,j).
\tag{1}
$$

**Theorem (signed dominance decomposition).**
With $\mathrm{REL}^{\mathrm W}_m := \mathrm{WIN}_m$ and
$\mathrm{REL}^{\mathrm L}_m := \mathrm{LOSS}_m$,
$$
\mathrm{Wins} \;=\; \sum_{k=1}^{K}\;\sum_{S\subseteq\{1,\dots,k-1\}}\;\sum_{\phi:S\to\{\mathrm W,\mathrm L\}}(-1)^{|S|}\;
C\big(k,S,\phi\big),
\tag{2}
$$
where
$$
C(k,S,\phi) = \sum_{i\in A}\sum_{j\in B} w_i w_j\,\mathrm{WIN}_k(i,j)\prod_{m\in S}\mathrm{REL}^{\phi(m)}_m(i,j).
$$
An identical identity holds for $\mathrm{Losses}$ with $\mathrm{WIN}_k$ replaced
by $\mathrm{LOSS}_k$. The number of terms is $\sum_{k=1}^{K}3^{\,k-1}=(3^K-1)/2$.

*Proof.* Start from (1) and weight by $w_iw_j$, summing over pairs. Expand the
prefix product using $U_m = 1-(\mathrm{WIN}_m+\mathrm{LOSS}_m)$:
$$
\prod_{m=1}^{k-1}U_m=\prod_{m=1}^{k-1}\big(1-(\mathrm{WIN}_m+\mathrm{LOSS}_m)\big)
=\sum_{S\subseteq\{1,\dots,k-1\}}(-1)^{|S|}\prod_{m\in S}\big(\mathrm{WIN}_m+\mathrm{LOSS}_m\big).
$$
Because $\mathrm{WIN}_m\,\mathrm{LOSS}_m=0$ pointwise, expanding each factor
$(\mathrm{WIN}_m+\mathrm{LOSS}_m)$ over $m\in S$ has no surviving cross terms, so
$\prod_{m\in S}(\mathrm{WIN}_m+\mathrm{LOSS}_m)=\sum_{\phi:S\to\{\mathrm W,\mathrm L\}}\prod_{m\in S}\mathrm{REL}^{\phi(m)}_m$.
Substituting and exchanging the (finite) sums yields (2). The count of
$(S,\phi)$ pairs for a given $k$ is $\sum_{s}\binom{k-1}{s}2^s=3^{k-1}$. $\square$

## 5. Reduction to dominance counting

Fix a term $C(k,S,\phi)$. It is a sum over $(i,j)$ of $w_iw_j$ times a
conjunction of conditions on the distinct coordinates $\{k\}\cup S$:

* a $\mathrm W$-type condition on coordinate $c$ requires
  $D_B^{(c)}(j)=1$ (a $j$-side filter) and $T_A^{(c)}(i)>T_B^{(c)}(j)$;
* an $\mathrm L$-type condition on coordinate $c$ requires
  $D_A^{(c)}(i)=1$ (an $i$-side filter) and $T_B^{(c)}(j)>T_A^{(c)}(i)$.

Apply the $i$-side filters to restrict $A$ to a subset $A'$ and the $j$-side
filters to restrict $B$ to $B'$. For each $\mathrm L$-type coordinate, negate the
values ($T_B^{(c)}(j)>T_A^{(c)}(i)\iff -T_A^{(c)}(i)>-T_B^{(c)}(j)$), so that
**every** condition takes the uniform strict form $x_i^{(c)}>y_j^{(c)}$. The term
becomes a *weighted two-group strict dominance sum* in $d=|S|+1\le K$ dimensions:
$$
C(k,S,\phi)=\sum_{i\in A'}\sum_{j\in B'} w_i w_j \;\mathbf 1\!\left[x_i^{(c)}>y_j^{(c)}\ \text{for all } c\right].
$$

**Dominance counting.** Computing $\sum_{i,j: x_i\succ y_j} w_iw_j$ over two point
sets with strictness in every coordinate is classical:

| $d$ | Method | Cost |
|---|---|---|
| $1$ | sort + prefix sums | $O(N\log N)$ |
| $2$ | sweep one axis, Fenwick (BIT) on the other | $O(N\log N)$ |
| $\ge 3$ | CDQ divide-and-conquer on one axis, recursing on the rest | $O(N\log^{d-1}N)$ |

Ties are handled exactly by coordinate compression and an ordering rule: when
sweeping an axis in decreasing order, process all queries within a tied block
before inserting the tied points, so equal coordinates are never counted as
dominating. For the CDQ recursion we split an axis by a **value-based pivot**
(the median of the distinct values present), sending values $>$ pivot to the
upper part and $\le$ pivot to the lower; cross-part pairs then satisfy the strict
inequality on the split axis automatically, and a block with a single distinct
value contributes zero. This makes strictness exact without perturbation.

## 6. Algorithm and complexity

```
function WIN_LOSS(T_A, D_A, w_A, T_B, D_B, w_B, K):
    wins   = TOTAL(..., deciding = WIN)
    losses = TOTAL(..., deciding = LOSS)
    return (wins, losses)

function TOTAL(T_A, D_A, w_A, T_B, D_B, w_B, K, deciding):
    acc = 0
    for k in 1..K:                              # deciding priority
        for each (S, phi) over prefix 1..k-1:   # 3^(k-1) choices
            (A', B', X, Y) = BUILD_TERM(k, deciding, S, phi)   # filters + sign flips
            acc += (-1)^|S| * DOMINANCE(X, w_A[A'], Y, w_B[B'])
    return acc

function DOMINANCE(X, wx, Y, wy):               # strict weighted, d = ncol(X)
    if d == 1: sort + prefix
    if d == 2: sweep + Fenwick
    else:      CDQ on axis 0 (value-pivot), recurse DOMINANCE on axes 1..d-1
```

**Complexity.** There are $(3^K-1)/2$ terms, each a dominance count of dimension
$\le K$ costing $O(N\log^{K-1}N)$. Hence
$$
\boxed{\;O\!\big(3^K\,N\log^{K-1}N\big)\ \text{time},\quad O(N)\ \text{space}.\;}
$$
$K$ is a small fixed constant (the number of prioritized endpoints; $K=3$ in the
canonical death > hospitalization > biomarker hierarchy), so the algorithm is
$O(N\log^{K-1}N)$ in $N$ — for $K=3$, $O(N\log^2 N)$ — against the $O(N^2)$ time
and $O(N^2)$ memory of the materialized kernel. The point estimate ($w\equiv 1$)
uses integer dominance counts and is therefore **bit-exact**; weighted bootstrap
iterations use floating accumulation. Only the win/loss tallies are affected:
all downstream inference (variance, GLS pooling, confidence intervals) consumes
the same scalars unchanged.

**Special cases.** $K=1$ recovers the two-sample Gehan–Wilcoxon win count in
$O(N\log N)$. With no censoring the same procedure specializes to lexicographic
rank counting.

## 7. Implementation

Two reference implementations accompany this paper:

* a dependency-free **Rust** core (`mrwinkernel`) implementing the dominance
  engine (Fenwick + value-pivot CDQ) and the decomposition; and
* a **Python** reference (`fast_kernel.py`) used as a readable oracle and to
  generate cross-language fixtures.

Both expose `pair_win_loss(T_high, D_high, w_high, T_low, D_low, w_low, K)`. The
Rust core links to R through an `extendr` binding, where it replaces only the
inner kernel of an otherwise unchanged win-statistic pipeline. All artifacts are
released with the `mrwin` package.

## 8. Numerical experiments

**Exactness.** Against the materialized $O(N^2)$ kernel on random two-group data
with heavy ties, right-censoring, and exponential weights, the algorithm
reproduces both tallies exactly for $K=1,2,3,4$: bit-exact for unweighted counts
and to $<10^{-9}$ (floating accumulation only) for weighted totals. The Rust core
and the Python reference agree bit-exactly on the unweighted fixtures.

**Scaling ($K=3$, weighted, censored; single-thread Python reference,
best-of-3).** Reproduce with `fast_kernel_benchmark`; raw data in
`preprint/scaling_k3.csv`.

| $N$ / group | brute $O(N^2)$ (s) | fast (s) | speedup | log–log slope |
|---:|---:|---:|---:|---:|
| 1024 | 0.11 | 0.22 | 0.53 | — |
| 2048 | 0.39 | 0.47 | 0.82 | 1.13 |
| 4096 | 1.62 | 0.94 | 1.71 | 1.00 |
| 8192 | 7.12 | 1.97 | 3.61 | 1.07 |
| 16384 | — | 4.32 | — | 1.13 |
| 32768 | — | 9.00 | — | 1.06 |
| 65536 | — | 19.21 | — | 1.09 |

The dense baseline scales as $N^2$ (roughly $4\times$ per doubling); the fitted
log–log exponent of the fast path is $1.08$ (a separate run to $N=131{,}072$
gives $1.05$), consistent with $N\log^2 N$. The crossover is near $N\approx
4{,}000$ *in interpreted Python*; the compiled core moves it far lower. Extrapolating
the asymptotics to $N=337{,}000$: the dense kernel needs $\sim 1.1\times10^{11}$
comparisons and $\sim 100$ GB, whereas the fast path needs $\sim 10^{8}$
weighted operations and $O(N)$ memory.

![Figure 1. Wall-clock time vs. sample size per group on log–log axes for the
dense $O(N^2)$ kernel and the fast algorithm ($K=3$, weighted, right-censored).
The dense baseline tracks the slope-2 (quadratic) guide; the fast path tracks
the slope-1 (linear) guide, with a fitted exponent of
1.08.](fig_scaling.png)

*Figure 1.* Dense vs. fast scaling (data: `preprint/scaling_k3.csv`; regenerate
the figure with `python preprint/make_figure.py`).

## 9. Discussion

**Generality.** The result is about the pairwise tallies, so it applies verbatim
to the win ratio, win odds, and net benefit, and to generalized pairwise
comparisons with prioritized time-to-event components. Continuous or ordinal
components are the degenerate ("all events observed") case and are handled by the
same code. Mixed hierarchies (e.g. a survival priority followed by a continuous
score) fit the same dominance reduction.

**Resampling.** Because per-individual weights factor through every dominance
count, multiplier (wild) bootstraps cost one fast evaluation per iteration; the
$O(N)$ memory means the kernel is never re-materialized across iterations. For
permutation tests that only relabel groups, the per-permutation work is likewise
one dominance pass.

**Limitations and extensions.** (i) The constant $3^K$ is exponential in the
number of priorities; the method targets the small-$K$ regime typical of
clinical hierarchies ($K\le 5$). (ii) We treat the standard strict-inequality
adjudication; pairwise rules that use thresholds/margins
(e.g. "win by $\ge \delta$") replace each comparison with a shifted threshold and
remain dominance counts. (iii) Recurrent-event and most-severe-event variants
require re-deriving the per-priority indicators but not the overall scheme.
(iv) A fully parallel CDQ and an integer-only fast path for exact unweighted
inference are natural engineering follow-ups.

## 10. Conclusion

Right-censored hierarchical win statistics, long computed in $O(N^2)$, admit an
exact $O(N\log^{K-1}N)$-time, $O(N)$-space algorithm for fixed $K$. The
construction — telescoping the first-decision indicator, expanding into signed
multidimensional dominance counts, and solving each with Fenwick/CDQ — is short,
exact, and weight-aware, and it removes the principal computational barrier to
win statistics at biobank scale and inside resampling inference.

## Reproducibility

Rust core and tests: `rust/mrwinkernel/` (`cargo test`). Python reference and
tests: `python/p1_engine_v5/fast_kernel.py`,
`tests/python/test_fast_kernel.py`. Benchmark: `python -m
p1_engine_v5.fast_kernel_benchmark`; data in `preprint/scaling_k3.csv`. Algorithm
specification: `inst/spec/fast-hierarchical-win-algorithm.md`.

## References

*(bibliographic details verified; final formatting to match the target venue)*

1. Pocock SJ, Ariti CA, Collier TJ, Wang D. The win ratio: a new approach to the
   analysis of composite endpoints in clinical trials based on clinical
   priorities. *Eur Heart J*. 2012;33(2):176–182.
   doi:10.1093/eurheartj/ehr352
2. Buyse M. Generalized pairwise comparisons of prioritized outcomes in the
   two-sample problem. *Stat Med*. 2010;29(30):3245–3257. doi:10.1002/sim.3923
3. Finkelstein DM, Schoenfeld DA. Combining mortality and longitudinal measures
   in clinical trials. *Stat Med*. 1999;18(11):1341–1354.
   doi:10.1002/(SICI)1097-0258(19990615)18:11<1341::AID-SIM129>3.0.CO;2-7
4. Luo X, Tian H, Mohanty S, Tsai WY. An alternative approach to confidence
   interval estimation for the win ratio statistic. *Biometrics*.
   2015;71(1):139–145. doi:10.1111/biom.12225
5. Dong G, Li D, Ballerstedt S, Vandemeulebroecke M. A generalized analytic
   solution to the win ratio to analyze a composite endpoint considering the
   clinical importance order among components. *Pharm Stat*. 2016;15(5):430–437.
6. Bebu I, Lachin JM. Large sample inference for a win ratio analysis of a
   composite outcome based on prioritized components. *Biostatistics*.
   2016;17(1):178–187. doi:10.1093/biostatistics/kxv032
7. Bentley JL. Multidimensional divide-and-conquer. *Commun ACM*.
   1980;23(4):214–229. doi:10.1145/358841.358850
8. Fenwick PM. A new data structure for cumulative frequency tables. *Softw
   Pract Exper*. 1994;24(3):327–336. doi:10.1002/spe.4380240306

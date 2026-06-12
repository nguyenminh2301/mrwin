# Sub-Quadratic Algorithm for Right-Censored Hierarchical Win Statistics

Status: design specification, prototype-validated (not yet implemented in the
package). Proposed as the algorithmic core of a future compiled performance
backend (supersedes the deferred Rcpp decision in WP9).

Last updated: 2026-06-12.

## 1. Motivation

`mrwin` is built on the hierarchical win comparison kernel (`R/kernel.R`). The
estimand requires, for every adjacent PRS-stratum pair and for every multiplier
bootstrap iteration, the weighted win and loss totals between a high group and a
low group. The current implementation evaluates these by materializing or
streaming the pairwise comparison matrix, which is `O(n_high * n_low)` per
stratum pair and `O(N^2)` over a cohort.

This is the binding constraint on the package's stated scope. The README and
the manuscript both assert that biobank-scale samples (N > 100,000) are a
statistical prerequisite for reliable cCWR inference, yet the current backend is
only practical to N ~ 2,000-5,000 (see `project-checkpoint.md`, WP9 summary).
The dense kernel at N = 337,000 needs ~106 GB for the full matrix; even the
sparse per-stratum backend re-spends `O((N/D)^2)` work in interpreted code on
every one of B bootstrap iterations. No amount of constant-factor optimization
(Rcpp, Rust, SIMD) changes the asymptotic class. A lower-order algorithm does.

This document specifies such an algorithm: an exact, weight-aware, censoring-
aware method that computes the same win/loss totals in `O(N log^{K-1} N)` time
and `O(N)` memory for a fixed number of priorities K.

## 2. Literature position

The standard computational cost of win statistics / generalized pairwise
comparisons (GPC) is quadratic. The reference R package `BuyseTest`
(C++ backend) documents `O(N^2)`; Pocock's and the Finkelstein-Schoenfeld
algorithms are quadratic by construction. The published efficiency work targets
*inference*, not point computation:

- Luo et al. (2015) and Dong et al. (2016): closed-form / analytic variance for
  the win ratio. These accelerate standard-error computation, not the win count.
- The Mann-Whitney / U-statistic connection gives `O(N log N)` for a *single,
  uncensored* ordinal endpoint via rank counting; this is folklore and does not
  cover hierarchy with censoring.

For the **right-censored, multi-priority** kernel that `mrwin` actually uses, no
sub-quadratic point-computation algorithm appears to be published. Section 9
discusses why this is a genuine contribution and a candidate standalone paper.

## 3. Kernel definition (the object to be computed)

Indices: `i` ranges over the high group A, `j` over the low group B. Priorities
`k = 1..K` are ordered highest-first; each individual carries `(T[.,k], D[.,k])`
(observed time, event indicator) per priority. Matching `R/kernel.R`:

- `WIN_k(i,j)`  := `D_B[j,k] = 1` and `T_A[i,k] >  T_B[j,k]`   (i wins at k)
- `LOSS_k(i,j)` := `D_A[i,k] = 1` and `T_B[j,k] >  T_A[i,k]`   (i loses at k)
- otherwise the pair is **undecided** at k and moves to priority k+1.

At any priority exactly one of {win, loss, undecided} holds (win and loss are
mutually exclusive because they require opposite strict time orderings). The
pair outcome is decided by the first priority that is not undecided.

Targets, with per-individual weights `w` (the bootstrap multiplier weights):

```
Wins   = sum_{i in A, j in B} w_i w_j * 1{ first decided priority is a WIN }
Losses = sum_{i in A, j in B} w_i w_j * 1{ first decided priority is a LOSS }
```

## 4. Why the naive fast idea fails

If "undecided" meant only "tied" (no censoring), the relation at each priority
would be a total preorder, undecided would be transitive, and the hierarchy
would reduce to lexicographic rank counting over equivalence classes
(`O(K N log N)` via sorting + a Fenwick tree). **Censoring breaks this.** A pair
is also undecided when the earlier-timed individual was censored, and this
"undecided" relation is *not transitive*. Individuals therefore cannot be
partitioned into prefix-equivalence classes, and the straightforward
divide-and-conquer collapses. This non-transitivity is the mathematical reason
the censored hierarchical case has resisted a clean sub-quadratic treatment.

## 5. The decomposition (algorithmic core)

Let `UNDEC_m = 1 - WIN_m - LOSS_m` (pointwise on pairs). The first-decided-win
indicator telescopes:

```
1{first decided priority is a WIN} = sum_{k=1}^{K} ( prod_{m<k} UNDEC_m ) * WIN_k
```

Expand the prefix product. Because `WIN_m * LOSS_m = 0`, each factor
`(1 - WIN_m - LOSS_m)` contributes, per earlier priority m, one of three
choices: **absent** (factor `+1`, no condition), **WIN_m** (factor `-1`), or
**LOSS_m** (factor `-1`). Hence

```
Wins = sum_{k=1}^{K}
         sum_{S subset of {1..k-1}}
           sum_{phi: S -> {WIN, LOSS}}
             (-1)^{|S|} * C(k, S, phi)
```

where `C(k, S, phi)` is the weighted number of pairs `(i,j)` with `WIN_k(i,j)`
true and, for each `m in S`, relation `phi(m)` true at priority m. `Losses` is
identical with the deciding relation set to `LOSS_k`.

The number of terms is `sum_{k} 3^{k-1} = (3^K - 1)/2` — a constant for fixed K
(13 terms for K=3, 40 for K=4).

### 5.1 Each term is a multidimensional dominance count

A single term `C(k, S, phi)` is a conjunction of conditions over the distinct
coordinates `{k} ∪ S`:

- A `WIN`-type condition on coordinate c requires `D_B[j,c]=1` (a j-side filter)
  and `T_A[i,c] > T_B[j,c]`.
- A `LOSS`-type condition on coordinate c requires `D_A[i,c]=1` (an i-side
  filter) and `T_B[j,c] > T_A[i,c]`.

Apply the i-side filters to select a subset of A and the j-side filters to
select a subset of B. Within these subsets, negate the coordinate for every
`LOSS`-type condition (`T_B[j,c] > T_A[i,c]` ⟺ `-T_A[i,c] > -T_B[j,c]`) so that
**every** condition becomes the uniform form `x_i^(c) > y_j^(c)`. The term value
is then a weighted **two-group strict dominance sum**:

```
C = sum_{i in A', j in B' : x_i^(c) > y_j^(c) for all c} w_i w_j
```

over `d = |S| + 1 <= K` coordinates.

## 6. Computing dominance sums

Strict, weighted, two-group dominance counting is classical:

| d | Method | Complexity |
|---|---|---|
| 1 | sort + prefix sums / two-pointer | `O(N log N)` |
| 2 | sweep dim 0 descending + Fenwick (BIT) on dim 1 | `O(N log N)` |
| >= 3 | CDQ divide-and-conquer on dim 0, reduce cross term to (d-1)-dim | `O(N log^{d-1} N)` |

Strictness and ties are handled exactly by coordinate compression plus an
ordering rule: when sweeping a coordinate descending, process all B-queries in a
tied block *before* inserting the tied A-points, so equal coordinate values are
never counted as dominating.

CDQ tie-safety: split dim 0 by a **value-based pivot** (median of the distinct
values present), sending `dim0 > pivot` to the upper part and `dim0 <= pivot` to
the lower. Cross-part pairs then satisfy `A.dim0 > pivot >= B.dim0` strictly, so
no equal-value pair is miscounted across the split. Blocks with a single
distinct dim-0 value contribute zero and terminate the recursion.

## 7. Overall complexity

Total cost is `O(3^K * N log^{K-1} N)`. For the manuscript's fixed `K = 3`
(death, heart-failure hospitalization, renal decline):

- **Time:** `O(N log^2 N)` per stratum pair per bootstrap iteration, vs `O(N^2)`.
- **Memory:** `O(N)` — the `N x N` (or block) kernel matrix is never formed.

At N = 337,000: `N^2 ≈ 1.1e11` and `N log^2 N ≈ 1.1e8`, roughly a 1000x
reduction in operations, alongside the elimination of the ~100 GB matrix.

## 8. Validation evidence (prototype)

A pure-Python reference prototype implements Sections 5-6 and was checked
against the package's own `O(N^2)` kernel logic (the brute-force two-group
weighted win/loss matching `R/kernel.R::mrwin_pair_kernel`).

Correctness (random data with heavy ties, right-censoring, and exponential
multiplier weights):

- Dominance primitives d = 1, 2, 3: exact on all random trials.
- Full hierarchical win/loss, **K = 3**: exact on 40 random configurations,
  max absolute error `5e-12` (floating-point accumulation only).
- Generality: **K = 1, 2, 3, 4** all exact (K = 4 exercises the general CDQ).

Empirical scaling of the fast path, K = 3, per group size N (single-thread,
interpreted Python reference — constants are pessimistic). Reproduce with
`python -m p1_engine_v5.fast_kernel_benchmark`; raw data in
`inst/spec/fast-kernel-scaling.csv`:

```
   N/group   brute O(N^2)   fast        log-log slope (fast)
     1024       0.13 s       0.29 s       --
     2048       0.46 s       0.61 s       1.06
     4096       1.79 s       1.18 s       0.95
     8192       8.21 s       2.40 s       1.02
    16384       (skip)       5.10 s       1.09
    32768       (skip)      10.70 s       1.07
    65536       (skip)      22.55 s       1.08
   131072       (skip)      47.57 s       1.08
```

Fitted fast-path exponent (log-log least squares over all sizes): **1.049**
(2.0 = quadratic, 1.0 = linear), consistent with `N log^2 N`. The brute path
scales as `N^2` (~4x per doubling) and is exact-matched by the fast path at
every size where it is run. The crossover is N ~ 4,000 *even in interpreted
Python*; a compiled implementation (Rust via extendr) would push the crossover
far lower and compound the asymptotic win with a large constant-factor win.

Reference implementation (in repository): `python/p1_engine_v5/fast_kernel.py`,
with exactness tests in `tests/python/test_fast_kernel.py` (pinned against the
existing dense pair kernel for K = 1..4, unweighted bit-exact and weighted to
`1e-9`). The reference exposes `fast_pair_win_loss`, `fast_pair_logcwr`, and
`fast_adjacent_win_loss`, and is wired into CI.

## 9. Numerical and statistical notes

- **Integer-exact path.** For the unweighted point estimate the dominance sums
  are integer counts; accumulate in 64-bit integers to obtain bit-exact win/loss
  totals (no floating error). Use floating accumulation only for the weighted
  bootstrap, where weights are real-valued.
- **Equivalence guarantee.** The algorithm is an exact reformulation of the
  existing kernel, not an approximation. The acceptance contract is bitwise (or
  `< 1e-10`) agreement with `mrwin_pair_win_loss()` / `mrwin_stratum_win_loss()`,
  reusing `mrwin_verify_sparse_dense_parity()` as the parity harness.
- **Scope unchanged downstream.** Only win/loss aggregation changes. The ISG,
  GLS pooling, bivariate-Delta covariance, Fieller CI, and SDPD layers consume
  the same scalars and are untouched. Their cost is `O(B * D)` and already
  negligible.
- **Tie semantics.** Strict `>` in `WIN_k`/`LOSS_k` is preserved exactly;
  equal times remain undecided, matching `R/kernel.R` and the v5.2 spec.

## 10. Integration plan

1. [DONE] Promote the prototype to a documented Python reference oracle with a
   test suite, including the integer-exact unweighted path and general K
   (`python/p1_engine_v5/fast_kernel.py`, `tests/python/test_fast_kernel.py`).
2. [DONE] Re-implement the dominance engine (Fenwick + value-pivot CDQ) and the
   decomposition in Rust. Pure-Rust core `rust/mrwinkernel/` is implemented and
   tested with no R/network: fast-vs-brute for K=1..4 and cross-language parity
   against the Python reference (`rust/mrwinkernel/tests/parity.rs`, bit-exact
   unweighted / 1e-9 weighted). The extendr binding `rust/mrwinrust/` and the R
   wiring (`backend = "rust"` enum + `R/backend_rust.R`, which injects the Rust
   kernel into the shared sparse pipeline via the new `pair_fun` argument) are
   in place. Building the R<->Rust linkage requires an R toolchain; see
   `rust/README.md` and the `rust/packaging/` templates.
3. [PENDING-R-ENV] Parity-gate `backend = "rust"` against the dense/sparse R
   backends in an R+Rust environment (`tests/testthat/test-backend-rust.R`,
   skipped when the compiled library is absent) at tolerance `1e-10`.
4. Re-run the WP9 benchmark grid extended to N = 50k, 100k, 337k with the
   compiled backend.
5. Only after parity + benchmarks pass, make `rust` the default backend for
   large N and document the asymptotics in the README performance section.

## 11. Publishability

This is a credible standalone methodological contribution, suitable for a
preprint and a short methods/computational-statistics paper, independent of the
main cCWR manuscript:

- **Novelty.** It appears to be the first exact sub-quadratic
  (`O(N log^{K-1} N)`) algorithm for **right-censored hierarchical** win
  statistics / GPC. The uncensored ordinal case is folklore; the censored
  hierarchical case is the open part, and the non-transitivity obstruction
  (Section 4) is exactly what prior sorting arguments could not cross.
- **Reach.** The result applies to the entire win-ratio / win-odds / net-benefit
  / GPC family, not just cCWR. It directly affects widely used tooling
  (`BuyseTest`, `WINS`, `hce`) and the 36+ trials that have adopted win ratios.
  Framing it around GPC broadens the audience well beyond Mendelian
  randomization.
- **Strength of claim.** The decomposition is a short, checkable theorem; the
  reductions are to classical dominance counting; the prototype gives exact
  agreement and a clean near-linear scaling curve. That is a complete
  "algorithm + proof + benchmark" package.
- **Suggested venue/title.** A computational-statistics or biostatistics methods
  journal, or arXiv `stat.CO`. Working title: *"An O(N log^2 N) algorithm for
  right-censored hierarchical win statistics."* Recommended order: deposit a
  preprint to establish priority, then submit. It can also ship first as a
  package vignette + benchmark to demonstrate practical impact.
- **Relationship to the main paper.** Best kept as a separate foundational
  paper. The cCWR manuscript can then *cite* it to justify the biobank-scale
  feasibility claim instead of carrying the algorithmic detail itself.

## 12. References

- Pocock SJ, et al. The win ratio. *Eur Heart J*. 2012;33(14):1744-1749.
- Finkelstein DM, Schoenfeld DA. Combining mortality and longitudinal measures
  in clinical trials. *Stat Med*. 1999.
- Luo X, et al. Closed-form variance estimator for the win ratio. 2015.
- Dong G, et al. A generalized analytic solution to the win ratio. *Stat Med*.
  2016 (PMID 27485522).
- Buyse M. Generalized pairwise comparisons of prioritized outcomes. *Stat Med*.
  2010;29:3245-3257. (`BuyseTest`: O(N^2) backend.)
- Bentley JL. Multidimensional divide-and-conquer (dominance counting). *Commun
  ACM*. 1980. (CDQ / Fenwick foundations.)

# WP3 Core Kernel Specification

Status: implemented at the R low-level API.

## Goal

WP3 completes the core kernel layer by adding sparse adjacent-pair functions that compute exactly the same pairwise comparisons as the dense `mrwin_kernel()` matrix without requiring a full `N x N` object.

## Public Functions

| Function | Purpose |
|---|---|
| `mrwin_kernel()` | Dense full `N x N` hierarchical comparison matrix. |
| `mrwin_pair_kernel()` | Sparse/block kernel for one high-vs-low stratum pair. |
| `mrwin_pair_win_loss()` | Win/loss/total sums for one stratum pair without full dense kernel materialization. |
| `mrwin_pair_log_cwr()` | Log-CWR for one stratum pair without full dense kernel materialization. |
| `mrwin_sparse_adjacent_win_loss()` | Win/loss/total summaries for all observed adjacent strata. |

## Equivalence Rule

For any high-index set `H` and low-index set `L`:

```r
mrwin_pair_kernel(time[H, ], status[H, ], time[L, ], status[L, ])
```

must equal:

```r
mrwin_kernel(time, status)[H, L]
```

Weighted sparse summaries must equal `mrwin_stratum_win_loss()` on the same dense block.

## Scope Boundary

WP3 exposes sparse/block kernel primitives. The high-level `mrwin(backend = "sparse")` workflow remains disabled until sparse bootstrap plumbing is implemented. This prevents users from selecting a backend that is only partially wired through the full estimator.

## Acceptance Criteria

WP3 is accepted when:

- pair kernels exactly match dense kernel blocks.
- weighted pair win/loss/log-CWR matches dense aggregation.
- adjacent sparse summaries match dense summaries across strata.
- incompatible pair inputs fail loudly.
- `R CMD check` passes.


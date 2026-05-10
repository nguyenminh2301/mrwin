\name{mrwin_precompute_pair_kernels}
\alias{mrwin_precompute_pair_kernels}
\title{Pre-compute Pair Kernels for Sparse Backend}
\description{
Computes the hierarchical comparison kernel for each adjacent stratum pair.
The pre-computed kernels can be reused across multiple calls with different weights.
}
\usage{
mrwin_precompute_pair_kernels(time, status, strata, contrast_plan)
}
\arguments{
\item{time}{Numeric matrix of endpoint times (N x K).}
\item{status}{Integer matrix of event indicators (N x K).}
\item{strata}{Integer vector of stratum assignments, length N.}
\item{contrast_plan}{Data frame with columns \code{high}, \code{low} defining adjacent contrasts.}
}
\value{
A list of length \code{nrow(contrast_plan)}, each element containing:
\item{high, low}{Stratum labels.}
\item{idx_high, idx_low}{Row indices for each stratum.}
\item{kernel}{Integer matrix of pairwise comparisons (-1, 0, +1).}
}
\seealso{
\code{\link{mrwin_sparse_estimate}}, \code{\link{mrwin_pair_kernel}}
}

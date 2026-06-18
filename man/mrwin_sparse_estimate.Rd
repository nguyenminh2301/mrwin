\name{mrwin_sparse_estimate}
\alias{mrwin_sparse_estimate}
\title{Sparse Backend Point Estimator}
\description{
Computes adjacent-stratum CWR, ISG, and GLS-pooling inputs using
the sparse backend (no full N x N kernel matrix).
}
\usage{
mrwin_sparse_estimate(
  time, status, G, X, beta_hat,
  n_strata = 10L, weights = NULL, floor = 1e-12,
  active_strata = NULL, pair_kernels = NULL, fast = FALSE
)
}
\arguments{
\item{time}{Numeric matrix of endpoint times (N x K).}
\item{status}{Integer matrix of event indicators (N x K), values 0/1.}
\item{G}{Numeric genotype matrix (N x M).}
\item{X}{Numeric exposure vector of length N.}
\item{beta_hat}{Numeric vector of GWAS weights, length M.}
\item{n_strata}{Integer, number of PRS strata (at least 2).}
\item{weights}{Optional numeric vector of observation weights, length N.}
\item{floor}{Numeric floor for log-ratio stability.}
\item{active_strata}{Optional integer vector of active stratum labels.}
\item{pair_kernels}{Optional pre-computed pair kernels from \code{mrwin_precompute_pair_kernels}.}

\item{fast}{Logical; if \code{TRUE}, use the compiled subquadratic win/loss kernel for K in \{1,2,3\}. Results are identical to the dense path.}
}
\value{
An object of class \code{mrwin_estimate} containing point estimates,
stratum assignments, win/loss counts, log-CWR, Delta-X, and ISG values.
}
\seealso{
\code{\link{mrwin_estimate}}, \code{\link{mrwin_precompute_pair_kernels}}
}

\name{mrwin_verify_sparse_dense_parity}
\alias{mrwin_verify_sparse_dense_parity}
\alias{print.mrwin_parity_check}
\title{Verify Sparse and Dense Backend Parity}
\description{
Runs both dense and sparse backends on the same simulated data
and checks that all numerical outputs match within tolerance.
}
\usage{
mrwin_verify_sparse_dense_parity(
  n = 200L, m_snps = 10L, n_strata = 4L,
  B = 20L, seed = 20260510L, tol = 1e-10
)

\method{print}{mrwin_parity_check}(x, ...)
}
\arguments{
\item{n}{Integer, sample size.}
\item{m_snps}{Integer, number of SNPs.}
\item{n_strata}{Integer, number of PRS strata.}
\item{B}{Integer, number of bootstrap iterations.}
\item{seed}{Integer random seed.}
\item{tol}{Numeric tolerance for comparison.}
\item{x}{A \code{mrwin_parity_check} object.}
\item{...}{Additional arguments (unused).}
}
\value{
An object of class \code{mrwin_parity_check} containing:
\item{checks}{Named list of maximum absolute differences.}
\item{all_close}{Logical, TRUE if all differences are below tolerance.}
\item{dense, sparse}{The two \code{mrwin_fit} objects.}
}
\seealso{
\code{\link{mrwin_benchmark}}, \code{\link{mrwin}}
}

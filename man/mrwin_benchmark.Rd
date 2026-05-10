\name{mrwin_benchmark}
\alias{mrwin_benchmark}
\alias{print.mrwin_benchmark}
\alias{summary.mrwin_benchmark}
\title{Performance Benchmark for mrwin Backends}
\description{
Runs the mrwin workflow across different sample sizes and backends,
recording runtime and numerical results.
}
\usage{
mrwin_benchmark(
  n_vec = c(200, 500, 1000), m_snps = 20L, n_strata = 5L,
  B = 50L, backends = c("dense", "sparse"),
  seed = 20260510L, n_iter = 1L
)

\method{print}{mrwin_benchmark}(x, ...)
\method{summary}{mrwin_benchmark}(object, ...)
}
\arguments{
\item{n_vec}{Integer vector of sample sizes to benchmark.}
\item{m_snps}{Integer, number of SNPs.}
\item{n_strata}{Integer, number of PRS strata.}
\item{B}{Integer, number of bootstrap iterations.}
\item{backends}{Character vector of backends to test: \code{"dense"}, \code{"sparse"}.}
\item{seed}{Integer random seed.}
\item{n_iter}{Integer, repetitions per configuration.}
\item{x, object}{A \code{mrwin_benchmark} object.}
\item{...}{Additional arguments (unused).}
}
\value{
An object of class \code{mrwin_benchmark} containing:
\item{results}{Data frame with columns: n, m_snps, n_strata, B, backend, iteration, success, time_seconds, peak_memory_mb, delta_gls, dscwr, n_valid_boot.}
\item{n_vec, m_snps, n_strata, B, backends, seed, n_iter}{Configuration metadata.}
}
\seealso{
\code{\link{mrwin}}, \code{\link{mrwin_verify_sparse_dense_parity}}
}

\name{mrwin_sparse_bootstrap}
\alias{mrwin_sparse_bootstrap}
\title{Sparse Backend Multiplier Bootstrap}
\description{
Runs the multiplier bootstrap using the sparse backend.
Computes per-pair kernels once, then re-applies bootstrap weights each iteration.
}
\usage{
mrwin_sparse_bootstrap(
  time, status, G, X, beta_hat, sigma_beta = 0,
  n_strata = 10L, B = 200L, seed = NULL,
  beta_draws = NULL, multiplier_weights = NULL,
  floor = 1e-12, covariates = NULL,
  adjustment = c("none", "ordinal_iptw"),
  iptw_truncation = c(0.01, 0.99), ess_fraction = 0.5
)
}
\arguments{
\item{time}{Numeric matrix of endpoint times (N x K).}
\item{status}{Integer matrix of event indicators (N x K).}
\item{G}{Numeric genotype matrix (N x M).}
\item{X}{Numeric exposure vector of length N.}
\item{beta_hat}{Numeric vector of GWAS weights, length M.}
\item{sigma_beta}{Numeric vector of GWAS standard errors, length M or 1.}
\item{n_strata}{Integer, number of PRS strata.}
\item{B}{Integer, number of bootstrap iterations (at least 2).}
\item{seed}{Optional random seed for reproducibility.}
\item{beta_draws}{Optional B x M matrix of GWAS weight draws for testing.}
\item{multiplier_weights}{Optional B x N matrix of multiplier weights for testing.}
\item{floor}{Numeric floor for log-ratio stability.}
\item{covariates}{Optional covariate matrix for IPTW adjustment.}
\item{adjustment}{Adjustment method: \code{"none"} or \code{"ordinal_iptw"}.}
\item{iptw_truncation}{Length-2 numeric vector of truncation quantiles.}
\item{ess_fraction}{Numeric ESS fraction threshold.}
}
\value{
An object of class \code{mrwin_bootstrap} containing bootstrap inference
results, point estimates, covariance, Fieller CI, Q statistic, and diagnostics.
}
\seealso{
\code{\link{mrwin_multiplier_bootstrap}}, \code{\link{mrwin}}
}

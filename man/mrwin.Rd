\name{mrwin}
\alias{mrwin}
\alias{mrwin_endpoint}
\alias{mrwin_gwas}
\alias{mrwin_controls}
\alias{mrwin_validate_data}
\alias{mrwin_config}
\alias{mrwin_simulate}
\alias{mrwin_kernel}
\alias{mrwin_pair_kernel}
\alias{mrwin_pair_win_loss}
\alias{mrwin_pair_log_cwr}
\alias{mrwin_sparse_adjacent_win_loss}
\alias{mrwin_stratum_win_loss}
\alias{mrwin_stratum_log_cwr}
\alias{mrwin_prs_strata}
\alias{mrwin_estimate}
\alias{mrwin_gls_pool}
\alias{mrwin_multiplier_bootstrap}
\alias{mrwin_aalen_per_snp}
\alias{mrwin_cox_per_snp}
\alias{mrwin_mr_egger}
\alias{mrwin_pleiotropy_bounded_ci}
\title{Causal win statistics for hierarchical composite endpoints}
\description{
R implementation of causal win-statistic methods for hierarchical composite endpoints.
The high-level \code{mrwin()} workflow validates endpoint, genotype, exposure, and
GWAS inputs before estimating adjacent-stratum causal win ratios, dose-standardized
GLS summaries, bootstrap intervals, heterogeneity diagnostics, and SDPD diagnostics.
}
\details{
The package is under active development. The current user-facing API is centered on
\code{mrwin()}, with lower-level functions exported for testing and method development.
See \code{inst/spec/algorithm-spec.md} for the implementation specification.
}
\value{
Return values depend on the function. \code{mrwin()} returns an object of class
\code{mrwin_fit}; \code{mrwin_validate_data()} returns an object of class
\code{mrwin_validated_data}; lower-level functions return lists, matrices, or numeric
vectors as appropriate.
}
\keyword{package}

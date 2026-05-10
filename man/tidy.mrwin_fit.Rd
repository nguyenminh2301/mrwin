\name{tidy.mrwin_fit}
\alias{tidy.mrwin_fit}
\title{Tidy an mrwin_fit Object}
\description{
Converts an \code{mrwin_fit} object into a tidy data frame with one row,
compatible with the broom tidy() convention.
}
\usage{
\method{tidy}{mrwin_fit}(x, ...)
}
\arguments{
\item{x}{An \code{mrwin_fit} object from \code{\link{mrwin}}.}
\item{...}{Additional arguments (unused).}
}
\value{
A data frame with columns:
\item{term}{Character, always \code{"DS-CWR"}.}
\item{estimate}{DS-CWR estimate.}
\item{delta}{Log-scale GLS estimate.}
\item{se_delta}{Standard error of delta_GLS.}
\item{statistic}{Wald statistic.}
\item{p_value}{Wald p-value.}
\item{ci_low, ci_high}{95\% CI on DS-CWR scale.}
\item{ci_delta_low, ci_delta_high}{95\% CI on log scale.}
\item{q_stat, q_p_value}{Heterogeneity Q statistic.}
\item{n, m_snps, n_strata, n_valid_bootstrap}{Data and bootstrap metadata.}
}
\seealso{
\code{\link{mrwin}}, \code{\link{mrwin_report}}
}

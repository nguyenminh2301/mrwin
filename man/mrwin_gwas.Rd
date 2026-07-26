\name{mrwin_gwas}
\alias{mrwin_gwas}
\title{Specify GWAS instrument weights}
\description{
Assemble the per-SNP GWAS effect sizes (and their uncertainty) used to build the
polygenic score that instruments the exposure in \code{\link{mrwin}}. Supplying
standard errors (or a full covariance) lets the estimator propagate
GWAS-weight uncertainty into the confidence interval.
}
\usage{
mrwin_gwas(beta, se = NULL, snp = NULL, covariance = NULL)
}
\arguments{
\item{beta}{Numeric vector of length \eqn{M}: the GWAS effect of each SNP on the
exposure (the instrument weights).}
\item{se}{Optional numeric vector of standard errors for \code{beta}. Defaults
to zeros (treat the weights as fixed/known). Must be finite and non-negative.}
\item{snp}{Optional character vector of SNP identifiers. Defaults to
\code{"snp_1"}, \code{"snp_2"}, \dots}
\item{covariance}{Optional \eqn{M \times M} covariance matrix for \code{beta}
(use when the weights are correlated, e.g.\ from LD). If given and \code{se} is
omitted, \code{se} is taken from its diagonal.}
}
\value{
An object of class \code{"mrwin_gwas_spec"}: a list with elements \code{beta},
\code{se}, \code{snp} and \code{covariance}, suitable for the \code{gwas}
argument of \code{\link{mrwin}}.
}
\seealso{\code{\link{mrwin}}, \code{\link{mrwin_endpoint}}, \code{\link{mrwin_gwas_resample_covariance}}}
\examples{
# 10 SNPs, weights known up to a standard error
g <- mrwin_gwas(beta = rnorm(10, 0, 0.05), se = rep(0.01, 10))
str(g)
}

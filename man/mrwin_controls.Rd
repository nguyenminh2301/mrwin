\name{mrwin_controls}
\alias{mrwin_controls}
\title{Estimation, inference and diagnostic settings for \code{mrwin}}
\description{
Collect the tuning knobs for a fit in one object: how many polygenic-score
strata to form, how to compute confidence intervals, which compiled backend to
use, whether to run pleiotropy diagnostics, and optional covariate adjustment.
Pass the result as the \code{controls} argument of \code{\link{mrwin}}.
}
\usage{
mrwin_controls(n_strata = 10L, bootstrap = 200L, seed = NULL,
               block_size = 4000L, run_sdpd = TRUE,
               sdpd_scale = c("aalen", "cox", "both", "none"),
               sdpd_alpha = 0.05, sdpd_min_snps = 10L,
               pleiotropy_bias_radius = NULL,
               adjustment = c("none", "ordinal_iptw", "gps"),
               backend = c("dense", "sparse", "fast", "rcpp"),
               inference = c("bootstrap", "analytic"),
               stratification = c("prs_rank", "doubly_ranked"),
               delta_x_tol = 1e-8, iptw_truncation = c(0.01, 0.99),
               ess_fraction = 0.5)
}
\arguments{
\item{n_strata}{Number of polygenic-score strata \eqn{D}; the estimator pools
\eqn{D-1} adjacent-stratum gradients. Report a sensitivity curve over a few
values rather than a single \code{n_strata}.}
\item{bootstrap}{Number of multiplier-bootstrap resamples (used when
\code{inference = "bootstrap"}).}
\item{seed}{Optional integer seed for reproducible resampling.}
\item{block_size}{Block size for the dense kernel's blocked evaluation.}
\item{run_sdpd}{Logical; run the step-down pleiotropy diagnostic (SDPD).}
\item{sdpd_scale}{Per-SNP effect scale for SDPD: Aalen, Cox, both, or none.}
\item{sdpd_alpha}{Significance level for the SDPD step-down.}
\item{sdpd_min_snps}{Minimum SNPs required to run SDPD.}
\item{pleiotropy_bias_radius}{Optional bias radius for the pleiotropy-bounded
confidence interval (\code{\link{mrwin_pleiotropy_bounded_ci}}).}
\item{adjustment}{Covariate adjustment: \code{"none"}, ordinal
inverse-probability-of-treatment weighting, or generalized propensity score.}
\item{backend}{Win/loss engine: \code{"dense"} (exact reference),
\code{"fast"} (subquadratic compiled, biobank-scale), \code{"sparse"}, or
\code{"rcpp"}.}
\item{inference}{\code{"bootstrap"} (default) or \code{"analytic"}
influence-function variance.}
\item{stratification}{\code{"prs_rank"} (rank on the polygenic score) or
\code{"doubly_ranked"} (see \code{\link{mrwin_doubly_ranked_strata}}).}
\item{delta_x_tol}{Tolerance below which an adjacent exposure shift
\eqn{\Delta X} is treated as too small to standardize by.}
\item{iptw_truncation}{Lower/upper quantiles for truncating IPTW weights.}
\item{ess_fraction}{Minimum effective-sample-size fraction required of a
stratum before it is used.}
}
\value{
An object of class \code{"mrwin_controls"}: a validated list of settings.
}
\seealso{\code{\link{mrwin}}, \code{\link{mrwin_doubly_ranked_strata}}, \code{\link{mrwin_sdpd}}}
\examples{
# Biobank-scale defaults: compiled kernel + analytic variance
ctrl <- mrwin_controls(n_strata = 10, backend = "fast", inference = "analytic")
str(ctrl, max.level = 1)
}

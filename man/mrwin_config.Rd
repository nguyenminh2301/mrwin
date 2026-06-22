\name{mrwin_config}
\alias{mrwin_config}
\title{Configure the simulation data-generating process}
\description{
Build the parameter object that drives \code{\link{mrwin_simulate}}. It defines a
cardiorenal-style data-generating process with three priority-ranked endpoints
(death, heart failure, renal), a polygenic instrument, an unmeasured confounder,
shared frailty, and right-censoring. Defaults give a protective exposure with
confounding, suitable for studying the estimator's calibration.
}
\usage{
mrwin_config(n_outcome = 6000L, m_snps = 40L, alpha_s = 0.4, alpha_u = 0.8,
             alpha_x = c(-0.4, -0.4, -0.4), gamma_direct = c(0, 0, 0),
             nu_u = c(0.5, 0.5, 0.5), baseline_haz = c(0.02, 0.05, 0.10),
             shape_weibull = 1.2, censoring_rate = 0.05, max_follow_up = 10,
             theta_f = 0.8, maf_low = 0.10, maf_high = 0.40,
             sigma_beta = 0.05, seed = 20260506L)
}
\arguments{
\item{n_outcome}{Number of subjects to simulate.}
\item{m_snps}{Number of independent SNPs in the instrument.}
\item{alpha_s}{Strength of the polygenic score's effect on the exposure
(instrument relevance).}
\item{alpha_u}{Effect of the unmeasured confounder on the exposure.}
\item{alpha_x}{Length-3 vector: causal effect of the exposure on each endpoint's
hazard (negative = protective).}
\item{gamma_direct}{Length-3 vector: direct (pleiotropic) effect of the score on
each endpoint, bypassing the exposure. Zero = a valid instrument.}
\item{nu_u}{Length-3 vector: effect of the confounder on each endpoint's hazard.}
\item{baseline_haz}{Length-3 vector of Weibull baseline rates per endpoint.}
\item{shape_weibull}{Weibull shape parameter shared across endpoints.}
\item{censoring_rate}{Rate of the exponential censoring distribution.}
\item{max_follow_up}{Administrative censoring time (maximum follow-up).}
\item{theta_f}{Gamma-frailty variance inducing within-subject correlation across
endpoints (0 = independent).}
\item{maf_low, maf_high}{Range of minor-allele frequencies (must satisfy
\code{maf_low < maf_high}).}
\item{sigma_beta}{Standard deviation of the per-SNP true effect sizes.}
\item{seed}{Integer seed (or \code{NULL}) for reproducibility.}
}
\value{
An object of class \code{"mrwin_config"}: a validated list of all parameters
above, consumed by \code{\link{mrwin_simulate}}.
}
\seealso{\code{\link{mrwin_simulate}}, \code{\link{mrwin_scenarios}}, \code{\link{mrwin_run_simulation_grid}}}
\examples{
cfg <- mrwin_config(n_outcome = 2000, m_snps = 20, seed = 1)
cfg$alpha_x   # protective exposure on each endpoint
}

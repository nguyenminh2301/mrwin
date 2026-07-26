\name{mrwin_ar_onesample_supci_test}
\alias{mrwin_ar_onesample_supci_test}
\title{Degeneracy-robust Anderson-Rubin point-null test (least-favorable / Andrews-Cheng)}
\description{
When \code{\link{mrwin_ar_onesample}} flags \code{degenerate = TRUE}, the basic
chi-squared Anderson-Rubin test can over-reject: the pairwise U-statistic
moment can degenerate (its first-order/Hajek variance vanishes) and acquire a
non-Gaussian limit, \code{n*U_n -> N(0,4c^2) + sum_k lambda_k(Z_k^2-1)}, so a
single chi-squared critical value is no longer valid. This function implements
the least-favorable (Andrews-Cheng) remedy: because the limiting quantile is
monotone increasing in \code{c^2}, the critical value is taken at the upper
endpoint of a bootstrap confidence interval for \code{c^2}, combined with the
estimated spectrum of the degenerate part via Monte Carlo simulation of the
limit law. This is \strong{uniformly valid but conservative near the degenerate
boundary} (a genuine validity-power frontier, not a tuning artifact): \code{c^2}
is weakly identified there, so the valid upper limit is necessarily wide.

Unlike \code{\link{mrwin_ar_onesample}}, this function tests a single specified
null value \code{beta0} (as validated) rather than returning a confidence
interval: the O(N^2) pairwise kernel matrix and its eigendecomposition needed
for the spectrum are too expensive to recompute at every point of a confidence-set
grid scan. It is intended for the flagged, extreme-tail case, not as the
default inference (the ordinary \code{\link{mrwin_ar_onesample}} CI is exact
and cheap in the non-degenerate regime that covers ordinary practice).
}
\usage{
mrwin_ar_onesample_supci_test(time, status, X, Z, beta0, alpha = 0.05,
                               alpha1 = 0.025, boot_reps = 300L,
                               mc_reps = 4000L, r_eigen = 6L, seed = NULL)
}
\arguments{
\item{time}{N x K matrix of event times, in priority order.}
\item{status}{N x K matrix of 0/1 event status, matching \code{time}.}
\item{X}{Numeric exposure (length N).}
\item{Z}{Numeric scalar instrument (length N).}
\item{beta0}{The null value of the causal gradient to test.}
\item{alpha}{Overall test level (default 0.05); split \code{alpha = alpha1 + alpha2}
(Bonferroni) between the bootstrap CI for \code{c^2} and the simulated critical value.}
\item{alpha1}{One-sided level for the bootstrap CI on \code{c^2} (default 0.025).}
\item{boot_reps}{Bootstrap resamples for the \code{c^2} CI (default 300).}
\item{mc_reps}{Monte Carlo draws for the limiting-law quantile (default 4000).}
\item{r_eigen}{Number of leading eigenvalues (by magnitude) of the degenerate
part retained (default 6; the tail is dropped, a documented approximation).}
\item{seed}{Optional RNG seed (reproducibility).}
}
\value{
A list: \code{stat} (\code{(n*U_n(beta0))^2}), \code{crit} (the
least-favorable critical value), \code{reject} (logical), \code{c2_hi} (the
bootstrap upper bound driving the critical value), \code{n}. Cost is O(N^2) in
memory (the pairwise kernel matrix at \code{beta0}) and O(N^3) for the
eigendecomposition; a subquadratic spectrum estimator is an open research
question, so this function is recommended for moderate N (a few thousand), not
biobank-scale data.
}
\seealso{\code{\link{mrwin_ar_onesample}}}
\examples{
\donttest{
cfg <- mrwin_config(n_outcome = 800, m_snps = 1, seed = 1)
sim <- mrwin_simulate(cfg)
mrwin_ar_onesample_supci_test(sim$time, sim$status, sim$X,
                               as.numeric(sim$G[, 1]), beta0 = 0)
}
}

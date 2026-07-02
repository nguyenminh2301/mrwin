\name{mrwin_ar_onesample}
\alias{mrwin_ar_onesample}
\title{One-sample identification-robust (Anderson-Rubin) inference for win-ratio MR}
\description{
Weak-instrument-robust inference for the causal net-benefit gradient from
individual-level one-sample data, inverting the pairwise Anderson-Rubin moment
\code{AR(b) = n*(a - b*d)^2 / (4*zeta1_hat(b))}, a degree-2 U-statistic analogue
of the linear-IV AR test, valid under any instrument strength (it never divides
by \code{Cov(Z,X)}). Because both the numerator and the plug-in variance
\code{zeta1_hat(b)} are quadratic in \code{b}, the confidence set has a closed
form (the same quadratic-inversion geometry as the Fieller interval), returned
exactly rather than by grid search.

A degeneracy diagnostic is also reported: the pairwise moment can become a
degenerate U-statistic when its first (Hajek) projection variance
\code{zeta1(beta_hat)} is small, in which case the naive chi-squared critical
value over-rejects (see \code{\link{mrwin_ar_onesample_supci_test}}). This
function estimates a percentile-bootstrap confidence interval for the
concentration parameter \code{c^2 = n*zeta1(beta_hat)} (a point estimate alone
is unreliable near the boundary) and flags \code{degenerate = TRUE} when that
interval's lower bound is not clearly bounded away from zero.
}
\usage{
mrwin_ar_onesample(time, status, X, Z, alpha = 0.05,
                    degeneracy_check = TRUE, boot_reps = 200L,
                    alpha1 = 0.025, seed = NULL)
}
\arguments{
\item{time}{N x K matrix of event times, in priority order (column 1 = highest priority).}
\item{status}{N x K matrix of 0/1 event status, matching \code{time}.}
\item{X}{Numeric exposure (length N).}
\item{Z}{Numeric scalar instrument, e.g. a polygenic score (length N).}
\item{alpha}{Test level (default 0.05, i.e. a 95\% confidence set).}
\item{degeneracy_check}{Logical; if \code{TRUE} (default) also compute the
bootstrap degeneracy diagnostic described above.}
\item{boot_reps}{Number of bootstrap resamples for the degeneracy diagnostic
(default 200); each resample re-runs the per-subject kernel, so this is the
dominant cost when \code{degeneracy_check = TRUE}.}
\item{alpha1}{One-sided level for the bootstrap CI on \code{c^2} (default
0.025; the CI's lower endpoint, not the point estimate, drives the flag).}
\item{seed}{Optional RNG seed for the bootstrap (reproducibility).}
}
\value{
A list: \code{gamma} (point estimate), \code{ci} (length-2 bounds, may be
\code{c(-Inf, Inf)}), \code{ci_level}, \code{bounded}, \code{empty} (TRUE if no
\code{b} satisfies the test -- a numerical edge case, not expected in
practice), \code{degenerate} (logical, \code{NA} if
\code{degeneracy_check = FALSE}), \code{c2_hat}, \code{c2_boot_ci} (length-2
bootstrap CI for \code{c^2 = n*zeta1(gamma)}), \code{n}.
}
\seealso{\code{\link{mrwin_ar_onesample_overid}}, \code{\link{mrwin_ar_onesample_supci_test}},
\code{\link{mrwin_winmr_ar}} (the two-sample / summary-data analogue)}
\examples{
\donttest{
cfg <- mrwin_config(n_outcome = 2000, m_snps = 1, seed = 1)
sim <- mrwin_simulate(cfg)
mrwin_ar_onesample(sim$time, sim$status, sim$X, as.numeric(sim$G[, 1]))
}
}

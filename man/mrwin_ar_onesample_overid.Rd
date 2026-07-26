\name{mrwin_ar_onesample_overid}
\alias{mrwin_ar_onesample_overid}
\title{One-sample over-identification (pleiotropy) test for multi-variant win-ratio MR}
\description{
The vector-instrument (e.g. multi-SNP) analogue of
\code{\link{mrwin_ar_onesample}}'s point estimate, giving a
weak-instrument-robust over-identification test: the minimized Anderson-Rubin
statistic \code{min_b AR(b)}, computed with an efficient (two-step GMM) weight
matrix, which is \eqn{\chi^2_{L-1}} under a valid instrument (a
degree-2-U-statistic analogue of the Sargan/Hansen-J test). Rejection signals a
violated exclusion restriction (pleiotropy); the test is, by construction,
blind to pleiotropy proportional to instrument strength (InSIDE-violating,
absorbed into the point estimate) and powered against pleiotropy independent of
instrument strength (InSIDE-satisfying).
}
\usage{
mrwin_ar_onesample_overid(time, status, X, G)
}
\arguments{
\item{time}{N x K matrix of event times, in priority order.}
\item{status}{N x K matrix of 0/1 event status, matching \code{time}.}
\item{X}{Numeric exposure (length N).}
\item{G}{N x L instrument matrix (e.g. L SNP dosages).}
}
\value{
A list: \code{Q} (the minimized AR statistic), \code{Q_df} (\code{L-1}),
\code{Q_p}, \code{gamma0} (the preliminary point estimate used to form the
weight matrix), \code{n_snp}.
}
\section{Many-instrument warning}{
This construction estimates an L x L weight matrix from N observations and is
validated only for L small relative to N. A stress simulation
(\code{tools/validation-scripts/p3-manyweak-instrument-probe.R}) shows type-I
error under a valid instrument is calibrated at L/N = 0.016 (0.02) and
L/N = 0.04 (0.09, mildly inflated), but breaks down sharply beyond that: 0.71
at L/N = 0.2 and 1.00 (always rejects) at L/N = 0.5. This mirrors the
many-weak-instrument problem known for the linear-IV analogue of this test
(Kleibergen, 2002, Econometrica), for which a correction exists in the linear
case; no such correction is implemented here. This function warns when
L/N > 0.02 and errors when L/N > 0.15, thresholds set from the simulation
above, not a formal bound.
}
\seealso{\code{\link{mrwin_ar_onesample}}, \code{\link{mrwin_winmr_ar}} (the
two-sample analogue, \code{Q = min_b AR(b)} there is the same construction on
summary statistics)}
\examples{
\donttest{
cfg <- mrwin_config(n_outcome = 2000, m_snps = 20, seed = 1)
sim <- mrwin_simulate(cfg)
mrwin_ar_onesample_overid(sim$time, sim$status, sim$X, sim$G)
}
}

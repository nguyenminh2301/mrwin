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
function computes a lower confidence bound for the concentration parameter
\code{c^2 = n*zeta1(beta_hat)} (a point estimate alone is unreliable near the
boundary) and flags \code{degenerate = TRUE} when that bound is not clearly
bounded away from zero. The bound is obtained either by percentile bootstrap
(default) or in closed form (\code{degeneracy_method = "analytic"}, about 100x
faster); see the Degeneracy diagnostic section below.
}
\usage{
mrwin_ar_onesample(time, status, X, Z, alpha = 0.05,
                    degeneracy_check = TRUE, boot_reps = 200L,
                    alpha1 = 0.025, seed = NULL,
                    degeneracy_method = c("bootstrap", "analytic"))
}
\arguments{
\item{time}{N x K matrix of event times, in priority order (column 1 = highest priority).}
\item{status}{N x K matrix of 0/1 event status, matching \code{time}.}
\item{X}{Numeric exposure (length N).}
\item{Z}{Numeric scalar instrument, e.g. a polygenic score (length N).}
\item{alpha}{Test level (default 0.05, i.e. a 95\% confidence set).}
\item{degeneracy_check}{Logical; if \code{TRUE} (default) also compute the
degeneracy diagnostic described above.}
\item{boot_reps}{Number of bootstrap resamples for the degeneracy diagnostic
(default 200); each resample re-runs the per-subject kernel, so this is the
dominant cost when \code{degeneracy_check = TRUE} and
\code{degeneracy_method = "bootstrap"}.}
\item{alpha1}{One-sided level for the bound on \code{c^2} (default 0.025; the
lower endpoint, not the point estimate, drives the flag).}
\item{seed}{Optional RNG seed for the bootstrap (reproducibility). Ignored
when \code{degeneracy_method = "analytic"}, which is deterministic.}
\item{degeneracy_method}{How to compute the lower bound on \code{c^2}:
\code{"bootstrap"} (default) resamples \code{boot_reps} times;
\code{"analytic"} uses the closed-form standard error described below, at
about 1/100th of the cost.}
}
\value{
A list: \code{gamma} (point estimate), \code{ci} (length-2 bounds, may be
\code{c(-Inf, Inf)}), \code{ci_level}, \code{bounded}, \code{empty} (TRUE if no
\code{b} satisfies the test -- a numerical edge case, not expected in
practice), \code{degenerate} (logical, \code{NA} if
\code{degeneracy_check = FALSE}), \code{c2_hat}, \code{c2_boot_ci} (length-2
bound for \code{c^2 = n*zeta1(gamma)}; with
\code{degeneracy_method = "analytic"} only the lower endpoint is computed and
the upper is \code{NA}), \code{c2_se} (the analytic standard error, \code{NA}
under the bootstrap), \code{degeneracy_method}, \code{n}.
}
\section{Degeneracy diagnostic -- the two methods}{
Because the resampling loop holds \code{gamma} \strong{fixed}, the bootstrapped
quantity is exactly \code{Var_i(g_i)} with \code{g_i = gh_i - gamma*gx_i}, so
no delta-method term for \code{gamma_hat}'s own randomness is needed and a
closed-form standard error follows from the influence function of a variance,
\code{IF(g) = (g - mu)^2 - sigma^2}:
\deqn{SE(c^2) = \sqrt{n (m_4 - m_2^2)}}
with \code{m_2}, \code{m_4} the second and fourth central moments of
\code{g_i}. The \code{"analytic"} bound is
\code{c2_hat - qnorm(1-alpha1)*SE}.

Validated against the 200-replication bootstrap over 840 simulated cells
(\code{tools/validation-scripts/p3-analytic-c2-bound.R}; N in {800, 2500} x 7
censoring rates x weak/strong instruments): the two bounds correlate at 0.9999
(Pearson and Spearman), the analytic bound is very slightly conservative
(median ratio 0.983, erring toward flagging degenerate -- the safe direction),
the resulting \code{degenerate} flags \strong{agree on 99.64\% of cells} (1
false-safe call in 840, 0.3\% of its "safe" calls), and it is about
\strong{101x faster}.

Caveat: \code{g_i} are U-statistic projections rather than iid draws, which the
bootstrap captures exactly and the influence-function approximation does not.
The agreement above is empirical evidence that this gap is immaterial in the
tested regimes, not a proof that it always is; the bootstrap remains the
default for that reason.
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

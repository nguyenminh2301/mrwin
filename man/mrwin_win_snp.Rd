\name{mrwin_win_snp}
\alias{mrwin_win_snp}
\title{Per-SNP win-odds coefficient with a closed-form standard error}
\description{
Estimates the log-win-odds slope of a priority-ranked (hierarchical) outcome on
allele dosage for a single SNP -- the outcome-side summary statistic for
two-sample win-ratio Mendelian randomization. The point estimate uses the
subquadratic hierarchical kernel; the standard error is the influence-function
variance of the log-win-odds (the same one used by
\code{\link{mrwin_analytic_covariance}}), evaluated on a capped subject subsample
so the whole routine stays near-linear in sample size.
}
\usage{
mrwin_win_snp(time, status, dosage, weights = NULL,
              split = c("median", "threshold"), threshold = NULL,
              max_subjects = 2000L, seed = NULL, floor = 1e-9)
}
\arguments{
\item{time, status}{N x K matrices of event times and 0/1 status in priority
order (column 1 = highest priority).}
\item{dosage}{Numeric allele dosage (length N); any scale is allowed.}
\item{weights}{Optional non-negative subject weights (length N).}
\item{split}{How to dichotomize dosage into the high/low groups whose win-odds is
contrasted: \code{"median"} (default) or \code{"threshold"}.}
\item{threshold}{Cut value for \code{split = "threshold"} (high = dosage >
threshold).}
\item{max_subjects}{Subject-subsample cap per group used for the variance
(default 2000); the point estimate always uses all subjects.}
\item{seed}{Optional RNG seed for the variance subsample.}
\item{floor}{Minimum weighted win/loss count; below it the estimate is \code{NA}.}
}
\value{
A list with \code{delta} (win-odds slope per dosage unit) and its \code{se}, plus
\code{log_theta}, \code{se_log_theta}, \code{d_diff} (mean-dosage gap between
groups), \code{wins}, \code{losses}, \code{n_high}, \code{n_low}.
}
\seealso{\code{\link{mrwin_win_gwas}}, \code{\link{mrwin_twosample_ivw}}}
\examples{
\donttest{
cfg <- mrwin_config(n_outcome = 3000, m_snps = 10, seed = 1)
sim <- mrwin_simulate(cfg)
mrwin_win_snp(sim$time, sim$status, sim$G[, 1])
}
}

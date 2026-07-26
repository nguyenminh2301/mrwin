\name{mrwin_win_gwas}
\alias{mrwin_win_gwas}
\title{Win-odds GWAS: per-SNP win-odds coefficients across a SNP panel}
\description{
Applies \code{\link{mrwin_win_snp}} to each column of a genotype/dosage matrix,
producing the shareable outcome-side summary statistics (per-SNP win-odds
coefficients and their standard errors) for two-sample win-ratio Mendelian
randomization.
}
\usage{
mrwin_win_gwas(time, status, G, weights = NULL, max_subjects = 2000L,
               seed = NULL, ...)
}
\arguments{
\item{time, status}{N x K matrices of event times and 0/1 status in priority
order (column 1 = highest priority).}
\item{G}{N x M genotype/dosage matrix (one column per SNP).}
\item{weights}{Optional non-negative subject weights (length N).}
\item{max_subjects}{Subject-subsample cap per group used for the variance.}
\item{seed}{Optional RNG seed for the variance subsample.}
\item{...}{Passed to \code{\link{mrwin_win_snp}} (e.g. \code{split}).}
}
\value{
A data frame with one row per SNP: \code{snp}, \code{delta}, \code{se},
\code{log_theta}, \code{se_log_theta}, \code{d_diff}, \code{n_high}, \code{n_low}.
}
\seealso{\code{\link{mrwin_win_snp}}, \code{\link{mrwin_twosample_ivw}}}
\examples{
\donttest{
cfg <- mrwin_config(n_outcome = 3000, m_snps = 10, seed = 1)
sim <- mrwin_simulate(cfg)
mrwin_win_gwas(sim$time, sim$status, sim$G)
}
}

\name{mrwin_winmr_ar}
\alias{mrwin_winmr_ar}
\title{Identification-robust (Anderson-Rubin) inference for two-sample win-ratio MR}
\description{
Weak-instrument-robust inference for the causal win-odds gradient from two-sample
summary statistics, inverting the pairwise Anderson-Rubin moment
\code{AR(b) = sum_l (delta_l - b*beta_gx_l)^2 / (se_gy_l^2 + b^2 se_gx_l^2)}, which is
\eqn{\chi^2_L} at the true gradient regardless of instrument strength (\code{L} =
number of SNPs). The confidence set inverts the test
(\code{{b : AR(b) <= qchisq(1-alpha, L)}}) and stays valid under a weak instrument,
returning an honest unbounded set rather than a falsely tight one (it never divides
by \code{Cov(Z,X)}). With one SNP it reduces to the Fieller interval; with several it
also yields an over-identification / pleiotropy test (\code{min_b AR(b)} \eqn{\sim
\chi^2_{L-1}}). The \code{delta_gy} are per-SNP win-odds coefficients from
\code{\link{mrwin_win_gwas}}.
}
\usage{
mrwin_winmr_ar(beta_gx, delta_gy, se_gy, se_gx = NULL,
               alpha = 0.05, grid = NULL, n_grid = 4000L)
}
\arguments{
\item{beta_gx}{Per-SNP exposure associations.}
\item{delta_gy}{Per-SNP win-odds coefficients (e.g. \code{mrwin_win_gwas()$delta}).}
\item{se_gy}{Standard errors of \code{delta_gy}.}
\item{se_gx}{Optional standard errors of \code{beta_gx} (default: exposure treated as
fixed, the usual two-sample convention).}
\item{alpha}{Test level (default 0.05, i.e. a 95\% confidence set).}
\item{grid}{Optional numeric vector of gradient values to scan; by default a grid
spanning about 12 IVW standard errors around the IVW estimate.}
\item{n_grid}{Grid size when \code{grid} is not supplied (default 4000).}
}
\value{
A list: \code{gamma} (AR point estimate = argmin AR), \code{ci} (length-2 bounds of
the accepted set), \code{bounded} (\code{FALSE} = weak instrument: the set touches a
grid edge or is disjoint), \code{contiguous}, \code{n_snp}, and the
over-identification test \code{Q} = min AR, \code{Q_df} = L-1, \code{Q_p}
(pleiotropy / invalid-instrument diagnostic).
}
\seealso{\code{\link{mrwin_win_gwas}}, \code{\link{mrwin_twosample_ivw}}}
\examples{
\donttest{
cfg <- mrwin_config(n_outcome = 3000, m_snps = 10, seed = 1)
sim <- mrwin_simulate(cfg)
wg <- mrwin_win_gwas(sim$time, sim$status, sim$G)
beta_gx <- apply(sim$G, 2, function(d) stats::cov(sim$X, d) / stats::var(d))
mrwin_winmr_ar(beta_gx, wg$delta, wg$se)
}
}

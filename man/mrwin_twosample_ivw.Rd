\name{mrwin_twosample_ivw}
\alias{mrwin_twosample_ivw}
\title{Two-sample win-ratio MR by inverse-variance weighting}
\description{
Pools per-SNP exposure coefficients \code{beta_gx} (from an exposure GWAS) with
per-SNP win-odds coefficients \code{delta_gy} (from \code{\link{mrwin_win_gwas}})
into a single causal win-odds gradient, using fixed-effect inverse-variance
weighting -- the win-ratio analogue of ordinary two-sample MR. By default the
exposure coefficients are treated as fixed (the usual two-sample convention);
supplying \code{se_gx} adds the second-order Wald-ratio variance term.
}
\usage{
mrwin_twosample_ivw(beta_gx, delta_gy, se_gy, se_gx = NULL)
}
\arguments{
\item{beta_gx}{Per-SNP exposure associations.}
\item{delta_gy}{Per-SNP win-odds coefficients (e.g. \code{mrwin_win_gwas()$delta}).}
\item{se_gy}{Standard errors of \code{delta_gy} (e.g. \code{mrwin_win_gwas()$se}).}
\item{se_gx}{Optional standard errors of \code{beta_gx}.}
}
\value{
A list: \code{gamma} (pooled win-odds gradient per exposure unit), \code{se},
\code{ci95}, \code{n_snp}, and Cochran's \code{Q}, \code{Q_df}, \code{Q_p}
heterogeneity test.
}
\seealso{\code{\link{mrwin_win_gwas}}, \code{\link{mrwin_win_snp}}}
\examples{
\donttest{
cfg <- mrwin_config(n_outcome = 3000, m_snps = 10, seed = 1)
sim <- mrwin_simulate(cfg)
wg <- mrwin_win_gwas(sim$time, sim$status, sim$G)
beta_gx <- apply(sim$G, 2, function(d) stats::cov(sim$X, d) / stats::var(d))
mrwin_twosample_ivw(beta_gx, wg$delta, wg$se)
}
}

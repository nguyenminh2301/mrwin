\name{mrwin_simulate}
\alias{mrwin_simulate}
\title{Simulate a hierarchical-endpoint MR cohort}
\description{
Draw one cohort from the data-generating process defined by
\code{\link{mrwin_config}}: genotypes and a polygenic score, an exposure
confounded by an unmeasured variable, and three priority-ranked, possibly
censored, frailty-correlated time-to-event endpoints. Use it to try the
estimator without real data and to study its calibration.
}
\usage{
mrwin_simulate(config = mrwin_config(), seed = config$seed)
}
\arguments{
\item{config}{A configuration object from \code{\link{mrwin_config}}.}
\item{seed}{Optional integer seed (defaults to \code{config$seed}); set
\code{NULL} to leave the RNG state untouched.}
}
\value{
A list with the simulated cohort, including \code{time} and \code{status}
(\eqn{n \times 3} endpoint matrices), \code{G} (genotype matrix), \code{X} (the
exposure), \code{true_betas} (per-SNP instrument weights), \code{S_true} (the
true polygenic score), \code{U} (the confounder), \code{W} (frailty), and the
\code{config} used.
}
\seealso{\code{\link{mrwin_config}}, \code{\link{mrwin}}, \code{\link{mrwin_run_simulation_grid}}}
\examples{
sim <- mrwin_simulate(mrwin_config(n_outcome = 500, m_snps = 20, seed = 1))
fit <- mrwin(
  endpoint = mrwin_endpoint(sim$time, sim$status, c("death", "hf", "renal")),
  genotype = sim$G, exposure = sim$X,
  gwas = mrwin_gwas(sim$true_betas, rep(0.01, 20)),
  controls = mrwin_controls(n_strata = 5, bootstrap = 100, backend = "fast")
)
fit
}

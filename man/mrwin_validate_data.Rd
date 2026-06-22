\name{mrwin_validate_data}
\alias{mrwin_validate_data}
\title{Validate and assemble \code{mrwin} inputs}
\description{
Check the endpoint, genotype, exposure and GWAS inputs for shape, missingness,
event coding and the terminal-event ordering, and return a clean, aligned data
object. \code{\link{mrwin}} calls this internally; run it yourself first to catch
data problems before estimation.
}
\usage{
mrwin_validate_data(data = NULL, endpoint, genotype, exposure, gwas = NULL,
                    beta_gwas = NULL, se_gwas = NULL, covariates = NULL,
                    controls = mrwin_controls(),
                    complete_cases = c("error", "drop"),
                    allow_time_zero = TRUE, terminal_priority = 1L,
                    terminal_tolerance = 0, check_terminal_order = TRUE,
                    check_strata = TRUE)
}
\arguments{
\item{data}{Optional data frame; when given, \code{endpoint} columns may be
referenced by name.}
\item{endpoint}{An \code{\link{mrwin_endpoint}} specification.}
\item{genotype}{An \eqn{n \times M} genotype (or polygenic-component) matrix.}
\item{exposure}{Numeric exposure vector of length \eqn{n}.}
\item{gwas}{An \code{\link{mrwin_gwas}} specification (or supply \code{beta_gwas}
and \code{se_gwas} directly).}
\item{beta_gwas, se_gwas}{Alternative to \code{gwas}: per-SNP weights and SEs.}
\item{covariates}{Optional covariate matrix for adjustment.}
\item{controls}{An \code{\link{mrwin_controls}} object (used for stratum checks).}
\item{complete_cases}{\code{"error"} (default) to stop on missing rows, or
\code{"drop"} to remove them.}
\item{allow_time_zero}{Logical; permit event/censoring times equal to zero.}
\item{terminal_priority}{Index of the terminal (death) priority level.}
\item{terminal_tolerance}{Numeric tolerance when checking the terminal-event
ordering.}
\item{check_terminal_order}{Logical; verify that the terminal event dominates
lower-priority times as required by the hierarchy.}
\item{check_strata}{Logical; verify that \code{n_strata} is feasible given \eqn{n}.}
}
\value{
An object of class \code{"mrwin_validated_data"}: the aligned, checked inputs
ready for estimation. Has a \code{print} method.
}
\seealso{\code{\link{mrwin}}, \code{\link{mrwin_endpoint}}, \code{\link{mrwin_gwas}}}
\examples{
sim <- mrwin_simulate(mrwin_config(n_outcome = 300, m_snps = 10, seed = 1))
vd <- mrwin_validate_data(
  endpoint = mrwin_endpoint(sim$time, sim$status, c("death", "hf", "renal")),
  genotype = sim$G, exposure = sim$X,
  gwas = mrwin_gwas(sim$true_betas, rep(0.01, 10))
)
vd
}

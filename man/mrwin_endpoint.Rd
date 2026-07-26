\name{mrwin_endpoint}
\alias{mrwin_endpoint}
\title{Specify hierarchical (priority-ranked) endpoints}
\description{
Bundle the time-to-event and event-indicator data for a set of clinically
ranked composite endpoints into the object expected by \code{\link{mrwin}}.
Outcomes are listed from highest to lowest priority (e.g.\ death first); the
hierarchical win/loss kernel decides each pair at the first priority level that
separates the two subjects.
}
\usage{
mrwin_endpoint(time, status, priority = NULL)
}
\arguments{
\item{time}{An \eqn{n \times K} numeric matrix of follow-up / event times (one
column per priority level, ordered most to least severe), or a length-\eqn{K}
character vector of column names when a \code{data} frame is passed to
\code{\link{mrwin}}.}
\item{status}{An \eqn{n \times K} 0/1 event-indicator matrix (1 = event
observed, 0 = censored) with the same shape as \code{time}, or the
corresponding column names.}
\item{priority}{Optional character labels for the \eqn{K} levels, most severe
first. Defaults to the column names of \code{time}, otherwise \code{"p1"},
\code{"p2"}, \dots}
}
\value{
An object of class \code{"mrwin_endpoint_spec"}: a list with elements
\code{time}, \code{status} and \code{priority}, suitable for the \code{endpoint}
argument of \code{\link{mrwin}}.
}
\seealso{\code{\link{mrwin}}, \code{\link{mrwin_gwas}}, \code{\link{mrwin_validate_data}}}
\examples{
sim <- mrwin_simulate(mrwin_config(n_outcome = 200, m_snps = 10, seed = 1))
ep <- mrwin_endpoint(sim$time, sim$status, priority = c("death", "hf", "renal"))
str(ep)
}

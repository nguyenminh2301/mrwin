\name{mrwin_report}
\alias{mrwin_report}
\title{Generate mrwin Analysis Report}
\description{
Generates a text or markdown report from an \code{mrwin_fit} object,
summarizing data, estimates, diagnostics, and warnings.
}
\usage{
mrwin_report(x, file = NULL, format = c("text", "markdown"))
}
\arguments{
\item{x}{An \code{mrwin_fit} object from \code{\link{mrwin}}.}
\item{file}{Optional file path to write the report. If NULL, prints to console.}
\item{format}{Output format: \code{"text"} or \code{"markdown"}.}
}
\value{
Invisibly returns the report as a character string.
}
\seealso{
\code{\link{mrwin}}, \code{\link{tidy.mrwin_fit}}
}

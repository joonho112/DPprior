% Please edit documentation in R/14_diagnostics.R
\name{compare_diagnostics}
\alias{compare_diagnostics}
\title{Compare Diagnostics Across Multiple Fits}
\usage{
compare_diagnostics(..., M = .QUAD_NODES_DEFAULT)
}
\arguments{
\item{...}{Named DPprior fit objects.}

\item{M}{Reserved historical argument. It must remain at the canonical
default because each fit's retained selected order is authoritative.}
}
\value{
Data frame of estimand-labelled diagnostics.
}
\description{
Compare Diagnostics Across Multiple Fits
}
\keyword{internal}

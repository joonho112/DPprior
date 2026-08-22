% Please edit documentation in R/14_diagnostics.R
\name{.get_K_pmf_support}
\alias{.get_K_pmf_support}
\title{Get K_J PMF and status metadata}
\usage{
.get_K_pmf_support(
  J,
  a,
  b,
  M = .QUAD_NODES_DEFAULT,
  M_verify = NULL,
  abs_tol = 1e-10,
  rel_tol = 1e-08
)
}
\arguments{
\item{J}{Integer sample size.}

\item{a, b}{Numeric Gamma hyperparameters.}

\item{M}{Selected quadrature order.}

\item{M_verify}{Optional verification order.}

\item{abs_tol, rel_tol}{Comparison tolerances.}
}
\value{
PMF, support, optional verification PMF, and provenance.
}
\description{
Get K_J PMF and status metadata
}
\keyword{internal}

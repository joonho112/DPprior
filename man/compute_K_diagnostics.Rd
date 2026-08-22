% Please edit documentation in R/14_diagnostics.R
\name{compute_K_diagnostics}
\alias{compute_K_diagnostics}
\title{K Distribution Diagnostics}
\usage{
compute_K_diagnostics(
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
A canonical K_J diagnostic component.
}
\description{
K Distribution Diagnostics
}
\keyword{internal}

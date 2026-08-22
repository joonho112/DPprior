% Please edit documentation in R/14_diagnostics.R
\name{compute_coclustering_diagnostics}
\alias{compute_coclustering_diagnostics}
\title{Co-Clustering Diagnostics (rho)}
\usage{
compute_coclustering_diagnostics(
  a,
  b,
  M = .QUAD_NODES_DEFAULT,
  M_verify = NULL,
  abs_tol = 1e-10,
  rel_tol = 1e-08
)
}
\arguments{
\item{a, b}{Numeric Gamma hyperparameters.}

\item{M}{Selected quadrature order.}

\item{M_verify}{Optional verification order.}

\item{abs_tol, rel_tol}{Comparison tolerances.}
}
\value{
Status-aware rho component without a qualitative category.
}
\description{
Co-Clustering Diagnostics (rho)
}
\keyword{internal}

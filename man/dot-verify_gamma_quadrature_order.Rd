% Please edit documentation in R/02_quadrature.R
\name{.verify_gamma_quadrature_order}
\alias{.verify_gamma_quadrature_order}
\title{Independently Verify a Gamma Quadrature Integral at Higher Order}
\usage{
.verify_gamma_quadrature_order(
  f,
  a,
  b,
  M = .QUAD_NODES_DEFAULT,
  M_verify = NULL,
  abs_tol = 1e-12,
  rel_tol = 1e-10
)
}
\arguments{
\item{f}{Function to integrate.}

\item{a, b}{Gamma shape and rate.}

\item{M}{Selected quadrature order.}

\item{M_verify}{Independent verification order. The minimum admissible order
is \code{max(2*M, M+40)}. When \code{NULL}, that minimum is used. Because
the supported ceiling is 512 nodes, verification is unavailable when
\code{M > 256}; the function then raises a typed verification error even
if an explicit \code{M_verify} is supplied.}

\item{abs_tol, rel_tol}{Non-negative finite tolerances.}
}
\value{
A list containing selected and verification estimates, absolute
difference, tolerance, pass flag, and versioned metadata.
}
\description{
Internal helper that evaluates the same contracted integrand at selected
and higher node counts, records both computation metadata objects, and
applies a combined absolute/relative tolerance to their difference.
}
\keyword{internal}

% Please edit documentation in R/07_jacobian.R
\name{verify_jacobian}
\alias{verify_jacobian}
\title{Verify Jacobian Against Finite Differences}
\usage{
verify_jacobian(
  J,
  a,
  b,
  eps = 1e-06,
  M = .QUAD_NODES_VERIFICATION,
  verbose = TRUE,
  abs_tol = .JACOBIAN_DERIV_ABS_TOL,
  rel_tol = .JACOBIAN_DERIV_REL_TOL
)
}
\arguments{
\item{J}{Integer; sample size.}

\item{a}{Numeric; shape parameter.}

\item{b}{Numeric; rate parameter.}

\item{eps}{Positive relative step for five-point finite differences
(default: 1e-6).}

\item{M}{Integer; number of quadrature nodes (default: 200 for verification).}

\item{verbose}{Logical; if TRUE, print detailed comparison.}

\item{abs_tol, rel_tol}{Non-negative mixed-error tolerances; at least one
must be positive.}
}
\value{
A named list with components:
\describe{
\item{\code{analytic}}{The score-identity Jacobian (legacy component
name retained for compatibility)}
\item{\code{numeric}}{The numerically computed Jacobian (finite differences)}
\item{\code{abs_error}}{Matrix of absolute errors}
\item{\code{rel_error}}{Matrix of relative errors}
\item{\code{max_rel_error}}{Maximum relative error across all entries}
\item{\code{status}}{One of \code{converged}, \code{approximate}, or
\code{failed}.}
\item{\code{component_status}}{Named statuses for order refinement and
finite-difference agreement.}
\item{\code{pass}}{Compatibility logical equal to
\code{status == "converged"}.}
\item{\code{conditioning}}{The structured conditioning result from the
score-identity Jacobian.}
}
}
\description{
Compares the score-identity Jacobian against numerical finite differences to
validate the implementation.
}
\details{
Uses the fourth-order, five-point central formula at a fixed recorded
quadrature order. The algebra is independent of the score identity, although
both computations use the package's marginal-moment integration layer.

\strong{Important:} Release tests additionally compare both methods against
adaptive log-Gamma integration with separately coded conditional moments.
}
\examples{
# Verify Jacobian for a specific case
result <- verify_jacobian(J = 50, a = 2.0, b = 1.0, verbose = TRUE)

}

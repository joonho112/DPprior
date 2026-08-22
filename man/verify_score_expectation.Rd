% Please edit documentation in R/07_jacobian.R
\name{verify_score_expectation}
\alias{verify_score_expectation}
\title{Verify Score Function Zero Expectation Property}
\usage{
verify_score_expectation(
  a,
  b,
  M = .QUAD_NODES_VERIFICATION,
  verbose = TRUE,
  tolerance = 0.001
)
}
\arguments{
\item{a}{Numeric; shape parameter.}

\item{b}{Numeric; rate parameter.}

\item{M}{Integer; number of quadrature nodes.}

\item{verbose}{Logical; if TRUE, print results.}

\item{tolerance}{Non-negative absolute tolerance for the raw quadrature
expectation of each score.}
}
\value{
A named list with components:
\describe{
\item{\code{E_score_a}}{Expectation of \eqn{s_a}}
\item{\code{E_score_b}}{Expectation of \eqn{s_b}}
\item{\code{status}}{\code{converged} when both raw score expectations
meet \code{tolerance}; otherwise \code{approximate}.}
\item{\code{abs_error}}{Named absolute raw score-expectation errors.}
\item{\code{passed}}{Compatibility logical equal to
\code{status == "converged"}.}
}
}
\description{
Verifies the fundamental property that \eqn{E[s_\theta(\alpha)] = 0}
for both score functions.
}
\details{
This is a fundamental property of score functions. Due to quadrature
approximation error, the computed expectations may not be exactly zero.

\strong{Expected behavior:}
\itemize{
\item \code{E[s_b]} should be very close to zero (typically < 1e-14)
because s_b is linear in alpha.
\item \code{E[s_a]} can show material error for small shape because the
\code{log(alpha)} singularity converges slowly under a finite rule.
}

The production Jacobian uses an exact zero-limit control variate, so raw
score-expectation error is reported rather than silently added to every
derivative. Release tests use adaptive log-Gamma integration as an
independent score-identity reference.
}
\examples{
\dontrun{
verify_score_expectation(a = 2.0, b = 1.0, verbose = TRUE)

}
}
\keyword{internal}

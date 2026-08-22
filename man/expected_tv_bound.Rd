% Please edit documentation in R/13_error_bounds.R
\name{expected_tv_bound}
\alias{expected_tv_bound}
\title{Expected TV Bound Under Gamma Prior}
\usage{
expected_tv_bound(
  J,
  a,
  b,
  cJ = log(J),
  M = .QUAD_NODES_DEFAULT,
  M_verify = NULL,
  abs_tol = 1e-10,
  rel_tol = 1e-08,
  strict = FALSE
)
}
\arguments{
\item{J}{Integer; sample size.}

\item{a, b}{Numeric; Gamma hyperparameters.}

\item{cJ}{Numeric; scaling constant (default: log(J)).}

\item{M}{Integer; number of nodes for a secondary Gauss-Laguerre audit.}

\item{M_verify}{Optional independent order at least
\code{max(2*M, M+40)} and no greater than 512 for the secondary audit.}

\item{abs_tol, rel_tol}{Non-negative adaptive selected-versus-verification
tolerances. They are also used and reported for the Gauss-Laguerre audit.}

\item{strict}{Logical; if \code{TRUE}, require successful agreement between
two independently controlled adaptive integrations and otherwise raise a
typed
\code{dpprior_tv_bound_convergence_error}.}
}
\value{
Numeric; \eqn{E[d_{TV} \text{ bound} | a, b]}. The scalar carries
a \code{"tv_bound_metadata"} attribute with status, numerical error and
omitted-tail allowances, discrepancies, a Gauss-Laguerre audit, and
theorem provenance. The returned conservative adaptive upper bound is
never silently replaced by either verification value.
}
\description{
Integrates the conditional TV bound over \eqn{\alpha \sim \text{Gamma}(a, b)}
to obtain the marginal error bound.
}
\details{
From the submitted manuscript Appendix D, Theorem
\code{thm:marginal-tv} (Equation D17), the TV error between the exact prior
predictive \eqn{p(S_J | a, b)} and the A1 Negative-Binomial proxy is bounded by:
\deqn{d_{TV}(P^{\text{exact}}, Q^{A1}) \le E_{\alpha \sim \Gamma(a,b)}[B_{\text{Pois}} + B_{\text{lin}}]}

This follows from the mixture contraction property of TV distance. The
expectation is evaluated on log-alpha by adaptive quadrature over central
Gamma quantiles. The returned value adds both the integrator's reported
absolute error and the full omitted Gamma-tail probability to the numerical
estimate. A tighter independently controlled integration determines the
convergence status. Fixed-order Gauss-Laguerre values are retained in the
metadata as a reproducibility audit because their convergence can be
non-monotone for this capped, non-smooth integrand.
}
\examples{
\dontrun{
# Marginal TV bound for J=100, Gamma(1, 1)
expected_tv_bound(J = 100, a = 1, b = 1)

}
}
\seealso{
\code{\link{compute_total_tv_bound}}, \code{\link{integrate_gamma}}
}
\keyword{internal}

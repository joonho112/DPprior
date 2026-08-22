% Please edit documentation in R/09_weights_rho.R
\name{rrho}
\alias{rrho}
\title{Random Generation from rho Distribution}
\usage{
rrho(
  n,
  a,
  b,
  n_sticks = 500L,
  remainder_tol = 1e-10,
  max_sticks = 100000L,
  strict = TRUE
)
}
\arguments{
\item{n}{Integer; number of samples to generate.}

\item{a}{Numeric; shape parameter of the Gamma prior on alpha (a > 0).}

\item{b}{Numeric; rate parameter of the Gamma prior on alpha (b > 0).}

\item{n_sticks}{Integer; minimum number of sticks generated for each draw.
The sampler continues adaptively when the remaining squared mass is above
\code{remainder_tol}. Default is 500.}

\item{remainder_tol}{Numeric; maximum allowed deterministic upper bound on
the omitted contribution to rho. Because all ungenerated weights sum to
remainder \eqn{R}, their contribution satisfies
\eqn{0 <= \sum_{\mathrm{tail}} w_h^2 <= R^2}. Default is \code{1e-10}.}

\item{max_sticks}{Integer; hard per-draw stick ceiling. Default is 100,000.}

\item{strict}{Logical; if \code{TRUE} (the default), reaching
\code{max_sticks} before the remainder contract is met raises a typed
\code{dpprior_rho_truncation_error}. If \code{FALSE}, the finite lower
approximation is returned with status \code{"approximate"} and its
deterministic remainder bound in the \code{"rrho_diagnostics"}
attribute.}
}
\value{
Numeric vector of length n. The \code{"rrho_diagnostics"} attribute
records the adaptive method, status, tolerance, sticks used, and the
per-draw deterministic omitted-rho bounds on both ordinary and log scales.
When a positive bound is smaller than the floating-point range, its
ordinary-scale representation is rounded upward to the smallest positive
double rather than silently reported as zero.
}
\description{
Generates random samples from the rho = sum w_h^2 distribution by
stick-breaking simulation.
}
\details{
Uses the hierarchical representation:
\enumerate{
\item alpha ~ Gamma(a, b)
\item v_h | alpha ~ Beta(1, alpha) independently
\item w_1 = v_1, w_h = v_h * prod(1 - v_l) for l < h
\item Generate at least \code{n_sticks} terms and continue until the
unallocated remainder \eqn{R} obeys \eqn{R^2 <= remainder_tol}
\item Return the partial sum; its omitted contribution is in
\eqn{[0,R^2]}
}

Thus finite-stick error is controlled draw by draw rather than being
silently discarded. Useful for Monte Carlo validation of analytical
formulas.
}
\examples{
set.seed(42)
rho_samples <- rrho(1000, a = 2, b = 1)
mean(rho_samples)
mean_rho(a = 2, b = 1)

}
\seealso{
\code{\link{rw1}} for w1 random generation
}

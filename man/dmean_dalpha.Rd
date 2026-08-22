% Please edit documentation in R/03_moments_conditional.R
\name{dmean_dalpha}
\alias{dmean_dalpha}
\title{Derivative of Conditional Mean w.r.t. Alpha}
\usage{
dmean_dalpha(J, alpha)
}
\arguments{
\item{J}{Sample size (integer >= 1).}

\item{alpha}{Concentration parameter (positive numeric, vectorized).}
}
\value{
Numeric vector of derivatives (same length as \code{alpha}).
}
\description{
Computes \eqn{\frac{d}{d\alpha} \mathbb{E}[K_J \mid \alpha]}.
This is used internally when building Jacobians for Newton-type solvers.
}
\details{
A special-function closed form is:
\deqn{\frac{d}{d\alpha}\mu_J(\alpha) =
      \{\psi(\alpha+J)-\psi(\alpha)\} + \alpha\{\psi_1(\alpha+J)-\psi_1(\alpha)\}}

This derivative is positive for \eqn{\alpha > 0} when \code{J > 1},
confirming that \eqn{E[K_J | \alpha]} is then strictly increasing in
\eqn{\alpha}. At \code{J = 1}, the mean is identically one and the
derivative is zero.
The implementation uses the stable finite sum
\deqn{\sum_{r=1}^{J-1} r/(\alpha+r)^2}. In particular, both the mean and
variance derivatives approach \eqn{H_{J-1}} as \eqn{\alpha\to0^+}; they do
not approach zero.
}
\examples{
\dontrun{
dmean_dalpha(50, 2.0)

# Verify with finite difference
J <- 50; alpha <- 2.0; eps <- 1e-6
fd <- (mean_K_given_alpha(J, alpha + eps) -
       mean_K_given_alpha(J, alpha - eps)) / (2 * eps)
abs(fd - dmean_dalpha(J, alpha)) < 1e-5  # TRUE

}
}
\keyword{internal}

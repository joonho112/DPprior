% Please edit documentation in R/03_moments_conditional.R
\name{var_K_given_alpha}
\alias{var_K_given_alpha}
\title{Conditional Variance of K_J Given Alpha}
\usage{
var_K_given_alpha(J, alpha)
}
\arguments{
\item{J}{Sample size (integer >= 1).}

\item{alpha}{Concentration parameter (positive numeric, vectorized).}
}
\value{
Numeric vector of conditional variances (same length as \code{alpha}).
}
\description{
Computes \eqn{\mathrm{Var}(K_J \mid \alpha)} under a Dirichlet process prior.
}
\details{
The classical trigamma closed form is:
\deqn{v_J(\alpha) = \mu_J(\alpha) - \alpha^2\{\psi_1(\alpha) - \psi_1(\alpha+J)\}}
where \eqn{\psi_1(\cdot)} is the trigamma function.

This is equivalent to the direct summation:
\deqn{v_J(\alpha) = \sum_{i=1}^{J} \frac{\alpha(i-1)}{(\alpha + i - 1)^2}}

\strong{Key property:} \eqn{0 \leq v_J(\alpha) < \mu_J(\alpha)} for all
\eqn{\alpha > 0}. The lower inequality is strict when \code{J > 1}; at
\code{J = 1}, \eqn{K_J=1} deterministically and the variance is zero.

\strong{Limiting behavior:}
\itemize{
\item \eqn{\alpha \to 0^+}: \eqn{v_J(\alpha) \to 0}
\item \eqn{\alpha \to \infty}: \eqn{v_J(\alpha) \to 0}
}

The implementation evaluates each Bernoulli variance contribution as a
non-negative ratio term. It remains stable at both small and large
\code{alpha} without clipping cancellation artifacts to zero.
}
\examples{
var_K_given_alpha(50, 2.0)
var_K_given_alpha(50, c(0.5, 1, 2, 5))

# Verify underdispersion
J <- 50; alpha <- 2.0
mean_K_given_alpha(J, alpha) > var_K_given_alpha(J, alpha)  # TRUE

}
\seealso{
\code{\link{mean_K_given_alpha}}, \code{\link{moments_K_given_alpha}}

Other conditional_K:
\code{\link[=cdf_K_given_alpha]{cdf_K_given_alpha()}},
\code{\link[=cv_K_given_alpha]{cv_K_given_alpha()}},
\code{\link[=dispersion_K_given_alpha]{dispersion_K_given_alpha()}},
\code{\link[=log_pmf_K_given_alpha]{log_pmf_K_given_alpha()}},
\code{\link[=mean_K_given_alpha]{mean_K_given_alpha()}},
\code{\link[=mode_K_given_alpha]{mode_K_given_alpha()}},
\code{\link[=moments_K_given_alpha]{moments_K_given_alpha()}},
\code{\link[=pmf_K_given_alpha]{pmf_K_given_alpha()}},
\code{\link[=quantile_K_given_alpha]{quantile_K_given_alpha()}}
}
\concept{conditional_K}

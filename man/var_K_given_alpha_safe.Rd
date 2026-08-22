% Please edit documentation in R/07_jacobian.R
\name{var_K_given_alpha_safe}
\alias{var_K_given_alpha_safe}
\title{Conditional Variance of K Given Alpha (Enhanced)}
\usage{
var_K_given_alpha_safe(J, alpha)
}
\arguments{
\item{J}{Integer; sample size.}

\item{alpha}{Numeric vector; concentration parameter values.}
}
\value{
Numeric vector of conditional variances.
}
\description{
Compatibility wrapper for the stable conditional variance implementation.
}
\details{
Delegates to \code{var_K_given_alpha()}, whose positive-term implementation
is stable without clipping or a discontinuous small-alpha substitution.
}
\keyword{internal}

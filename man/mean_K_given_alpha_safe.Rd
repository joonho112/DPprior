% Please edit documentation in R/07_jacobian.R
\name{mean_K_given_alpha_safe}
\alias{mean_K_given_alpha_safe}
\title{Conditional Mean of K Given Alpha (Enhanced)}
\usage{
mean_K_given_alpha_safe(J, alpha)
}
\arguments{
\item{J}{Integer; sample size.}

\item{alpha}{Numeric vector; concentration parameter values.}
}
\value{
Numeric vector of conditional means.
}
\description{
Compatibility wrapper for the stable conditional mean implementation.
}
\details{
Delegates to \code{mean_K_given_alpha()}, whose finite-sum implementation is
stable without a discontinuous small-alpha substitution.
}
\keyword{internal}

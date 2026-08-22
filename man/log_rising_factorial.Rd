% Please edit documentation in R/04_pmf_conditional.R
\name{log_rising_factorial}
\alias{log_rising_factorial}
\title{Log Rising Factorial (Pochhammer Symbol)}
\usage{
log_rising_factorial(alpha, J)
}
\arguments{
\item{alpha}{Numeric; concentration parameter (must be positive scalar).}

\item{J}{Integer; number of terms (must be >= 1).}
}
\value{
Numeric; \eqn{\log(\alpha)_J}.
}
\description{
Computes the logarithm of the rising factorial (Pochhammer symbol)
\eqn{(\alpha)_J = \alpha(\alpha+1)\cdots(\alpha+J-1)}.
}
\details{
Although the identity
\deqn{\log(\alpha)_J = \log\Gamma(\alpha+J) - \log\Gamma(\alpha)}
is exact, subtracting the two log-gamma values loses accuracy when
\eqn{\alpha} is large relative to \eqn{J}. The implementation instead uses
finite sums with \code{log1p()} corrections, choosing the algebraically
equivalent branch whose ratios are small.
}
\examples{
\dontrun{
# (2)_3 = 2 * 3 * 4 = 24
exp(log_rising_factorial(2, 3))

# (1)_5 = 5!
exp(log_rising_factorial(1, 5))
factorial(5)

}
}
\keyword{internal}

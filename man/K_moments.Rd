% Please edit documentation in R/05_moments_marginal.R
\name{K_moments}
\alias{K_moments}
\title{Convenience Wrapper for Marginal Moments}
\usage{
K_moments(
  J,
  a,
  b,
  M = .QUAD_NODES_DEFAULT,
  M_verify = NULL,
  abs_tol = 1e-10,
  rel_tol = 1e-08,
  strict = FALSE
)
}
\arguments{
\item{J}{Integer; sample size (positive integer >= 1).}

\item{a}{Numeric; shape parameter of Gamma prior (> 0).}

\item{b}{Numeric; rate parameter of Gamma prior (> 0).}

\item{M}{Integer; number of quadrature nodes (default: 80).}

\item{M_verify}{Optional quadrature order satisfying the same verification
contract as \code{exact_K_moments()}: at least \code{max(2*M, M+40)} and
no greater than 512.}

\item{abs_tol, rel_tol}{Non-negative selected-versus-verification tolerances.}

\item{strict}{Logical; require successful higher-order verification.}
}
\value{
Named numeric vector \code{c(mean = ..., var = ...)} carrying a
\code{"marginal_metadata"} attribute with status, variance decomposition,
and quadrature verification details.
}
\description{
Returns marginal mean and variance as a named numeric vector.
}
\examples{
K_moments(50, 2.0, 1.0)

}
\seealso{
\code{\link{exact_K_moments}} for full output

Other marginal_K:
\code{\link[=cdf_K_marginal]{cdf_K_marginal()}},
\code{\link[=exact_K_moments]{exact_K_moments()}},
\code{\link[=mode_K_marginal]{mode_K_marginal()}},
\code{\link[=pmf_K_marginal]{pmf_K_marginal()}},
\code{\link[=quantile_K_marginal]{quantile_K_marginal()}},
\code{\link[=summary_K_marginal]{summary_K_marginal()}}
}
\concept{marginal_K}

% Please edit documentation in R/06_pmf_marginal.R
\name{cdf_K_marginal}
\alias{cdf_K_marginal}
\title{CDF of Marginal K Distribution}
\usage{
cdf_K_marginal(
  J,
  a,
  b,
  logS,
  M = .QUAD_NODES_DEFAULT,
  M_verify = NULL,
  abs_tol = 1e-10,
  rel_tol = 1e-08,
  strict = FALSE
)
}
\arguments{
\item{J}{Integer; sample size.}

\item{a}{Numeric; shape parameter of Gamma prior (> 0).}

\item{b}{Numeric; rate parameter of Gamma prior (> 0).}

\item{logS}{Matrix; pre-computed log-Stirling matrix.}

\item{M}{Integer; number of quadrature nodes (default: 80).}

\item{M_verify}{Optional independent order satisfying the package
verification rule (at least \code{max(2*M, M+40)}).}

\item{abs_tol, rel_tol}{Non-negative selected-versus-verification tolerances.}

\item{strict}{Logical; require successful higher-order verification.}
}
\value{
Numeric vector of length \eqn{J+1} containing
\eqn{P(K_J \leq k \mid a, b)} for \eqn{k = 0, 1, \ldots, J}.
The vector carries the PMF's \code{"marginal_metadata"} attribute.
}
\description{
Computes the cumulative distribution function \eqn{P(K_J \leq k \mid a, b)}
for \eqn{k = 0, 1, \ldots, J}.
}
\details{
The CDF satisfies:
\itemize{
\item \eqn{F(0) = 0} (since \eqn{P(K_J = 0) = 0})
\item \eqn{F(J) = 1}
\item \eqn{F(k)} is non-decreasing in \eqn{k}
}
}
\examples{
logS <- compute_log_stirling(50)
cdf <- cdf_K_marginal(50, 1.5, 0.5, logS)

# Verify CDF ends at 1
cdf[51]

# P(K <= 10)
cdf[11]

}
\seealso{
\code{\link{pmf_K_marginal}}, \code{\link{quantile_K_marginal}}

Other marginal_K:
\code{\link[=K_moments]{K_moments()}},
\code{\link[=exact_K_moments]{exact_K_moments()}},
\code{\link[=mode_K_marginal]{mode_K_marginal()}},
\code{\link[=pmf_K_marginal]{pmf_K_marginal()}},
\code{\link[=quantile_K_marginal]{quantile_K_marginal()}},
\code{\link[=summary_K_marginal]{summary_K_marginal()}}
}
\concept{marginal_K}

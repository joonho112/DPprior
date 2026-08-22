% Please edit documentation in R/06_pmf_marginal.R
\name{quantile_K_marginal}
\alias{quantile_K_marginal}
\title{Quantile of Marginal K Distribution}
\usage{
quantile_K_marginal(
  p,
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
\item{p}{Numeric; probability level(s) in \eqn{[0, 1]}. Can be scalar or vector.}

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
Integer vector of quantiles (same length as \code{p}).
Each element is the smallest \eqn{k} such that \eqn{P(K_J \leq k) \geq p}.
The vector carries the CDF's \code{"marginal_metadata"} attribute so an
approximate, unverified, or discrepant quadrature result is not silent.
}
\description{
Computes the \eqn{p}-th quantile of the marginal distribution of \eqn{K_J}.
}
\details{
This is the standard quantile definition for discrete distributions:
\eqn{Q(p) = \min\{k : F(k) \geq p\}}.
The lower endpoint follows the mathematical support: \eqn{Q(0)=1}, while
\eqn{Q(1)=J}. The structural \eqn{k=0} vector entry is never returned.

The function is vectorized over \code{p}, allowing efficient computation
of multiple quantiles in a single call.
When \code{M_verify} is supplied, the selected- and verification-order
integer quantiles are also compared exactly. A changed quantile downgrades
metadata to \code{"approximate"}; \code{strict = TRUE} raises a typed
convergence error. The selected-order quantile is always retained.
}
\examples{
logS <- compute_log_stirling(50)

# Single quantile (median)
quantile_K_marginal(0.5, 50, 1.5, 0.5, logS)

# Multiple quantiles at once
quantile_K_marginal(c(0.1, 0.25, 0.5, 0.75, 0.9), 50, 1.5, 0.5, logS)

# Interquartile range
qs <- quantile_K_marginal(c(0.25, 0.75), 50, 1.5, 0.5, logS)
diff(qs)

}
\seealso{
\code{\link{cdf_K_marginal}}, \code{\link{pmf_K_marginal}}

Other marginal_K:
\code{\link[=K_moments]{K_moments()}},
\code{\link[=cdf_K_marginal]{cdf_K_marginal()}},
\code{\link[=exact_K_moments]{exact_K_moments()}},
\code{\link[=mode_K_marginal]{mode_K_marginal()}},
\code{\link[=pmf_K_marginal]{pmf_K_marginal()}},
\code{\link[=summary_K_marginal]{summary_K_marginal()}}
}
\concept{marginal_K}

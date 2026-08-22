% Please edit documentation in R/06_pmf_marginal.R
\name{log_pmf_K_marginal}
\alias{log_pmf_K_marginal}
\title{Log Marginal PMF of K_J under Gamma Hyperprior}
\usage{
log_pmf_K_marginal(
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
\item{J}{Integer; sample size (must be >= 1).}

\item{a}{Numeric; shape parameter of Gamma prior (must be > 0).}

\item{b}{Numeric; rate parameter of Gamma prior (must be > 0).}

\item{logS}{Matrix; pre-computed log-Stirling matrix from
\code{\link{compute_log_stirling}}.}

\item{M}{Integer; number of quadrature nodes (default: \code{.QUAD_NODES_DEFAULT}).}

\item{M_verify}{Optional independent quadrature order satisfying the
package verification rule (at least \code{max(2*M, M+40)}, within the
supported ceiling). When supplied, the selected PMF is compared with a
fresh higher-order PMF.}

\item{abs_tol, rel_tol}{Non-negative absolute and relative tolerances for the
selected-versus-verification L1 discrepancy and the PMF-versus-direct-
moment checks at both quadrature orders.}

\item{strict}{Logical; require successful higher-order verification when
\code{TRUE}; otherwise return the selected PMF with explicit approximation
metadata.}
}
\value{
Numeric vector of length \eqn{J+1} containing log-probabilities
for \eqn{k = 0, 1, \ldots, J}. Entry \code{[1]} corresponds to \eqn{k=0}
and is always \code{-Inf}.
}
\description{
Computes \eqn{\log P(K_J = k \mid a, b)} for \eqn{k = 0, 1, \ldots, J}
using log-space mixing for numerical stability.
}
\details{
This routine normalizes \eqn{P(K_J = \cdot \mid \alpha_m)} at each quadrature
node before mixing, and then mixes in log-space via \code{logsumexp_vec}:
\deqn{\log p_k \approx \log\sum_m \exp\{\log w_m + \log p_{k\mid m}\}.}

The log-space computation is essential for numerical stability when:
\itemize{
\item J is large (tail probabilities become very small)
\item Alpha values span a wide range (extreme quadrature nodes)
\item Parameters lead to concentrated distributions
}

A verified result must pass both the selected-versus-verification PMF L1
budget and an independent cross-representation audit: at each order, the
mean and variance represented by the PMF must agree with the direct
marginal-moment quadrature. Thus two mutually agreeing but scientifically
incorrect PMFs cannot be classified as converged.
}
\keyword{internal}

% Please edit documentation in R/06_pmf_marginal.R
\name{mode_K_marginal}
\alias{mode_K_marginal}
\title{Mode of Marginal K Distribution}
\usage{
mode_K_marginal(
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
Integer; the value \eqn{k} that maximizes \eqn{P(K_J = k \mid a, b)},
carrying the PMF's \code{"marginal_metadata"} attribute.
}
\description{
Computes the mode (most likely value) of the marginal distribution of \eqn{K_J}.
}
\details{
The mode is always >= 1 since \eqn{P(K_J = 0) = 0}.
When \code{M_verify} is supplied, the selected and verification modes must
be identical for converged mode metadata. Disagreement is explicit and
never causes substitution of the verification-order mode.
}
\examples{
logS <- compute_log_stirling(50)
mode_K_marginal(50, 1.5, 0.5, logS)

}
\seealso{
\code{\link{pmf_K_marginal}}, \code{\link{summary_K_marginal}}

Other marginal_K:
\code{\link[=K_moments]{K_moments()}},
\code{\link[=cdf_K_marginal]{cdf_K_marginal()}},
\code{\link[=exact_K_moments]{exact_K_moments()}},
\code{\link[=pmf_K_marginal]{pmf_K_marginal()}},
\code{\link[=quantile_K_marginal]{quantile_K_marginal()}},
\code{\link[=summary_K_marginal]{summary_K_marginal()}}
}
\concept{marginal_K}

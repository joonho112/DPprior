% Please edit documentation in R/13_error_bounds.R
\name{DPprior_error_bounds}
\alias{DPprior_error_bounds}
\title{Compute A1 Approximation Error Bounds}
\usage{
DPprior_error_bounds(
  J,
  a,
  b,
  cJ = log(J),
  M = .QUAD_NODES_DEFAULT,
  M_verify = NULL,
  abs_tol = 1e-10,
  rel_tol = 1e-08,
  strict = FALSE
)
}
\arguments{
\item{J}{Integer; sample size.}

\item{a, b}{Numeric; Gamma hyperparameters (shape, rate).}

\item{cJ}{Numeric; scaling constant (default: log(J)).}

\item{M}{Integer; number of quadrature nodes (default: 80).}

\item{M_verify}{Optional independent quadrature order. When omitted and an
admissible order exists, the package-wide minimum
\code{max(2*M, M+40)} is used automatically.}

\item{abs_tol, rel_tol}{Non-negative selected-versus-verification tolerances.}

\item{strict}{Logical; require successful numerical verification of the
current exact moments, every adequacy-threshold scan point, and the
marginal TV bound.}
}
\value{
An S3 object of class "DPprior_error_bounds" with components:
\describe{
\item{J, a, b, cJ}{Input parameters}
\item{moment_errors}{List of moment error metrics from \code{a1_moment_error}}
\item{tv_bounds}{List with conditional and marginal TV bounds}
\item{recommendation}{"A1_sufficient" or "A2_recommended"}
\item{threshold_J}{First scan-grid J meeting the A1 adequacy criteria
(or NA), with verified/provisional status in \code{verification}}
\item{status, verified, usable}{Numerical verification state for this
error-analysis computation. These fields do not assert that A1 itself
is accurate; that is described by the reported errors and bounds.}
\item{verification, provenance}{Numerical checks, including threshold
scan provenance, and the theorem crosswalk.}
}
}
\description{
Comprehensive error analysis for the A1 large-J approximation.
Implements the error quantification framework from Lee (2026, Section 3.3).
}
\details{
The recommendation is based on:
\itemize{
\item A1 sufficient if: mean relative error < 5\\% AND variance relative error < 10\\%
\item A2 recommended otherwise
}

This function provides:
\enumerate{
\item Moment errors: Exact vs A1 approximation for mean and variance
\item TV bounds: Conditional bounds at multiple alpha values, plus marginal bound
\item Recommendation: Whether to use A1 or refine with A2
\item Threshold: Minimum J for A1 adequacy with current prior
}
}
\examples{
# Check A1 adequacy for J=50, typical prior
bounds <- DPprior_error_bounds(J = 50, a = 1.6, b = 1.2)
print(bounds)

# For larger J, A1 becomes more adequate
bounds_200 <- DPprior_error_bounds(J = 200, a = 1.6, b = 1.2)
print(bounds_200)

}
\references{
Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
\emph{arXiv preprint} arXiv:2602.06301.
}
\seealso{
\code{\link{a1_moment_error}}, \code{\link{compute_total_tv_bound}},
\code{\link{DPprior_a1}}, \code{\link{DPprior_a2_newton}}

Other diagnostics:
\code{\link[=DPprior_diagnostics]{DPprior_diagnostics()}},
\code{\link[=compute_error_landscape]{compute_error_landscape()}},
\code{\link[=compute_weight_diagnostics]{compute_weight_diagnostics()}}
}
\concept{diagnostics}

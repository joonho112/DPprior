% Please edit documentation in R/06_pmf_marginal.R
\name{verify_pmf_marginal_convergence}
\alias{verify_pmf_marginal_convergence}
\title{Report Marginal PMF Quadrature Discrepancies}
\usage{
verify_pmf_marginal_convergence(
  J,
  a,
  b,
  logS,
  M_values = c(20L, 40L, 80L, 120L),
  verbose = TRUE,
  abs_tol = 1e-10,
  rel_tol = 1e-08
)
}
\arguments{
\item{J}{Integer; sample size.}

\item{a}{Numeric; shape parameter of Gamma prior.}

\item{b}{Numeric; rate parameter of Gamma prior.}

\item{logS}{Matrix; pre-computed log-Stirling matrix.}

\item{M_values}{Integer vector; quadrature node counts to test.}

\item{verbose}{Logical; if TRUE, print results.}

\item{abs_tol, rel_tol}{Non-negative L1 error-budget components.}
}
\value{
Data frame of successive-order discrepancies, L1 budgets, and
pairwise budget classifications.
}
\description{
Reports successive-order PMF discrepancies and classifies each selected/
refined pair against an explicit L1 budget. Successive discrepancies are
not assumed to decrease monotonically.
}
\examples{
\dontrun{
logS <- compute_log_stirling(50)
verify_pmf_marginal_convergence(50, 1.5, 0.5, logS)

}
}
\keyword{internal}

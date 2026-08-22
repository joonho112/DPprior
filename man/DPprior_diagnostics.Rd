% Please edit documentation in R/14_diagnostics.R
\name{DPprior_diagnostics}
\alias{DPprior_diagnostics}
\title{Comprehensive Prior Diagnostics}
\usage{
DPprior_diagnostics(
  fit,
  thresholds = c(0.5, 0.9),
  warning_policy = NULL,
  M_verify = NULL,
  abs_tol = 1e-10,
  rel_tol = 1e-08,
  allow_approximate = FALSE
)
}
\arguments{
\item{fit}{A canonical \code{dpprior.result/1} fit. Legacy and flat lists
are rejected; call \code{upgrade_DPprior_object()} explicitly first.}

\item{thresholds}{Reserved fixed summary-view thresholds. The only accepted
value is \code{c(0.5, 0.9)}. Use \code{prob_wsb_exceeds()} for arbitrary
W_SB tails or \code{wmax_tail_bounds()} for certified W_max bounds.}

\item{warning_policy}{Optional ordinary list with exact fields
\code{estimand}, \code{direction}, \code{weight_threshold}, and
\code{action_threshold}. W_max policies are retained as indeterminate
because canonical diagnostics make no W_max point-probability claim.}

\item{M_verify}{Optional independent quadrature order.}

\item{abs_tol, rel_tol}{Comparison tolerances.}

\item{allow_approximate}{Return an approximate bundle when TRUE; otherwise
signal a typed condition carrying the complete bundle.}
}
\value{
A \code{DPprior_diagnostics} object using the shared result schema.
}
\description{
Comprehensive Prior Diagnostics
}
\seealso{
Other diagnostics:
\code{\link[=DPprior_error_bounds]{DPprior_error_bounds()}},
\code{\link[=compute_error_landscape]{compute_error_landscape()}},
\code{\link[=compute_weight_diagnostics]{compute_weight_diagnostics()}}
}
\concept{diagnostics}

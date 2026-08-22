% Please edit documentation in R/14_diagnostics.R
\name{compute_weight_diagnostics}
\alias{compute_weight_diagnostics}
\title{First size-biased and largest-weight diagnostics}
\usage{
compute_weight_diagnostics(
  a,
  b,
  thresholds = c(0.3, 0.5, 0.7, 0.9),
  M = .QUAD_NODES_DEFAULT,
  M_verify = NULL,
  abs_tol = 1e-10,
  rel_tol = 1e-08
)
}
\arguments{
\item{a, b}{Numeric Gamma hyperparameters.}

\item{thresholds}{Numeric vector in (0,1).}

\item{M}{Selected quadrature order for E(W_SB).}

\item{M_verify}{Optional verification order.}

\item{abs_tol, rel_tol}{Comparison tolerances.}
}
\value{
Separate \code{size_biased} and \code{maximum} components. Numeric
compatibility aliases refer only to W_SB; no risk category is returned.
}
\description{
First size-biased and largest-weight diagnostics
}
\seealso{
Other diagnostics:
\code{\link[=DPprior_diagnostics]{DPprior_diagnostics()}},
\code{\link[=DPprior_error_bounds]{DPprior_error_bounds()}},
\code{\link[=compute_error_landscape]{compute_error_landscape()}}
}
\concept{diagnostics}

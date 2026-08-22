% Please edit documentation in R/20_elicitation_sensitivity.R
\name{.dpprior_run_elicitation_sensitivity}
\alias{.dpprior_run_elicitation_sensitivity}
\title{Run a canonical elicitation-sensitivity grid}
\usage{
.dpprior_run_elicitation_sensitivity(
  J,
  scenarios,
  method = NULL,
  M = .QUAD_NODES_DEFAULT,
  check_diagnostics = TRUE,
  weight_thresholds = c(0.5, 0.9),
  .fit_fun = DPprior_fit,
  .key_fun = .dp_sensitivity_content_key
)
}
\arguments{
\item{J}{Common positive integer sample size.}

\item{scenarios}{Non-empty list of named scenario lists, or a data frame.}

\item{method}{Optional common A2-MN, A1, or A2-KL method.}

\item{M}{Common selected quadrature order in 10:256.}

\item{check_diagnostics}{Whether fresh diagnostic metrics are requested.}

\item{weight_thresholds}{Fixed c(0.5, 0.9) output thresholds.}

\item{.fit_fun}{Internal canonical fit adapter.}

\item{.key_fun}{Internal deterministic scenario-key adapter.}
}
\value{
A canonical dpprior.result/1 sensitivity object.
}
\description{
This internal API normalizes all scenarios before the first calibration,
evaluates actual \code{DPprior_fit()} routes, and returns one canonical
\code{dpprior.result/1} sensitivity object.  It describes prior-input
sensitivity only; posterior robustness requires a simulation study.
}
\keyword{internal}

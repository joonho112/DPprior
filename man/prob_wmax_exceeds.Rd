% Please edit documentation in R/08_weights_w1.R
\name{prob_wmax_exceeds}
\alias{prob_wmax_exceeds}
\title{Largest DP Population Weight Tail Probability}
\usage{
prob_wmax_exceeds(
  threshold,
  alpha = NULL,
  a = NULL,
  b = NULL,
  method = c("auto", "deterministic", "bounds", "monte_carlo"),
  rel_tol = 1e-10,
  abs_tol = 1e-12,
  subdivisions = 1000L,
  n = 250000L,
  seed = NULL,
  conf_level = 0.95,
  max_sticks = 100000L,
  warn_low_successes = TRUE
)
}
\arguments{
\item{threshold}{One finite probability in the closed unit interval.}

\item{alpha}{Optional fixed positive DP concentration.}

\item{a, b}{Optional positive Gamma shape and rate for alpha.}

\item{method}{One of \code{"auto"}, \code{"deterministic"},
\code{"bounds"}, or \code{"monte_carlo"}.}

\item{rel_tol, abs_tol}{Positive numerical tolerances. Deterministic public
results require \code{rel_tol <= 1e-10} and \code{abs_tol <= 1e-12}.}

\item{subdivisions}{Positive integration subdivision limit.}

\item{n}{Monte Carlo draw count; the default supports a worst-case normal
95-percent half-width of about 0.002.}

\item{seed}{Required non-negative integer for Monte Carlo.}

\item{conf_level}{Monte Carlo confidence level in the open unit interval.}

\item{max_sticks}{Maximum GEM iterations. Reaching it produces an explicit
failed result; an unresolved remainder is never discarded.}

\item{warn_low_successes}{Whether to warn when fewer than 25 exceedances
are observed.}
}
\value{
A \code{wmax_tail_result} with the named estimand, conditioning,
direct estimate (when available), direct numerical/sampling interval,
absolute-error statement, universal certified bounds, status,
verification flag, method, and provenance. A bounds-only result has
\code{estimate = NA} and \code{usable = FALSE}; neither certified endpoint
is used as an estimate. Its certified interval remains available with
\code{bounds_usable = TRUE}, or directly from \code{wmax_tail_bounds()}.
}
\description{
Computes or bounds \eqn{P(W_{max}>threshold)} without conflating it with
the first size-biased tail \eqn{P(W_{SB}>threshold)}. For thresholds at
least 0.5, \code{method = "auto"} uses a stable one-dimensional integral
and an independent calculation. Below 0.5, deterministic direct evaluation
is intentionally deferred and \code{"auto"} returns certified bounds only.
Explicit seeded Monte Carlo is available at every threshold and uses an
exact GEM stopping rule, so its only approximation is sampling error.
}
\examples{
prob_wmax_exceeds(0.5, alpha = 1) # log(2)
prob_wmax_exceeds(0.9, a = 2, b = 1)
prob_wmax_exceeds(0.4, alpha = 1) # certified bounds only

}
\seealso{
\code{\link{wmax_tail_bounds}},
\code{\link{prob_wsb_exceeds}}

Other weights_w1:
\code{\link[=cdf_w1]{cdf_w1()}},
\code{\link[=density_w1]{density_w1()}},
\code{\link[=mean_w1]{mean_w1()}},
\code{\link[=prob_w1_exceeds]{prob_w1_exceeds()}},
\code{\link[=prob_wsb_exceeds]{prob_wsb_exceeds()}},
\code{\link[=quantile_w1]{quantile_w1()}},
\code{\link[=summary_w1]{summary_w1()}},
\code{\link[=var_w1]{var_w1()}},
\code{\link[=wmax_tail_bounds]{wmax_tail_bounds()}}
}
\concept{weights_w1}

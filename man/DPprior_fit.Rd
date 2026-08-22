% Please edit documentation in R/16_DPprior_fit.R
\name{DPprior_fit}
\alias{DPprior_fit}
\title{Fit a Gamma hyperprior from an occupied-cluster target}
\usage{
DPprior_fit(
  J,
  mu_K = NULL,
  var_K = NULL,
  confidence = c("medium", "low", "high"),
  method = c("A2-MN", "A1", "A2-KL"),
  target_pmf = NULL,
  check_diagnostics = TRUE,
  warn_dominance = NULL,
  M = .QUAD_NODES_DEFAULT,
  verbose = FALSE,
  warning_policy = NULL,
  a1_projection = c("error", "nearest"),
  cv_K = NULL,
  K_interval = NULL,
  target_K = NULL,
  ...
)
}
\arguments{
\item{J}{Integer design size. When \code{target_K} is supplied, \code{J}
must match the target exactly.}

\item{mu_K, var_K}{Optional direct target mean and variance.}

\item{confidence}{Optional \code{"low"}, \code{"medium"}, or
\code{"high"} qualitative uncertainty route.}

\item{method}{One of \code{"A2-MN"}, \code{"A1"}, or \code{"A2-KL"}.}

\item{target_pmf}{Optional strict target PMF on \eqn{1,\ldots,J}.}

\item{check_diagnostics}{Whether to attach freshly verified canonical
diagnostics.}

\item{warn_dominance}{Deprecated first-size-biased-weight warning adapter.}

\item{M}{Selected quadrature order.}

\item{verbose}{Whether the numerical producer prints progress.}

\item{warning_policy}{Optional ordinary list with exact canonical fields
\code{estimand}, \code{direction}, \code{weight_threshold}, and
\code{action_threshold}. The legacy \code{threshold} field is accepted
here and normalized immediately for compatibility.}

\item{a1_projection}{A1 projection policy: \code{"error"} or an explicit
opt-in to \code{"nearest"}.}

\item{cv_K}{Optional coefficient of variation \eqn{SD(K_J)/E(K_J)}.}

\item{K_interval}{Optional canonical interval-target specification.}

\item{target_K}{Optional validated \code{dpprior.target/1} K target. It
cannot be combined with scalar target arguments.}

\item{...}{Must be empty.}
}
\value{
A canonical \code{dpprior.result/1} \code{DPprior_fit} object.
Use \code{parameters$a} and \code{parameters$b} for the Gamma shape and
rate; \code{target$K}, \code{achieved}, and \code{residuals} for the
scientific claim; and \code{computation}, \code{verification}, and
\code{provenance} for its audit trail. A result is decision-ready only
when its mode-specific contract and \code{status}, \code{usable}, and
\code{verified} fields permit that use. If a backend signals a typed
calibration condition, the canonical target or result is retained in
\code{condition$result}.
}
\description{
\code{DPprior_fit()} constructs one canonical \eqn{K_J} target, or accepts
one from \code{\link{DPprior_target_K}}, dispatches to A1, A2-MN, or A2-KL,
and returns the numerical producer's native canonical result without
flattening or reinterpreting its evidence.
}
\details{
Supply exactly one uncertainty route: direct \code{var_K},
\code{confidence}, \code{cv_K}, \code{K_interval}, or
\code{target_pmf}. A pre-built \code{target_K} is an alternative to all
scalar target arguments. Compatibility aliases are non-authoritative; new
code should use the canonical nested fields.
}
\examples{
target <- DPprior_target_K(J = 50, mu_K = 5, var_K = 8)
fit <- DPprior_fit(
  J = 50, target_K = target, method = "A2-MN",
  check_diagnostics = FALSE
)
fit$parameters
fit[c("status", "usable", "verified")]
}
\seealso{
Other elicitation:
\code{\link[=DPprior_a1]{DPprior_a1()}},
\code{\link[=DPprior_a2_kl]{DPprior_a2_kl()}},
\code{\link[=DPprior_a2_newton]{DPprior_a2_newton()}},
\code{\link[=DPprior_dual]{DPprior_dual()}},
\code{\link[=DPprior_dual_hard]{DPprior_dual_hard()}},
\code{\link[=DPprior_dual_soft]{DPprior_dual_soft()}},
\code{\link[=DPprior_target_K]{DPprior_target_K()}}
}
\concept{elicitation}

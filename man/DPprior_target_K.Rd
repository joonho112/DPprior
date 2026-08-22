% Please edit documentation in R/19_elicitation_targets.R
\name{DPprior_target_K}
\alias{DPprior_target_K}
\title{Construct a Canonical Bounded-Discrete Cluster-Count Target}
\usage{
DPprior_target_K(
  J,
  mu_K = NULL,
  var_K = NULL,
  confidence = NULL,
  cv_K = NULL,
  K_interval = NULL,
  target_pmf = NULL,
  tolerance = 1e-09,
  root_control = NULL
)
}
\arguments{
\item{J}{Integer design size and upper support point for \eqn{K_J}.}

\item{mu_K}{Optional target mean for \eqn{K_J}.}

\item{var_K}{Optional target variance for \eqn{K_J}.}

\item{confidence}{Optional qualitative confidence label.}

\item{cv_K}{Optional coefficient of variation, defined as
\eqn{SD(K_J)/E(K_J)}.}

\item{K_interval}{Optional explicit interval-target specification.}

\item{target_pmf}{Optional strict PMF on the support \eqn{1,\ldots,J}; a
structural zero for \eqn{K=0} may be supplied as the first of \code{J+1}
entries.}

\item{tolerance}{Positive numerical postcondition tolerance.}

\item{root_control}{Optional named controls for interval-target root solving.}
}
\value{
A validated \code{dpprior.target/1} \code{dpprior_K_target}
object. The exact support and original request are retained alongside
normalized and used values, derivation evidence, implied moments,
tolerances, verification, and provenance. Pass the object unchanged as
\code{target_K} to \code{\link{DPprior_fit}}.
}
\description{
Constructs and validates a canonical elicitation target for the occupied
cluster count \eqn{K_J} on the exact support \eqn{1,\ldots,J}. Supply exactly
one uncertainty source: \code{var_K}, \code{confidence}, \code{cv_K},
\code{K_interval}, or \code{target_pmf}.
}
\examples{
target <- DPprior_target_K(J = 50, mu_K = 5, var_K = 8)
target$used[c("mu_K", "var_K")]
fit <- DPprior_fit(
  J = 50, target_K = target, check_diagnostics = FALSE
)

}
\seealso{
Other elicitation:
\code{\link[=DPprior_a1]{DPprior_a1()}},
\code{\link[=DPprior_a2_kl]{DPprior_a2_kl()}},
\code{\link[=DPprior_a2_newton]{DPprior_a2_newton()}},
\code{\link[=DPprior_dual]{DPprior_dual()}},
\code{\link[=DPprior_dual_hard]{DPprior_dual_hard()}},
\code{\link[=DPprior_dual_soft]{DPprior_dual_soft()}},
\code{\link[=DPprior_fit]{DPprior_fit()}}
}
\concept{elicitation}

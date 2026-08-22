% Please edit documentation in R/10_a1_mapping.R
\name{DPprior_a1}
\alias{DPprior_a1}
\title{A1 Closed-Form Prior Elicitation}
\usage{
DPprior_a1(
  J,
  mu_K,
  var_K,
  scaling = c("log", "harmonic", "digamma"),
  epsilon = .TOL_PROJECTION_BUFFER,
  projection = c("error", "nearest")
)
}
\arguments{
\item{J}{Integer; number of items/sites (must be >= 2).}

\item{mu_K}{Numeric; target prior mean of \eqn{K_J} (must satisfy \eqn{1 < \mu_K < J}).}

\item{var_K}{Numeric; target prior variance of \eqn{K_J} (must be > 0).}

\item{scaling}{Character; scaling constant method: "log" (default),
"harmonic", or "digamma".}

\item{epsilon}{Numeric; preferred numerical-interior buffer and near-boundary
tolerance. The buffer is adaptively capped in narrow feasible domains and
floored at local floating-point spacing. A target is changed only when an
explicit projection is requested. Default is
\code{.TOL_PROJECTION_BUFFER} (1e-6).}

\item{projection}{Character; A1 target-projection policy. The default
\code{"error"} never changes the requested target. Set
\code{"nearest"} to explicitly opt into projection of an A1-infeasible
variance to the nearest numerically interior A1 target. The original and
projected targets, distance, and policy are retained in the result.}
}
\value{
A canonical \code{dpprior.result/1} \code{DPprior_fit} object.
Gamma shape and rate are in \code{parameters}; the immutable canonical
target is in \code{target$K}; proxy achievements and residuals are in
\code{achieved} and \code{residuals}. The \code{proxy},
\code{computation}, \code{verification}, and \code{provenance} records
retain the closed-form formula, scaling, projection, and algebraic
round-trip evidence. A1 is always labelled \code{status = "approximate"},
\code{usable = TRUE}, and \code{verified = FALSE}: its proxy identity is
not finite-design verification. Registered flat aliases are migration
views and are not authoritative scientific fields.
}
\description{
Maps target beliefs about the number of clusters \eqn{(\mu_K, \sigma^2_K)}
to Gamma hyperprior parameters \eqn{(a, b)} using the A1 closed-form
approximation based on Negative Binomial moment matching.
}
\details{
\subsection{Theory (TSMM Stage 1)}{

The A1 method uses a shifted Negative Binomial approximation:
\deqn{K_J - 1 \mid \alpha \approx \text{Poisson}(\alpha \cdot c_J)}

With \eqn{\alpha \sim \text{Gamma}(a, b)}, the marginal becomes:
\deqn{K_J - 1 \approx \text{NegBin}(a, b/(b + c_J))}
}

\subsection{Inverse Formulas (Theorem 1)}{

Let \eqn{m = \mu_K - 1} (shifted mean) and \eqn{D = \sigma^2_K - m}.
If \eqn{D > 0} (overdispersion):
\deqn{a = m^2 / D, \quad b = m \cdot c_J / D}

If \eqn{D \leq 0}, the target is infeasible for A1. The default is a typed
error. Projection occurs only with \code{projection = "nearest"} and emits
exactly one typed projection warning.
}

\subsection{Feasibility}{

The NegBin model requires overdispersion: \eqn{\sigma^2_K > \mu_K - 1}.
High-confidence specifications (low variance) may violate this constraint
under the A1 proxy, even though they may be feasible under the exact DP.
}
}
\examples{
# Basic usage with moment targets
fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)
print(fit)

# Using VIF specification
fit <- DPprior_a1(J = 50, mu_K = 5, var_K = vif_to_variance(5, 2))

# Using confidence-based specification
vif <- confidence_to_vif("medium")
fit <- DPprior_a1(J = 50, mu_K = 5, var_K = vif_to_variance(5, vif))

# Explicit projection of an A1-infeasible variance
fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 3,
                  projection = "nearest")

# Compare scaling methods
fit_log <- DPprior_a1(50, 5, 8, scaling = "log")
fit_harm <- DPprior_a1(50, 5, 8, scaling = "harmonic")

}
\references{
Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
\emph{arXiv preprint} arXiv:2602.06301.
}
\seealso{
\code{\link{vif_to_variance}} for VIF conversion,
\code{\link{confidence_to_vif}} for confidence mapping,
\code{\link{print.DPprior_fit}} for print method

Other elicitation:
\code{\link[=DPprior_a2_kl]{DPprior_a2_kl()}},
\code{\link[=DPprior_a2_newton]{DPprior_a2_newton()}},
\code{\link[=DPprior_dual]{DPprior_dual()}},
\code{\link[=DPprior_dual_hard]{DPprior_dual_hard()}},
\code{\link[=DPprior_dual_soft]{DPprior_dual_soft()}},
\code{\link[=DPprior_fit]{DPprior_fit()}},
\code{\link[=DPprior_target_K]{DPprior_target_K()}}
}
\concept{elicitation}

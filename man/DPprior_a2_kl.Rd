% Please edit documentation in R/12_a2_kl.R
\name{DPprior_a2_kl}
\alias{DPprior_a2_kl}
\title{A2-KL: KL Divergence Minimization for Prior Calibration}
\usage{
DPprior_a2_kl(
  J,
  target,
  method = c("pmf", "chisq"),
  max_iter = 100L,
  tol = 1e-06,
  M = .QUAD_NODES_DEFAULT,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{J}{Integer; sample size (number of observations). Must be >= 2.}

\item{target}{Either:
\itemize{
\item Numeric vector of length J or J+1: target PMF for k = 1, ..., J,
or for k = 0, ..., J with the k = 0 entry zero and dropped
(used when \code{method = "pmf"}).
\item Named list with \code{mu_K} and \code{var_K}: construct from moments
using a discretized scaled chi-square (used when \code{method = "chisq"}).
}}

\item{method}{Character; \code{"pmf"} or \code{"chisq"}. Default: \code{"pmf"}.}

\item{max_iter}{Integer; maximum optimization iterations. Default: 100.}

\item{tol}{Numeric; convergence tolerance for optimization. Default: 1e-6.}

\item{M}{Integer; number of quadrature nodes. Default: 80.}

\item{verbose}{Logical; if TRUE, print optimization progress. Default: FALSE.}

\item{...}{Optional tuning parameters:
\itemize{
\item \code{log_bounds}: numeric length-2 vector giving lower/upper bounds
on \code{log(a)} and \code{log(b)} (default \code{c(-15, 15)}).
\item \code{M_verify}: compliant higher quadrature order (by default
\code{max(2*M, M+40)} when supported).
\item \code{kl_tol}, \code{l1_tol}, \code{mean_scaled_tol}, and
\code{var_scaled_tol}: independent adequacy tolerances.
\item \code{pmf_abs_tol}, \code{pmf_rel_tol}: higher-order PMF
verification tolerances.
\item \code{boundary_tol}: log-parameter bound proximity tolerance.
\item \code{fallback_max_iter}: iteration budget for the predeclared
bounded \code{nlminb} fallback.
}}
}
\value{
A canonical \code{dpprior.result/1} \code{DPprior_fit} object. The
authoritative fields are:
\describe{
\item{\code{parameters}}{Selected Gamma shape \code{a}, rate \code{b},
and the \code{log_ab} parameterization.}
\item{\code{J}}{Integer; sample size}
\item{\code{target$K}}{Canonical requested, normalized, and used target.
A strict PMF remains authoritative. For \code{method = "chisq"}, the
requested moments remain authoritative and the conditioned chi-square
PMF is recorded only as \code{A2_KL_objective} derivation evidence.}
\item{\code{achieved}, \code{residuals}}{Selected-order induced PMF,
moments, and distinct KL, L1, mean, and variance residuals.}
\item{\code{computation}}{Closed controls, quadrature orders, trace,
typed initializer/optimizer attempts, fresh candidate evaluations,
fallback, candidate selection, and termination evidence.}
\item{\code{verification}}{A fixed-parameter verifier snapshot at the
required independent order, four distribution-adequacy gates, and a
separate PMF order-stability check.}
\item{\code{method}}{Character; \code{"A2-KL"}.}
\item{\code{status}}{One of \code{converged}, \code{boundary},
\code{approximate}, \code{infeasible}, or \code{failed}.}
\item{\code{usable}, \code{verified}}{Logical adequacy flags.}
\item{\code{provenance}}{Method, backend, approximation, projection,
and migration records.}
\item{\code{compatibility}}{A quarantined non-authoritative Phase 8
view. Top-level \code{a}, \code{b}, \code{fit}, \code{diagnostics},
\code{attempts}, and \code{trace} are registered migration aliases.}
}
}
\description{
Calibrates Gamma hyperprior parameters \eqn{(a, b)} by minimizing the
Kullback-Leibler divergence between a target distribution and the induced
marginal PMF of the number of clusters \eqn{K_J}.
}
\details{
The A2-KL algorithm finds \eqn{(a^*, b^*)} by solving:
\deqn{(a^*, b^*) = \arg\min_{a,b} D_{KL}(p^*(K_J) \| p_{a,b}(K_J))}

\strong{Algorithm:}
\enumerate{
\item Validate a normalized custom PMF, or explicitly construct and
support-condition a chi-square target PMF
\item Initialize \eqn{(a_0, b_0)} from A2-MN (exact moment matching)
\item Optimize KL divergence using L-BFGS-B in log-space with bounds
\item If the primary optimizer fails, try a recorded bounded
\code{nlminb} fallback for the same exact objective
\item Independently reconstruct target and induced PMFs at a compliant
higher quadrature order and assess KL, L1, and scaled moment errors
\item Classify status from verification and adequacy, not optimizer exit
}

\strong{Stability features:}
\itemize{
\item Log-parameterization ensures positivity of (a, b)
\item Explicit bounds prevent numerical overflow/underflow
\item Exact zero-support KL semantics without epsilon smoothing
\item Explicit fallback attempts; initialization is never labelled success
\item A2-MN initialization provides good starting point
}

\strong{When to use A2-KL vs A2-MN:}
\itemize{
\item \strong{A2-MN}: Target moments only (exact moment matching)
\item \strong{A2-KL}: Full target distribution shape (multi-modal,
skewed, expert-elicited)
}

The canonical \code{achieved} and \code{residuals} fields always contain the
selected-\code{M} computation. Higher-order PMF values are retained in the
independent \code{verification$verifier_snapshot}; they determine status but
never replace the selected-order result. The default adequacy limits are KL
0.015, L1 0.11, scaled mean 0.01, and scaled variance 0.065. The PMF order
stability tolerance is a separate numerical check and is never substituted
for any of those four scientific adequacy gates.
}
\examples{
# Example 1: Target from moments (method = "chisq")
fit <- DPprior_a2_kl(J = 50, target = list(mu_K = 5, var_K = 8),
                     method = "chisq")
print(fit)

# Example 2: Custom target PMF (method = "pmf")
target_pmf <- dbinom(1:50, size = 50, prob = 0.1)
target_pmf <- target_pmf / sum(target_pmf)
fit2 <- DPprior_a2_kl(J = 50, target = target_pmf, method = "pmf")

# Example 3: Compare A2-KL vs A2-MN
a2_mn <- DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8)
a2_kl <- DPprior_a2_kl(J = 50, target = list(mu_K = 5, var_K = 8),
                       method = "chisq")
cat(sprintf(
  "A2-MN: a=\%.4f, b=\%.4f\n",
  a2_mn$parameters$a, a2_mn$parameters$b
))
cat(sprintf(
  "A2-KL: a=\%.4f, b=\%.4f, KL=\%.4e\n",
  a2_kl$parameters$a, a2_kl$parameters$b,
  a2_kl$residuals$distribution$kl
))

}
\references{
Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
\emph{arXiv preprint} arXiv:2602.06301.
}
\seealso{
\code{\link{DPprior_a2_newton}} for exact moment matching,
\code{\link{DPprior_a1}} for closed-form initialization,
\code{\link{kl_divergence_K}} for KL divergence computation,
\code{\link{discretize_chisq}} for chi-square discretization

Other elicitation:
\code{\link[=DPprior_a1]{DPprior_a1()}},
\code{\link[=DPprior_a2_newton]{DPprior_a2_newton()}},
\code{\link[=DPprior_dual]{DPprior_dual()}},
\code{\link[=DPprior_dual_hard]{DPprior_dual_hard()}},
\code{\link[=DPprior_dual_soft]{DPprior_dual_soft()}},
\code{\link[=DPprior_fit]{DPprior_fit()}},
\code{\link[=DPprior_target_K]{DPprior_target_K()}}
}
\concept{elicitation}

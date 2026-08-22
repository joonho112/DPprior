% Please edit documentation in R/11_a2_newton.R
\name{DPprior_a2_newton}
\alias{DPprior_a2_newton}
\title{A2-MN Calibrated-Moment Newton Solver}
\usage{
DPprior_a2_newton(
  J,
  mu_K,
  var_K,
  a0 = NULL,
  b0 = NULL,
  tol_F = .TOL_NEWTON,
  tol_rel = 1e-08,
  tol_step = 1e-10,
  max_iter = 20L,
  damping = TRUE,
  use_fallback = TRUE,
  M = .QUAD_NODES_DEFAULT,
  M_verify = NULL,
  verification_abs_tol = 1e-10,
  verification_rel_tol = 1e-08,
  verbose = FALSE
)
}
\arguments{
\item{J}{Integer; sample size (number of observations/sites). Must be >= 2.}

\item{mu_K}{Numeric; target prior mean \eqn{E[K_J]}. Must satisfy \eqn{1 < \mu_K < J}.}

\item{var_K}{Numeric; target prior variance \eqn{\mathrm{Var}(K_J)}. Must be positive.}

\item{a0}{Numeric or NULL; initial shape parameter. If NULL, computed via
\code{\link{DPprior_a1}}.}

\item{b0}{Numeric or NULL; initial rate parameter. If NULL, computed via
\code{\link{DPprior_a1}}.}

\item{tol_F}{Positive numeric; absolute part of the componentwise moment
tolerance. Retained under its historical name for compatibility.}

\item{tol_rel}{Non-negative numeric; relative part of the componentwise
tolerance. Each component uses
\code{tol_F + tol_rel * max(abs(target), 1)}.}

\item{tol_step}{Numeric; stopping tolerance for Newton step size.
Default: 1e-10.}

\item{max_iter}{Integer; maximum Newton iterations. Default: 20.}

\item{damping}{Logical; if TRUE, use backtracking line search for
damped Newton updates. Default: TRUE.}

\item{use_fallback}{Logical; if TRUE, use Nelder-Mead fallback when Newton
fails to converge. Default: TRUE.}

\item{M}{Integer; number of quadrature nodes for moment computation.
Default: 80.}

\item{M_verify}{Optional integer; independent verification order. It must be
at least \code{max(2*M, M+40)} and no greater than 512. The required order
is selected automatically when available. If the required order exceeds
512, the function fails closed with a typed no-candidate condition. The
numerical source candidate is then available only in the condition's
non-authoritative compatibility view, with guidance to refit using
\code{M <= 256}.}

\item{verification_abs_tol, verification_rel_tol}{Non-negative absolute and
relative tolerances for selected-versus-verification order stability. At
least one must be positive.}

\item{verbose}{Logical; if TRUE, print iteration progress. Default: FALSE.}
}
\value{
A canonical \code{dpprior.result/1} \code{DPprior_fit} object with
authoritative nested components:
\describe{
\item{\code{parameters}}{Gamma shape, rate, and parameterization.}
\item{\code{J}}{Integer; sample size}
\item{\code{target}}{Canonical immutable \code{K} target.}
\item{\code{achieved},\code{residuals}}{Selected-order public
moment evidence.}
\item{\code{method}}{Character; "A2-MN" or "A2-MN+NM" if fallback was used}
\item{\code{status}}{One of \code{"converged"}, \code{"boundary"},
\code{"approximate"}, or \code{"failed"}.}
\item{\code{usable},\code{verified}}{Logical status fields derived from
the numerical and independent-verification gates.}
\item{\code{computation}}{Orders, fixed scaling, typed attempts,
candidate ledger, fallback, termination, and trace.}
\item{\code{verification}}{Independent higher-order recomputation,
componentwise adequacy, and order-stability evidence.}
\item{\code{provenance},\code{compatibility}}{Canonical provenance
plus a quarantined non-authoritative legacy view.}
}
Finite results also retain exact compatibility aliases \code{a}, \code{b},
\code{converged}, \code{iterations}, \code{termination}, \code{fit},
\code{diagnostics}, \code{trace}, and \code{attempts}. If no independent
verifier order is available, the function signals a
\code{dpprior_a2_no_candidate} condition containing a valid parameterless
failed result in its \code{result} field.
}
\description{
Seeks Gamma(a, b) hyperprior parameters that match target moments for the
number of clusters K_J under a Dirichlet process prior, and reports
convergence only after an independent higher-order check.
}
\details{
This implements TSMM Stage 2 (A2-MN) from Lee (2026).
The A2-MN algorithm uses Newton's method in log-scale to ensure positivity
of the Gamma parameters. The Jacobian uses score-function identities and
their order-refinement diagnostics. Residual components are scaled by
\code{max(abs(target), 1)} and checked individually under an explicit
absolute-plus-relative budget; a small aggregate norm cannot hide a failed
component.

\strong{Algorithm Steps:}
\enumerate{
\item Initialize: \eqn{(a_0, b_0)} from A1 closed-form or user-provided
\item Log-parameterize: \eqn{\eta = (\log a, \log b)}
\item For each iteration:
\itemize{
\item Compute moments \eqn{(M_1, V)} and Jacobian \eqn{J_F}
\item Compute natural-unit and scaled component residuals
\item Transform Jacobian to log-scale: \eqn{J_{\log} = J_F \cdot \text{diag}(a, b)}
\item Newton step: \eqn{\Delta = -J_{\log}^{-1} F}
\item Backtracking line search (if damping enabled)
\item Update: \eqn{\eta \leftarrow \eta + \lambda \Delta}
}
\item If needed, run the predeclared Nelder--Mead fallback on
\eqn{(\log a,\log b)}, never on raw positive parameters
\item Independently recompute the selected candidate at a distinct
higher quadrature order; retain the fit-order candidate unchanged
}

\strong{Termination and status:}
\itemize{
\item Optimizer exit, a small step, or line-search exhaustion is never by
itself convergence.
\item \code{converged} requires the selected- and higher-order residual
gates plus the quadrature-stability gate.
\item \code{approximate} retains a finite but unusable candidate when a
solver or available independent-verification gate is unmet.
\item An unavailable required verifier order produces a typed
\code{failed}/no-candidate result; it never exposes the numerical source
candidate as public scientific output.
\item Singular/ill-conditioned Jacobians, line-search failure,
stagnation, fallback outcome, and verification failure have stable
reason codes in \code{computation} and the compatibility audit view.
}

A1 projection may be requested explicitly to obtain an initialization only.
Its projection record is retained in the non-authoritative initialization
audit view, while the canonical requested A2 target is immutable and is
never silently projected.
}
\examples{
# Basic usage
fit <- DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8)
print(fit)

# Inspect the selected-order match and independent verification
achieved <- exact_K_moments(
  50, fit$parameters$a, fit$parameters$b
)
cat(sprintf("Target E[K]=5, Achieved E[K]=\%.10f\n", achieved$mean))
cat(sprintf("Target Var=8, Achieved Var=\%.10f\n", achieved$var))

# Compare A1 vs A2 accuracy
a1 <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)
a1_mom <- exact_K_moments(
  50, a1$parameters$a, a1$parameters$b
)
a2_mom <- exact_K_moments(
  50, fit$parameters$a, fit$parameters$b
)
cat(sprintf("A1 mean error: \%.6f\n", abs(a1_mom$mean - 5)))
cat(sprintf("A2 mean error: \%.2e\n", abs(a2_mom$mean - 5)))

# View iteration trace (includes step size and Jacobian determinant)
head(fit$computation$trace)

}
\references{
Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
\emph{arXiv preprint} arXiv:2602.06301.
}
\seealso{
\code{\link{DPprior_a1}} for closed-form initialization,
\code{\link{moments_with_jacobian}} for Jacobian computation,
\code{\link{exact_K_moments}} for moment verification

Other elicitation:
\code{\link[=DPprior_a1]{DPprior_a1()}},
\code{\link[=DPprior_a2_kl]{DPprior_a2_kl()}},
\code{\link[=DPprior_dual]{DPprior_dual()}},
\code{\link[=DPprior_dual_hard]{DPprior_dual_hard()}},
\code{\link[=DPprior_dual_soft]{DPprior_dual_soft()}},
\code{\link[=DPprior_fit]{DPprior_fit()}},
\code{\link[=DPprior_target_K]{DPprior_target_K()}}
}
\concept{elicitation}

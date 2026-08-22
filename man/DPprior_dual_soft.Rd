% Please edit documentation in R/22_dual_anchor_soft.R
\name{DPprior_dual_soft}
\alias{DPprior_dual_soft}
\title{Calibrate a fixed-scale soft Dual-Anchor trade-off}
\usage{
DPprior_dual_soft(
  fit,
  target,
  lambda,
  max_iter = 100L,
  M_fit = .QUAD_NODES_DEFAULT,
  M_verify = NULL,
  log_bounds = .LOG_BOUNDS_DEFAULT,
  control = list(),
  allow_approximate = FALSE,
  start = NULL,
  ...
)
}
\arguments{
\item{fit}{A K-only \code{DPprior_fit} object.}

\item{target}{Named weight target with \code{metric}, \code{relation}, and
\code{value}; threshold or probability is required when the metric needs
it. Relations are \code{"target"}, \code{"at_most"}, or
\code{"at_least"}.}

\item{lambda}{Mandatory finite scalar satisfying \code{0 < lambda <= 1}.}

\item{max_iter}{Maximum iterations for each optimizer attempt.}

\item{M_fit}{Selected quadrature order.}

\item{M_verify}{Independent verification order; by default at least twice
\code{M_fit} and at least \code{M_fit + 40}.}

\item{log_bounds}{Two finite log-parameter bounds.}

\item{control}{Optimizer and verification controls.}

\item{allow_approximate}{Return a finite but unusable approximate candidate
when \code{TRUE}; otherwise signal a typed condition retaining the result.}

\item{start}{Optional positive \code{c(a,b)} warm start.}

\item{...}{Must be empty. Unknown or legacy hard-mode arguments are rejected
as typed invalid input without evaluating their expressions.}
}
\value{
A canonical \code{dpprior.result/1} \code{DPprior_dual_soft}
object. The \code{tradeoff} extension records the fixed input-derived K
scales, lambda, weight target, component losses, and total objective.
Soft objects never contain \code{constraint_satisfied}. Unapproved
approximate, infeasible, or failed outcomes are signalled as typed
conditions whose \code{condition$result} retains the complete object;
\code{allow_approximate} changes return policy, not scientific status.
}
\description{
Calibrate a fixed-scale soft Dual-Anchor trade-off
}
\details{
A soft result is decision-ready only when its mode-specific
contract and \code{status}, \code{usable}, and \code{verified} fields
permit that use. Lambda is a trade-off weight, not a constraint
probability or a hard-satisfaction certificate.
}
\seealso{
Other elicitation:
\code{\link[=DPprior_a1]{DPprior_a1()}},
\code{\link[=DPprior_a2_kl]{DPprior_a2_kl()}},
\code{\link[=DPprior_a2_newton]{DPprior_a2_newton()}},
\code{\link[=DPprior_dual]{DPprior_dual()}},
\code{\link[=DPprior_dual_hard]{DPprior_dual_hard()}},
\code{\link[=DPprior_fit]{DPprior_fit()}},
\code{\link[=DPprior_target_K]{DPprior_target_K()}}
}
\concept{elicitation}

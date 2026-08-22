% Please edit documentation in R/21_dual_anchor_hard.R
\name{DPprior_dual_hard}
\alias{DPprior_dual_hard}
\title{Dual-anchor calibration with a verified hard inequality}
\usage{
DPprior_dual_hard(
  fit,
  constraint = list(metric = "wsb_tail", threshold = 0.5, relation = "<=", bound = 0.3),
  constraint_tol = list(abs = 1e-06, rel = 1e-06),
  M = .QUAD_NODES_DEFAULT,
  M_verify = NULL,
  log_bounds = .LOG_BOUNDS_DEFAULT,
  control = list(),
  allow_approximate = FALSE,
  ...
)
}
\arguments{
\item{fit}{A \code{DPprior_fit} object retaining positive \code{a}, positive \code{b}, \code{J},
and the requested K mean and variance.}

\item{constraint}{A named list with \code{metric}, \code{relation}, and \code{bound}.
Supported metrics are \code{"wsb_tail"} (requires \code{threshold}),
\code{"wsb_mean"}, \code{"wsb_quantile"} (requires \code{probability}), and
\code{"wmax_tail_upper"} (requires \code{threshold}). Relations may be \code{"<="} or
\code{"at_most"}, and \code{">="} or \code{"at_least"} except that a certified W-max
upper bound supports only an at-most safety constraint.}

\item{constraint_tol}{Named list containing non-negative \code{abs} and \code{rel}
tolerances. Independent satisfaction uses
\code{abs + rel * max(abs(bound), 1e-8)}.}

\item{M}{Fit-order Gauss--Laguerre quadrature order, no greater than 256 so
that a required higher-order verification remains available.}

\item{M_verify}{Independent verification order. It must be at least
\code{max(2*M, M+40)} and no greater than 512. The default is that required
minimum.}

\item{log_bounds}{Two finite increasing bounds shared by \code{log(a)} and
\code{log(b)}.}

\item{control}{Named solver-control list. Supported components are
\code{maxit}, \code{scan_points}, \code{scan_keep}, \code{profile_starts}, \code{root_tol},
\code{optim_reltol}, \code{penalty}, \code{boundary_tol},
\code{verification_abs_tol}, \code{verification_rel_tol}, \code{perturbation_step}, and
\code{perturbation_abs_tol}.}

\item{allow_approximate}{Whether to return an independently verified
feasible scan candidate without successful optimizer evidence. The
result remains visibly \code{status = "approximate"} and intrinsically
\code{usable = FALSE}; this flag changes return-versus-condition behavior only.}

\item{...}{Must be empty. In particular, supplying \code{lambda} raises a typed
invalid-input condition before fitting.}
}
\value{
A canonical \code{dpprior.result/1} \code{DPprior_dual_hard}
object for an approved solution. The \code{constraint} extension contains
the named inequality, achieved metric, signed residual, tolerance, and
independent satisfaction certificate. Certified infeasibility, unknown
feasibility, and an unapproved approximate result raise a typed condition
whose \code{condition$result} contains the complete canonical object.
}
\description{
Minimizes the fixed-input-scale cluster-count loss subject to one named
weight inequality. This hard method has no lambda: lambda belongs only to
soft trade-off calibration.
}
\details{
The signed residual is \code{achieved - bound} for an at-most constraint and
\code{bound - achieved} for an at-least constraint. Thus feasibility always has
the auditable form \code{constraint_residual <= constraint_tolerance$effective}.
\code{constraint_satisfied} is set only from a fresh higher-order calculation;
optimizer exit codes are never treated as constraint evidence.
Decision-ready use requires the returned object's mode-specific contract
together with \code{usable = TRUE} and \code{verified = TRUE}.
}
\seealso{
Other elicitation:
\code{\link[=DPprior_a1]{DPprior_a1()}},
\code{\link[=DPprior_a2_kl]{DPprior_a2_kl()}},
\code{\link[=DPprior_a2_newton]{DPprior_a2_newton()}},
\code{\link[=DPprior_dual]{DPprior_dual()}},
\code{\link[=DPprior_dual_soft]{DPprior_dual_soft()}},
\code{\link[=DPprior_fit]{DPprior_fit()}},
\code{\link[=DPprior_target_K]{DPprior_target_K()}}
}
\concept{elicitation}

% Please edit documentation in R/15_dual_anchor.R
\name{DPprior_dual}
\alias{DPprior_dual}
\title{Legacy Dual-Anchor equality-loss adapter}
\usage{
DPprior_dual(
  fit,
  w1_target,
  lambda = 0.5,
  max_iter = 100L,
  verbose = FALSE,
  M = .QUAD_NODES_DEFAULT,
  loss_type = c("relative", "adaptive", "absolute")
)
}
\arguments{
\item{fit}{A K-only \code{DPprior_fit}.}

\item{w1_target}{Historical W_SB target using \code{prob}, \code{mean}, or
\code{quantile} nesting.}

\item{lambda}{Historical mixing weight in \code{[0,1]}. The new soft API
rejects zero; this compatibility adapter retains it only to reproduce old
analyses.}

\item{max_iter}{Maximum iterations per historical optimizer.}

\item{verbose}{Print historical optimizer progress.}

\item{M}{Quadrature order.}

\item{loss_type}{Historical \code{"relative"}, \code{"adaptive"}, or
\code{"absolute"} scaling.}
}
\value{
A canonical \code{dpprior.result/1} legacy \code{DPprior_fit} with
\code{mode = "dual_legacy"}, \code{status = "approximate"},
\code{usable = TRUE}, and \code{verified = FALSE}. The authoritative
equality-loss, optimizer, scaling, and lifecycle records are
\code{legacy}, \code{provenance}, and \code{computation}. The
\code{dual_anchor} alias is a quarantined, non-authoritative compatibility
view only. Use \code{\link{DPprior_dual_hard}} for a verified inequality
or \code{\link{DPprior_dual_soft}} for a current fixed-scale trade-off.
}
\description{
\code{DPprior_dual()} is retained throughout v2.x, with no removal before
v3.0 and a migration review. It emits one typed lifecycle warning and
reproduces the historical equality-loss objective, including the old
path-derived adaptive scaling. Its result is always labelled legacy,
approximate, and unverified. It never represents a hard inequality or
target-satisfaction certificate.
}
\seealso{
Other elicitation:
\code{\link[=DPprior_a1]{DPprior_a1()}},
\code{\link[=DPprior_a2_kl]{DPprior_a2_kl()}},
\code{\link[=DPprior_a2_newton]{DPprior_a2_newton()}},
\code{\link[=DPprior_dual_hard]{DPprior_dual_hard()}},
\code{\link[=DPprior_dual_soft]{DPprior_dual_soft()}},
\code{\link[=DPprior_fit]{DPprior_fit()}},
\code{\link[=DPprior_target_K]{DPprior_target_K()}}
}
\concept{elicitation}

% Please edit documentation in R/15_dual_anchor.R
\name{dual_anchor_loss}
\alias{dual_anchor_loss}
\title{Legacy Dual-Anchor loss}
\usage{
dual_anchor_loss(
  a,
  b,
  J,
  K_target,
  w1_target,
  lambda,
  M = .QUAD_NODES_DEFAULT,
  loss_type = c("relative", "adaptive", "absolute"),
  L_K_scale = NULL,
  L_w_scale = NULL
)
}
\arguments{
\item{a, b}{Positive Gamma shape and rate.}

\item{J}{Design size.}

\item{K_target}{Legacy list with \code{mu_K} and \code{var_K}.}

\item{w1_target}{Legacy W_SB target.}

\item{lambda}{Legacy mixing weight in \code{[0,1]}.}

\item{M}{Quadrature order.}

\item{loss_type}{Historical scaling rule.}

\item{L_K_scale, L_w_scale}{Historical adaptive scale factors.}
}
\value{
Historical scalar equality loss. This helper is not the v2 soft
contract because adaptive scales are path-derived.
}
\description{
Legacy Dual-Anchor loss
}
\keyword{internal}

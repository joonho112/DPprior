% Please edit documentation in R/15_dual_anchor.R
\name{dual_anchor_diagnostics}
\alias{dual_anchor_diagnostics}
\title{Compare a Dual-Anchor fit with its K-only start}
\usage{
dual_anchor_diagnostics(fit_dual, fit_K_only = NULL, M = .QUAD_NODES_DEFAULT)
}
\arguments{
\item{fit_dual}{A legacy or v2 soft Dual-Anchor fit.}

\item{fit_K_only}{Optional K-only fit.}

\item{M}{Quadrature order.}
}
\value{
A comparison data frame with explicit W_SB labels.
}
\description{
Compare a Dual-Anchor fit with its K-only start
}

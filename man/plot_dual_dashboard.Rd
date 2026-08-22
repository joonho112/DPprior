% Please edit documentation in R/18_visualization.R
\name{plot_dual_dashboard}
\alias{plot_dual_dashboard}
\title{Plot Canonical Dual Comparison Dashboard}
\usage{
plot_dual_dashboard(
  fit_dual,
  tradeoff_data = NULL,
  engine = c("ggplot2", "base"),
  base_size = 10,
  title = NULL,
  show = TRUE
)
}
\arguments{
\item{fit_dual}{A canonical hard or soft dual fit from
\code{DPprior_dual_hard()} or \code{DPprior_dual_soft()}.}

\item{tradeoff_data}{Optional canonical \code{dpprior_tradeoff_curve}. When
supplied, its retained evidence is validated before the fit comparison is
drawn; it is not substituted for the fit's authoritative trade-off record.}

\item{engine}{"ggplot2" (default) or "base".}

\item{base_size}{Base font size.}

\item{title}{Optional title.}

\item{show}{If TRUE, draw the plot.}
}
\value{
A gtable grob, a list of plots when gtable assembly is unavailable,
or invisible(NULL) for the base engine.
}
\description{
Creates the canonical four-panel comparison view: alpha, \eqn{K_J}, and
\eqn{W_{SB}} distributions plus a mode-specific evidence table. Hard fits
show constraint evidence; soft fits show the recorded fixed-scale trade-off.
Retained legacy fits raise a typed comparison-lineage unavailable condition.
}
\references{
Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
\emph{arXiv preprint} arXiv:2602.06301.
}
\seealso{
\code{\link{DPprior_fit}} for fitting, \code{\link{plot.DPprior_fit}} for S3 plot method

Other visualization:
\code{\link[=DPprior_colors]{DPprior_colors()}},
\code{\link[=plot_K_prior]{plot_K_prior()}},
\code{\link[=plot_alpha_prior]{plot_alpha_prior()}},
\code{\link[=plot_dual_comparison]{plot_dual_comparison()}},
\code{\link[=plot_prior_dashboard]{plot_prior_dashboard()}},
\code{\link[=plot_tradeoff_curve]{plot_tradeoff_curve()}},
\code{\link[=plot_tradeoff_dashboard]{plot_tradeoff_dashboard()}},
\code{\link[=plot_w1_prior]{plot_w1_prior()}},
\code{\link[=theme_DPprior]{theme_DPprior()}}
}
\concept{visualization}

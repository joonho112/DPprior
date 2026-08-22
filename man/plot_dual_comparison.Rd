% Please edit documentation in R/18_visualization.R
\name{plot_dual_comparison}
\alias{plot_dual_comparison}
\title{Plot Dual-Anchor Comparison Dashboard}
\usage{
plot_dual_comparison(
  fit_dual,
  fit_K_only = NULL,
  engine = c("ggplot2", "base"),
  base_size = 10,
  title = NULL,
  show = TRUE
)
}
\arguments{
\item{fit_dual}{A canonical hard or soft dual fit from
\code{DPprior_dual_hard()} or \code{DPprior_dual_soft()}. Retained legacy
fits have no authoritative comparison lineage and raise a typed
unavailable condition.}

\item{fit_K_only}{Optional canonical K-only fit. If supplied, it must exactly
match the authoritative \code{provenance.input_fit} reference retained by
\code{fit_dual}.}

\item{engine}{"ggplot2" (default) or "base".}

\item{base_size}{Base font size.}

\item{title}{Optional title for the dashboard.}

\item{show}{If TRUE, draw the plot.}
}
\value{
A gtable grob or list of ggplot objects.
}
\description{
Creates a comparison dashboard showing K-only vs dual-anchor solutions.
Displays changes in alpha, K, and the first size-biased weight
\eqn{W_{SB}} side-by-side.
}
\examples{
fit_K <- DPprior_fit(J = 50, mu_K = 5, var_K = 8)
fit_dual <- DPprior_dual_soft(
  fit_K,
  target = list(
    metric = "wsb_tail", relation = "target",
    threshold = 0.5, value = 0.3
  ),
  lambda = 0.5
)
plot_dual_comparison(fit_dual)

}
\seealso{
\code{\link{DPprior_dual_hard}},
\code{\link{DPprior_dual_soft}}, and \code{\link{plot.DPprior_fit}}

Other visualization:
\code{\link[=DPprior_colors]{DPprior_colors()}},
\code{\link[=plot_K_prior]{plot_K_prior()}},
\code{\link[=plot_alpha_prior]{plot_alpha_prior()}},
\code{\link[=plot_dual_dashboard]{plot_dual_dashboard()}},
\code{\link[=plot_prior_dashboard]{plot_prior_dashboard()}},
\code{\link[=plot_tradeoff_curve]{plot_tradeoff_curve()}},
\code{\link[=plot_tradeoff_dashboard]{plot_tradeoff_dashboard()}},
\code{\link[=plot_w1_prior]{plot_w1_prior()}},
\code{\link[=theme_DPprior]{theme_DPprior()}}
}
\concept{visualization}

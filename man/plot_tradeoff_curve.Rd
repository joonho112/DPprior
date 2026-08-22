% Please edit documentation in R/18_visualization.R
\name{plot_tradeoff_curve}
\alias{plot_tradeoff_curve}
\title{Plot Trade-off Curve}
\usage{
plot_tradeoff_curve(
  tradeoff_data,
  metric = c("w1_prob_gt_50", "E_w1", "K_loss", "mu_K", "var_K"),
  target_value = NULL,
  engine = c("ggplot2", "base"),
  base_size = 11,
  title = NULL,
  show = TRUE
)
}
\arguments{
\item{tradeoff_data}{Data frame from compute_tradeoff_curve().}

\item{metric}{Which metric to plot on y-axis: \code{"w1_prob_gt_50"}
(default), \code{"E_w1"}, \code{"K_loss"}, \code{"mu_K"}, or
\code{"var_K"}. The two legacy column names refer specifically to
\eqn{W_{SB}}, the first size-biased DP weight.}

\item{target_value}{Optional target value to mark with horizontal line.}

\item{engine}{"ggplot2" (default) or "base".}

\item{base_size}{Base font size.}

\item{title}{Optional title.}

\item{show}{If TRUE, print the plot.}
}
\value{
A ggplot object or invisible(NULL).
}
\description{
Visualizes the Pareto trade-off between K_J fit and weight constraint
across different lambda values.
}
\examples{
curve <- compute_tradeoff_curve(
  J = 50,
  K_target = list(mu_K = 5, var_K = 8),
  w1_target = list(prob = list(threshold = 0.5, value = 0.25)),
  lambda_seq = seq(0.1, 1, by = 0.1)
)
plot_tradeoff_curve(curve, target_value = 0.25)

}
\seealso{
\code{\link{DPprior_fit}} for fitting, \code{\link{plot.DPprior_fit}} for S3 plot method

Other visualization:
\code{\link[=DPprior_colors]{DPprior_colors()}},
\code{\link[=plot_K_prior]{plot_K_prior()}},
\code{\link[=plot_alpha_prior]{plot_alpha_prior()}},
\code{\link[=plot_dual_comparison]{plot_dual_comparison()}},
\code{\link[=plot_dual_dashboard]{plot_dual_dashboard()}},
\code{\link[=plot_prior_dashboard]{plot_prior_dashboard()}},
\code{\link[=plot_tradeoff_dashboard]{plot_tradeoff_dashboard()}},
\code{\link[=plot_w1_prior]{plot_w1_prior()}},
\code{\link[=theme_DPprior]{theme_DPprior()}}
}
\concept{visualization}

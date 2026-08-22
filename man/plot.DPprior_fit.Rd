% Please edit documentation in R/17_s3_methods.R
\name{plot.DPprior_fit}
\alias{plot.DPprior_fit}
\title{Plot Method for DPprior_fit Objects}
\usage{
\method{plot}{DPprior_fit}(
  x,
  type = c("auto", "dashboard", "alpha", "K", "w1", "dual", "comparison"),
  engine = c("ggplot2", "base"),
  ...
)
}
\arguments{
\item{x}{A \code{DPprior_fit} object.}

\item{type}{Character; the type of plot to create:
\describe{
\item{"auto"}{(Default) Automatically selects the appropriate plot type.
Uses \code{"dual"} for current hard and soft dual fits. A1, A2, and
retained legacy fits use the descriptive single-fit dashboard.}
\item{"dashboard"}{4-panel dashboard showing alpha, K, the first
size-biased weight \eqn{W_{SB}}, and summary.}
\item{"alpha"}{Prior density of the concentration parameter alpha.}
\item{"K"}{Prior PMF of the number of clusters \eqn{K_J}.}
\item{"w1"}{Prior density of the first size-biased DP weight
\eqn{W_{SB}}.}
\item{"dual"}{Authoritative K-only comparison for current hard or soft
dual fits. Legacy fits raise a typed lineage-unavailable condition.}
\item{"comparison"}{Same as "dual".}
}}

\item{engine}{Character; graphics engine to use:
\code{"ggplot2"} (default) or \code{"base"}.}

\item{...}{Additional arguments passed to the underlying plot functions.
Common options include:
\describe{
\item{base_size}{Base font size (default: 11)}
\item{ci_level}{Credible interval level for alpha plot (default: 0.95)}
\item{title}{Optional title for the dashboard}
\item{show}{If TRUE, display the plot; if FALSE, return silently}
}}
}
\value{
Depends on the plot type and engine:
\itemize{
\item For ggplot2: Returns a ggplot object or gtable (for dashboards)
\item For base: Returns invisible(NULL)
}
}
\description{
Creates visualizations of a prior elicitation result. Multiple plot types
are available, including individual distribution plots and comprehensive
dashboards.
}
\details{
The \code{"auto"} type is recommended for most use cases. It automatically
uses an authoritative comparison only for current hard or soft calibration.
Retained legacy results remain available as status-labelled descriptive
single-fit plots; no compatibility initialization is used as comparison
science.

For dual-anchor fits, the comparison dashboard shows:
\itemize{
\item Alpha prior: K-only vs Dual-anchor
\item K distribution comparison
\item \eqn{W_{SB}} distribution comparison with named tail thresholds
\item Summary comparison table
}
}
\section{Plot Type Details}{

\describe{
\item{dashboard}{
A 2x2 grid showing:
(A) Alpha prior density with CI
(B) \eqn{K_J} prior PMF with mode and mean
(C) first size-biased weight density with threshold shading
(D) Summary statistics table
}
\item{alpha}{
Gamma(a, b) density with:
- Mean line (dashed)
- Credible interval (shaded region)
- Annotation with moments and CI
}
\item{K}{
Bar plot of \eqn{P(K_J = k)} with:
- Target mean line
- Achieved mean line
- Optional CDF overlay
}
\item{w1}{
Density plot with:
- Explicit threshold-region shading for \eqn{W_{SB}}
- Threshold lines
- Estimand-labelled exceedance probabilities
}
}
}

\examples{
# Create a fit object
fit <- DPprior_fit(J = 50, mu_K = 5, var_K = 8)

# Auto-detect best plot type
plot(fit)

# Specific plot types
plot(fit, type = "alpha")
plot(fit, type = "K")
plot(fit, type = "w1")
plot(fit, type = "dashboard")

# With custom options
plot(fit, type = "dashboard", title = "My Prior Analysis")

# Current soft dual comparison
fit_K <- DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8)
fit_dual <- DPprior_dual_soft(
  fit_K,
  target = list(
    metric = "wsb_tail", relation = "target",
    threshold = 0.5, value = 0.3
  ),
  lambda = 0.5
)
plot(fit_dual)  # Auto-selects dual comparison
plot(fit_dual, type = "comparison")  # Explicit

# Retained legacy results use a descriptive single-fit dashboard.
fit_legacy <- suppressWarnings(DPprior_dual(
  fit_K, list(prob = list(threshold = 0.5, value = 0.3))
))
plot(fit_legacy)
try(plot(fit_legacy, type = "comparison")) # typed lineage unavailable

}
\seealso{
\code{\link{plot_prior_dashboard}}, \code{\link{plot_alpha_prior}},
\code{\link{plot_K_prior}}, \code{\link{plot_w1_prior}},
\code{\link{plot_dual_comparison}}
}

% Please edit documentation in R/08_weights_w1.R
\name{summary_w1}
\alias{summary_w1}
\title{Summary Statistics for the First Size-Biased DP Weight}
\usage{
summary_w1(
  a,
  b,
  probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
  M = .QUAD_NODES_DEFAULT
)
}
\arguments{
\item{a}{Numeric; shape parameter of the Gamma prior on \eqn{\alpha}
(a > 0).}

\item{b}{Numeric; rate parameter of the Gamma prior on \eqn{\alpha}
(b > 0).}

\item{probs}{Numeric vector; quantile probabilities. Default is
\code{c(0.05, 0.25, 0.5, 0.75, 0.95)}.}

\item{M}{Integer; number of quadrature nodes for mean/variance. Default is 80.}
}
\value{
A list of class "w1_summary" containing:
\describe{
\item{mean}{\eqn{E(w_1)}}
\item{var}{\eqn{Var(w_1)}}
\item{sd}{\eqn{SD(w_1) = \sqrt{Var(w_1)}}}
\item{median}{Median of \eqn{w_1}}
\item{quantiles}{Named vector of quantiles}
\item{prob_gt_50}{P(W_SB > 0.5), a size-biased weight tail}
\item{prob_gt_90}{P(W_SB > 0.9), a size-biased weight tail}
\item{params}{List of input parameters (a, b)}
\item{estimand,label,conditioning,provenance}{Explicit scientific and
computational metadata}
}
}
\description{
Computes comprehensive summary statistics for the \eqn{w_1} distribution.
}
\examples{
# Standard summary
summary_w1(a = 2, b = 1)

# Lee et al. DP-inform prior
summary_w1(a = 1.6, b = 1.22)

}
\seealso{
\code{\link{cdf_w1}}, \code{\link{quantile_w1}}, \code{\link{mean_w1}}

Other weights_w1:
\code{\link[=cdf_w1]{cdf_w1()}},
\code{\link[=density_w1]{density_w1()}},
\code{\link[=mean_w1]{mean_w1()}},
\code{\link[=prob_w1_exceeds]{prob_w1_exceeds()}},
\code{\link[=prob_wmax_exceeds]{prob_wmax_exceeds()}},
\code{\link[=prob_wsb_exceeds]{prob_wsb_exceeds()}},
\code{\link[=quantile_w1]{quantile_w1()}},
\code{\link[=var_w1]{var_w1()}},
\code{\link[=wmax_tail_bounds]{wmax_tail_bounds()}}
}
\concept{weights_w1}

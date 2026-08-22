% Please edit documentation in R/08_weights_w1.R
\name{prob_wsb_exceeds}
\alias{prob_wsb_exceeds}
\title{Tail Probability of the First Size-Biased DP Weight}
\usage{
prob_wsb_exceeds(threshold, a, b)
}
\arguments{
\item{threshold}{Numeric vector of thresholds.}

\item{a, b}{Positive scalar shape and rate of the Gamma prior on alpha.}
}
\value{
Numeric vector containing \eqn{P(W_{SB}>threshold)}.
}
\description{
Explicitly named alias for \code{\link{prob_w1_exceeds}}. The legacy
\code{w1} name remains mathematically valid because the first GEM weight is
a size-biased pick, but \code{W_SB} makes the estimand unambiguous.
}
\examples{
prob_wsb_exceeds(0.5, a = 2, b = 1)

}
\seealso{
\code{\link{prob_w1_exceeds}},
\code{\link{prob_wmax_exceeds}}

Other weights_w1:
\code{\link[=cdf_w1]{cdf_w1()}},
\code{\link[=density_w1]{density_w1()}},
\code{\link[=mean_w1]{mean_w1()}},
\code{\link[=prob_w1_exceeds]{prob_w1_exceeds()}},
\code{\link[=prob_wmax_exceeds]{prob_wmax_exceeds()}},
\code{\link[=quantile_w1]{quantile_w1()}},
\code{\link[=summary_w1]{summary_w1()}},
\code{\link[=var_w1]{var_w1()}},
\code{\link[=wmax_tail_bounds]{wmax_tail_bounds()}}
}
\concept{weights_w1}

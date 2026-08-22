% Please edit documentation in R/08_weights_w1.R
\name{wmax_tail_bounds}
\alias{wmax_tail_bounds}
\title{Certified Bounds for the Largest DP Population Weight Tail}
\usage{
wmax_tail_bounds(threshold, alpha = NULL, a = NULL, b = NULL)
}
\arguments{
\item{threshold}{One finite probability in the closed unit interval.}

\item{alpha}{Optional positive scalar fixed DP concentration.}

\item{a, b}{Optional positive scalar Gamma shape and rate. Supply either
\code{alpha} or both \code{a} and \code{b}, never both specifications.}
}
\value{
A \code{wmax_tail_bounds} object. The fields
\code{lower_bound} and \code{upper_bound} are certified bounds, not point
estimates. \code{size_biased_tail} is \eqn{P(W_{SB}>threshold)}.
}
\description{
Computes distribution-free bounds on
\eqn{P(W_{max}>threshold)}. If \eqn{q_t=P(W_{SB}>threshold)}, then
\deqn{q_t \le P(W_{max}>threshold) \le \min(q_t/threshold,1).}
The identity remains valid after mixing alpha over a Gamma hyperprior.
}
\examples{
wmax_tail_bounds(0.4, alpha = 1)
wmax_tail_bounds(0.5, a = 2, b = 1)

}
\seealso{
\code{\link{prob_wmax_exceeds}},
\code{\link{prob_wsb_exceeds}}

Other weights_w1:
\code{\link[=cdf_w1]{cdf_w1()}},
\code{\link[=density_w1]{density_w1()}},
\code{\link[=mean_w1]{mean_w1()}},
\code{\link[=prob_w1_exceeds]{prob_w1_exceeds()}},
\code{\link[=prob_wmax_exceeds]{prob_wmax_exceeds()}},
\code{\link[=prob_wsb_exceeds]{prob_wsb_exceeds()}},
\code{\link[=quantile_w1]{quantile_w1()}},
\code{\link[=summary_w1]{summary_w1()}},
\code{\link[=var_w1]{var_w1()}}
}
\concept{weights_w1}

% Please edit documentation in R/08_weights_w1.R
\name{prob_w1_exceeds}
\alias{prob_w1_exceeds}
\title{Tail Probability of the First Size-Biased DP Weight}
\usage{
prob_w1_exceeds(t, a, b)
}
\arguments{
\item{t}{Numeric vector of thresholds. Values outside the unit interval are
allowed but are mapped to the boundary values (1 for \eqn{t \le 0}, 0 for
\eqn{t \ge 1}).}

\item{a}{Numeric; shape parameter of the Gamma prior on \eqn{\alpha}
(a > 0).}

\item{b}{Numeric; rate parameter of the Gamma prior on \eqn{\alpha}
(b > 0).}
}
\value{
Numeric vector of survival probabilities.
}
\description{
Computes P(W_SB > t | a, b) = 1 - F(t), the probability that a randomly
selected unit belongs to a population cluster whose DP mass exceeds t.
}
\details{
The survival function has the closed form:
\deqn{P(w_1 > t | a, b) = \left(\frac{b}{b - \log(1-t)}\right)^a}

This is a mass-weighted cluster-size probability. It is not
\eqn{P(W_{max}>t)}, the probability that the largest population weight
exceeds \eqn{t}. Use \code{\link{prob_wmax_exceeds}} for that estimand.
}
\examples{
# Probability that a random unit's population cluster has mass above 0.5
prob_w1_exceeds(0.5, a = 1.6, b = 1.22)  # ~0.487 (Lee et al. DP-inform)
prob_w1_exceeds(0.5, a = 2, b = 1)       # ~0.349

}
\references{
Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
\emph{arXiv preprint} arXiv:2602.06301.
}
\seealso{
\code{\link{cdf_w1}}, \code{\link{summary_w1}}

Other weights_w1:
\code{\link[=cdf_w1]{cdf_w1()}},
\code{\link[=density_w1]{density_w1()}},
\code{\link[=mean_w1]{mean_w1()}},
\code{\link[=prob_wmax_exceeds]{prob_wmax_exceeds()}},
\code{\link[=prob_wsb_exceeds]{prob_wsb_exceeds()}},
\code{\link[=quantile_w1]{quantile_w1()}},
\code{\link[=summary_w1]{summary_w1()}},
\code{\link[=var_w1]{var_w1()}},
\code{\link[=wmax_tail_bounds]{wmax_tail_bounds()}}
}
\concept{weights_w1}

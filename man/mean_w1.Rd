% Please edit documentation in R/08_weights_w1.R
\name{mean_w1}
\alias{mean_w1}
\title{Mean of the First Size-Biased DP Weight}
\usage{
mean_w1(a, b, M = .QUAD_NODES_DEFAULT)
}
\arguments{
\item{a}{Numeric; shape parameter of the Gamma prior on \eqn{\alpha}
(a > 0).}

\item{b}{Numeric; rate parameter of the Gamma prior on \eqn{\alpha}
(b > 0).}

\item{M}{Integer; number of quadrature nodes. Default is 80.}
}
\value{
Numeric; \eqn{E(w_1)}.
}
\description{
Computes \eqn{E(w_1 \mid a, b)} via Gauss-Laguerre quadrature.
}
\details{
The expectation is computed using the identity:
\deqn{E[w_1 | a, b] = E\left[\frac{1}{1+\alpha}\right] = I_1(a, b)}

where the integral is evaluated via Gauss-Laguerre quadrature.

\strong{Key identity:} \eqn{E(w_1 \mid a, b) = E(\rho \mid a, b)}, where
\eqn{\rho = \sum_h w_h^2} is the co-clustering probability.
}
\examples{
mean_w1(a = 2, b = 1)       # ~0.404
mean_w1(a = 1.6, b = 1.22)  # ~0.508

}
\references{
Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
\emph{arXiv preprint} arXiv:2602.06301.
}
\seealso{
\code{\link{var_w1}}, \code{\link{summary_w1}}

Other weights_w1:
\code{\link[=cdf_w1]{cdf_w1()}},
\code{\link[=density_w1]{density_w1()}},
\code{\link[=prob_w1_exceeds]{prob_w1_exceeds()}},
\code{\link[=prob_wmax_exceeds]{prob_wmax_exceeds()}},
\code{\link[=prob_wsb_exceeds]{prob_wsb_exceeds()}},
\code{\link[=quantile_w1]{quantile_w1()}},
\code{\link[=summary_w1]{summary_w1()}},
\code{\link[=var_w1]{var_w1()}},
\code{\link[=wmax_tail_bounds]{wmax_tail_bounds()}}
}
\concept{weights_w1}

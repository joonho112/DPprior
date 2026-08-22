% Please edit documentation in R/08_weights_w1.R
\name{var_w1}
\alias{var_w1}
\title{Variance of the First Size-Biased DP Weight}
\usage{
var_w1(a, b, M = .QUAD_NODES_DEFAULT)
}
\arguments{
\item{a}{Numeric; shape parameter of the Gamma prior on \eqn{\alpha}
(a > 0).}

\item{b}{Numeric; rate parameter of the Gamma prior on \eqn{\alpha}
(b > 0).}

\item{M}{Integer; number of quadrature nodes. Default is 80.}
}
\value{
Numeric; \eqn{Var(w_1)}.
}
\description{
Computes \eqn{Var(w_1 \mid a, b)} using the law of total variance.
}
\details{
Uses the law of total variance:
\deqn{Var(w_1) = E[Var(w_1 | \alpha)] + Var(E[w_1 | \alpha])}

where \eqn{w_1 \mid \alpha \sim Beta(1, \alpha)}, so:
\itemize{
\item \eqn{E(w_1 \mid \alpha) = 1/(1+\alpha)}
\item \eqn{Var(w_1 \mid \alpha) = \alpha / ((1+\alpha)^2(2+\alpha))}
}
}
\examples{
var_w1(a = 2, b = 1)       # ~0.090
var_w1(a = 1.6, b = 1.22)  # ~0.105

}
\references{
Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
\emph{arXiv preprint} arXiv:2602.06301.
}
\seealso{
\code{\link{mean_w1}}, \code{\link{summary_w1}}

Other weights_w1:
\code{\link[=cdf_w1]{cdf_w1()}},
\code{\link[=density_w1]{density_w1()}},
\code{\link[=mean_w1]{mean_w1()}},
\code{\link[=prob_w1_exceeds]{prob_w1_exceeds()}},
\code{\link[=prob_wmax_exceeds]{prob_wmax_exceeds()}},
\code{\link[=prob_wsb_exceeds]{prob_wsb_exceeds()}},
\code{\link[=quantile_w1]{quantile_w1()}},
\code{\link[=summary_w1]{summary_w1()}},
\code{\link[=wmax_tail_bounds]{wmax_tail_bounds()}}
}
\concept{weights_w1}

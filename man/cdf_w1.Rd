% Please edit documentation in R/08_weights_w1.R
\name{cdf_w1}
\alias{cdf_w1}
\title{CDF of the First Size-Biased DP Weight}
\usage{
cdf_w1(x, a, b)
}
\arguments{
\item{x}{Numeric vector. Values outside the unit interval are allowed but
are mapped to the boundary values of the CDF (0 for \eqn{x \le 0}, 1 for
\eqn{x \ge 1}).}

\item{a}{Numeric; shape parameter of the Gamma prior on \eqn{\alpha}
(a > 0).}

\item{b}{Numeric; rate parameter of the Gamma prior on \eqn{\alpha}
(b > 0).}
}
\value{
Numeric vector of CDF values F(x | a, b) with same length as x.
}
\description{
Computes \eqn{P(w_1 \le x \mid a, b)} using the closed-form expression
derived by marginalizing over \eqn{\alpha \sim Gamma(a, b)}.
}
\details{
The unconditional CDF is given by:
\deqn{F_{w_1}(x | a, b) = 1 - \left(\frac{b}{b - \log(1-x)}\right)^a}

The implementation uses \code{log1p} and \code{expm1} for numerical
stability, particularly when the CDF is close to 0 (small x).
}
\section{Interpretation}{

The weight \eqn{w_1} is in \strong{GEM (size-biased) order}, not ranked by size.
It represents the asymptotic cluster share of a randomly chosen unit,
\strong{not} the largest cluster proportion. See Lee (2026, Section 4) for details.
}

\examples{
# P(w1 <= 0.3) under standard prior
cdf_w1(0.3, a = 2, b = 1)

# Vectorized computation
cdf_w1(c(0.1, 0.3, 0.5, 0.7), a = 1.6, b = 1.22)

}
\references{
Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
\emph{arXiv preprint} arXiv:2602.06301.

Vicentini, C. and Jermyn, I. H. (2025). Prior selection for the precision
parameter of Dirichlet Process Mixtures. arXiv:2502.00864.
}
\seealso{
\code{\link{quantile_w1}}, \code{\link{density_w1}},
\code{\link{prob_w1_exceeds}}

Other weights_w1:
\code{\link[=density_w1]{density_w1()}},
\code{\link[=mean_w1]{mean_w1()}},
\code{\link[=prob_w1_exceeds]{prob_w1_exceeds()}},
\code{\link[=prob_wmax_exceeds]{prob_wmax_exceeds()}},
\code{\link[=prob_wsb_exceeds]{prob_wsb_exceeds()}},
\code{\link[=quantile_w1]{quantile_w1()}},
\code{\link[=summary_w1]{summary_w1()}},
\code{\link[=var_w1]{var_w1()}},
\code{\link[=wmax_tail_bounds]{wmax_tail_bounds()}}
}
\concept{weights_w1}

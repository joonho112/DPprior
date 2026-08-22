% Please edit documentation in R/08_weights_w1.R
\name{density_w1}
\alias{density_w1}
\title{Density of the First Size-Biased DP Weight}
\usage{
density_w1(x, a, b, log = FALSE)
}
\arguments{
\item{x}{Numeric vector; evaluation points.}

\item{a}{Numeric; shape parameter of the Gamma prior on \eqn{\alpha}
(a > 0).}

\item{b}{Numeric; rate parameter of the Gamma prior on \eqn{\alpha}
(b > 0).}

\item{log}{Logical; if \code{TRUE}, returns log-density. Default is \code{FALSE}.}
}
\value{
Numeric vector of density (or log-density if \code{log = TRUE}) values.
Returns 0 (or -Inf on log scale) for x outside (0, 1).
}
\description{
Computes the probability density \eqn{p(w_1 = x \mid a, b)}.
}
\details{
The marginal density of \eqn{w_1} is:
\deqn{p(w_1 | a, b) = \frac{a \cdot b^a}{(1-w_1) \cdot [b - \log(1-w_1)]^{a+1}}}

\strong{Important:} For small values of a (a < 1), the density has significant
mass concentrated very close to x = 1.
}
\examples{
# Density at several points
x <- seq(0.1, 0.9, by = 0.1)
density_w1(x, a = 2, b = 1)

# Log-density for numerical stability
density_w1(0.5, a = 2, b = 1, log = TRUE)

}
\seealso{
\code{\link{cdf_w1}}, \code{\link{quantile_w1}}

Other weights_w1:
\code{\link[=cdf_w1]{cdf_w1()}},
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

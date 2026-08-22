% Please edit documentation in R/08_weights_w1.R
\name{rw1}
\alias{rw1}
\title{Random Generation from the First Size-Biased DP Weight Distribution}
\usage{
rw1(n, a, b)
}
\arguments{
\item{n}{Integer; number of samples to generate.}

\item{a}{Numeric; shape parameter of the Gamma prior on \eqn{\alpha}
(a > 0).}

\item{b}{Numeric; rate parameter of the Gamma prior on \eqn{\alpha}
(b > 0).}
}
\value{
Numeric vector of length n; random samples from the \eqn{w_1}
distribution.
}
\description{
Generates random samples from the \eqn{w_1} distribution by first sampling
\eqn{\alpha \sim Gamma(a, b)}, then
\eqn{w_1 \mid \alpha \sim Beta(1, \alpha)}.
}
\details{
This function uses the hierarchical representation:
\enumerate{
\item \eqn{\alpha \sim Gamma(a, b)}
\item \eqn{w_1 \mid \alpha \sim Beta(1, \alpha)}
}

Useful for Monte Carlo validation of the closed-form functions.
}
\examples{
# Generate samples
set.seed(42)
samples <- rw1(10000, a = 2, b = 1)

# Compare empirical vs theoretical mean
mean(samples)          # ~0.404
mean_w1(a = 2, b = 1)  # 0.4037

}

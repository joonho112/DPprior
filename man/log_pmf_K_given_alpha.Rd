% Please edit documentation in R/04_pmf_conditional.R
\name{log_pmf_K_given_alpha}
\alias{log_pmf_K_given_alpha}
\title{Log-PMF of K Given Alpha (Antoniak Distribution)}
\usage{
log_pmf_K_given_alpha(J, alpha, logS)
}
\arguments{
\item{J}{Integer; sample size (number of observations, must be >= 1).}

\item{alpha}{Numeric; DP concentration parameter (must be positive scalar).}

\item{logS}{Matrix; pre-computed log-Stirling matrix from
\code{\link{compute_log_stirling}}.}
}
\value{
Numeric vector of length \eqn{J+1} containing
\eqn{\log P(K_J = k \mid \alpha)} for \eqn{k = 0, 1, \ldots, J}.
Entry \code{[1]} corresponds to \eqn{k=0} and always equals \code{-Inf}.
The mathematical support is \eqn{\{1,\ldots,J\}}; the leading structural
zero is retained in the vector representation for compatibility.
}
\description{
Computes \eqn{\log P(K_J = k \mid \alpha)} for \eqn{k = 0, 1, \ldots, J}.
}
\details{
Uses the Antoniak distribution formula in log-space:
\deqn{\log P(K_J = k \mid \alpha) = \log|s(J,k)| + k\log\alpha - \log(\alpha)_J}

where \eqn{|s(J,k)|} is the unsigned Stirling number of the first kind and
\eqn{(\alpha)_J} is the rising factorial. The implementation normalizes the
equivalent log-weights over \eqn{k=1,\ldots,J} with log-sum-exp. This avoids
cancellation in a separately evaluated log-gamma normalizer at very large
finite \eqn{\alpha}.

This log-space computation is numerically stable for large \eqn{J} where
direct computation would overflow.
}
\examples{
logS <- compute_log_stirling(50)
log_pmf <- log_pmf_K_given_alpha(50, 2.0, logS)

# Convert to probabilities (numerically stable softmax)
pmf <- exp(log_pmf - max(log_pmf))
pmf <- pmf / sum(pmf)
sum(pmf)  # Should be 1

# Or use pmf_K_given_alpha() directly
pmf2 <- pmf_K_given_alpha(50, 2.0, logS)
sum(pmf2)  # Should be 1

}
\seealso{
\code{\link{pmf_K_given_alpha}} for normalized PMF,
\code{\link{compute_log_stirling}} for Stirling computation

Other conditional_K:
\code{\link[=cdf_K_given_alpha]{cdf_K_given_alpha()}},
\code{\link[=cv_K_given_alpha]{cv_K_given_alpha()}},
\code{\link[=dispersion_K_given_alpha]{dispersion_K_given_alpha()}},
\code{\link[=mean_K_given_alpha]{mean_K_given_alpha()}},
\code{\link[=mode_K_given_alpha]{mode_K_given_alpha()}},
\code{\link[=moments_K_given_alpha]{moments_K_given_alpha()}},
\code{\link[=pmf_K_given_alpha]{pmf_K_given_alpha()}},
\code{\link[=quantile_K_given_alpha]{quantile_K_given_alpha()}},
\code{\link[=var_K_given_alpha]{var_K_given_alpha()}}
}
\concept{conditional_K}

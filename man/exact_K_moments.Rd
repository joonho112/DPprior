% Please edit documentation in R/05_moments_marginal.R
\name{exact_K_moments}
\alias{exact_K_moments}
\title{Marginal Moments of K_J under Gamma Prior}
\usage{
exact_K_moments(
  J,
  a,
  b,
  M = .QUAD_NODES_DEFAULT,
  M_verify = NULL,
  abs_tol = 1e-10,
  rel_tol = 1e-08,
  strict = FALSE
)
}
\arguments{
\item{J}{Integer; sample size (positive integer >= 1).}

\item{a}{Numeric; shape parameter of Gamma prior (> 0).}

\item{b}{Numeric; rate parameter of Gamma prior (> 0).}

\item{M}{Integer; number of quadrature nodes (default: 80).}

\item{M_verify}{Optional integer quadrature order at least
\code{max(2*M, M+40)} and no greater than 512. When supplied, the
selected-order moments are compared with a fresh higher-order computation.
When omitted, the result is explicitly marked as an unverified
approximation. For \code{M > 256}, no admissible verification order exists;
only omitted \code{M_verify} with \code{strict = FALSE} is allowed.}

\item{abs_tol, rel_tol}{Non-negative absolute and relative tolerances used
for selected-versus-verification comparisons.}

\item{strict}{Logical; if \code{TRUE}, require higher-order verification and
raise a typed convergence error when it is absent or fails.}
}
\value{
A named list with components:
\describe{
\item{\code{mean}}{Marginal mean \eqn{E[K_J | a, b]}}
\item{\code{var}}{Marginal variance \eqn{Var(K_J | a, b)}}
\item{\code{sd}}{Marginal standard deviation}
\item{\code{cv}}{Coefficient of variation (sd/mean)}
\item{\code{status}}{Either \code{"converged"} after successful
higher-order verification or \code{"approximate"}.}
\item{\code{decomposition}}{Within-alpha and between-alpha components
of the law of total variance.}
\item{\code{quadrature}}{Selected/verification order, tolerances,
discrepancies, and node-rule provenance.}
}
}
\description{
Computes the marginal mean \eqn{E[K_J]} and variance \eqn{Var(K_J)}
when the DP concentration parameter follows a Gamma(a, b) prior.
}
\details{
Uses Gauss-Laguerre quadrature to numerically evaluate:
\deqn{M_1(a,b) = E_{\alpha \sim \Gamma(a,b)}[\mu_J(\alpha)]}
\deqn{V(a,b) = E[v_J(\alpha)] + E[\mu_J(\alpha)^2] - M_1^2}

where \eqn{\mu_J(\alpha)} and \eqn{v_J(\alpha)} are the conditional
mean and variance from Module 03.

The Law of Total Variance decomposes the marginal variance into:
\itemize{
\item Within-alpha variance: \eqn{E[v_J(\alpha)]}
\item Between-alpha variance: \eqn{Var(\mu_J(\alpha))}
}

The function name is retained for compatibility. The returned values use a
fixed Gauss--Laguerre rule and are therefore numerical approximations.
Supply \code{M_verify} to obtain an explicit higher-order agreement check;
the selected-order result is always returned and is never silently replaced.
The quadrature metadata reports the required verification order and whether
that order is available under the 512-node implementation ceiling.

\strong{Key properties:}
\itemize{
\item The mean is bounded: \eqn{1 \leq E[K_J] \leq J}
\item Gamma mixing adds the non-negative between-\eqn{\alpha} component
to the average conditional variance; this identity alone does not imply
a universal marginal variance-to-mean ratio
\item The marginal variance equals the average conditional variance plus
the non-negative between-\eqn{\alpha} variance of the conditional mean
}
}
\examples{
# Example: J=50, Gamma(1.5, 0.5) prior
result <- exact_K_moments(50, 1.5, 0.5)
print(result)

# Verify the law-of-total-variance decomposition
sum(unlist(result$decomposition)) - result$var

# Verify mean bounds
1 <= result$mean && result$mean <= 50  # TRUE

}
\references{
Antoniak, C. E. (1974). Mixtures of Dirichlet Processes.
\emph{The Annals of Statistics}, 2(6), 1152-1174.
}
\seealso{
\code{\link{K_moments}} for convenience wrapper,
\code{\link{mean_K_given_alpha}}, \code{\link{var_K_given_alpha}}

Other marginal_K:
\code{\link[=K_moments]{K_moments()}},
\code{\link[=cdf_K_marginal]{cdf_K_marginal()}},
\code{\link[=mode_K_marginal]{mode_K_marginal()}},
\code{\link[=pmf_K_marginal]{pmf_K_marginal()}},
\code{\link[=quantile_K_marginal]{quantile_K_marginal()}},
\code{\link[=summary_K_marginal]{summary_K_marginal()}}
}
\concept{marginal_K}

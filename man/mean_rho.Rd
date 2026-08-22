% Please edit documentation in R/09_weights_rho.R
\name{mean_rho}
\alias{mean_rho}
\title{Marginal Mean of rho}
\usage{
mean_rho(a, b, M = .QUAD_NODES_DEFAULT)
}
\arguments{
\item{a}{Numeric; shape parameter of the Gamma prior on alpha (a > 0).}

\item{b}{Numeric; rate parameter of the Gamma prior on alpha (b > 0).}

\item{M}{Integer; number of quadrature nodes. Default is 80.}
}
\value{
Numeric; \code{E(rho | a, b)}.
}
\description{
Computes E(rho | a, b) when alpha ~ Gamma(a, b) (shape-rate).
}
\details{
Uses Gauss-Laguerre quadrature via \code{integrate_gamma}. A key identity is
\code{E(rho | alpha) = E(w1 | alpha) = 1/(1+alpha)}, so \code{E(rho | a, b)} equals
\code{E(w1 | a, b)} (but the full distributions differ).

The result is the marginal probability that two exchangeable units share a
cluster, averaged over both the random DP weights and the Gamma hyperprior.
No qualitative category or action threshold is imposed by this function.
}
\examples{
mean_rho(a = 2, b = 1)
mean_rho(a = 1.6, b = 1.22)

}
\references{
Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
\emph{arXiv preprint} arXiv:2602.06301.
}
\seealso{
\code{\link{var_rho}}, \code{\link{cv_rho}}, \code{\link{mean_w1}}

Other co_clustering:
\code{\link[=cv_rho]{cv_rho()}},
\code{\link[=mean_rho_given_alpha]{mean_rho_given_alpha()}},
\code{\link[=var_rho]{var_rho()}},
\code{\link[=var_rho_given_alpha]{var_rho_given_alpha()}}
}
\concept{co_clustering}

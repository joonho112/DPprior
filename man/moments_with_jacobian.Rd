% Please edit documentation in R/07_jacobian.R
\name{moments_with_jacobian}
\alias{moments_with_jacobian}
\title{Compute Marginal Moments and Jacobian Simultaneously}
\usage{
moments_with_jacobian(J, a, b, M = .QUAD_NODES_DEFAULT)
}
\arguments{
\item{J}{Integer; sample size (number of observations/sites).}

\item{a}{Numeric; shape parameter of the Gamma prior on \eqn{\alpha} (> 0).}

\item{b}{Numeric; rate parameter of the Gamma prior on \eqn{\alpha} (> 0).}

\item{M}{Integer; number of quadrature nodes (default: 80).}
}
\value{
A named list with components:
\describe{
\item{\code{mean}}{Marginal mean \eqn{E[K_J]}}
\item{\code{var}}{Marginal variance \eqn{Var(K_J)}}
\item{\code{jacobian}}{2x2 Jacobian matrix with structure:
\deqn{J_F = \begin{bmatrix} \partial M_1/\partial a & \partial M_1/\partial b \\
                                 \partial V/\partial a & \partial V/\partial b \end{bmatrix}}}
\item{\code{derivative_diagnostics}}{Order-refinement status, tolerances,
and quadrature orders used for the score derivatives.}
\item{\code{conditioning}}{Scale-aware singularity diagnostics for both
shape-rate and solver log-parameterizations.}
}
}
\description{
Computes quadrature approximations to the marginal moments
\eqn{M_1 = E[K_J]} and \eqn{V = Var(K_J)}
along with the Jacobian matrix of the moment map \eqn{F(a,b) = (M_1, V)}
using score function identities.
}
\details{
This function uses the score identity (Lee, 2026, Section 3.2, Corollary 1)
to compute derivatives without finite differences:
\deqn{\frac{\partial}{\partial\theta} E[f(\alpha)] = E[f(\alpha) \cdot s_\theta(\alpha)]}

The Jacobian components are computed as:
\itemize{
\item \eqn{\partial M_1/\partial \theta =
    E[(\mu_J(\alpha)-1) \cdot s_\theta(\alpha)]}
\item \eqn{\partial V/\partial \theta = \partial E[v_J]/\partial \theta +
             \partial E[\mu_J^2]/\partial \theta - 2 M_1 \partial M_1/\partial \theta}
}
Subtracting the exact zero-concentration limits is an exact control variate
because \eqn{E[s_\theta(\alpha)]=0}. It removes the spurious constant times
\code{log(alpha)} term that otherwise dominates finite-order shape-score
quadrature. Derivatives are recomputed at a distinct higher order when the
supported ceiling permits; the returned diagnostic is \code{converged} only
when a full scheduled refinement is available and the two orders meet the
release mixed error budget
\code{5e-5 + 2e-3 * max(abs(selected), abs(verification))} componentwise.

\strong{Numerical Considerations:}
\itemize{
\item The score function \code{s_a} contains \code{log(alpha)}, which causes
slower quadrature convergence compared to moment computation.
\item The score calculation uses at least 80 nodes and a distinct higher
order when available; both orders are retained in diagnostics.
\item Conditional moments use stable finite sums without a small-alpha
threshold or negative clipping.
\item The returned marginal variance uses the non-negative within/between
decomposition rather than cancellation-prone second-moment
subtraction.
\item Conditioning status is based on the solver-relevant log-parameter
Jacobian, not on an unscaled determinant alone.
}
}
\examples{
# Compute moments and Jacobian for J=50, a=2, b=1
result <- moments_with_jacobian(J = 50, a = 2.0, b = 1.0)
print(result$mean)      # E[K_J]
print(result$var)       # Var(K_J)
print(result$jacobian)  # 2x2 Jacobian matrix

# Use in Newton iteration
target <- c(5.0, 8.0)  # Target (E[K], Var(K))
current <- c(result$mean, result$var)
residual <- current - target
delta <- solve(result$jacobian, -residual)

}
\references{
Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
\emph{arXiv preprint} arXiv:2602.06301.
}
\seealso{
\code{\link{exact_K_moments}} for moments only,
\code{\link{score_a}}, \code{\link{score_b}} for score functions
}

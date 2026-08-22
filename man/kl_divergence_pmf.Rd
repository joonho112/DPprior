% Please edit documentation in R/12_a2_kl.R
\name{kl_divergence_pmf}
\alias{kl_divergence_pmf}
\title{KL Divergence Between Two PMFs}
\usage{
kl_divergence_pmf(p, q, eps = NULL)
}
\arguments{
\item{p}{Numeric vector; target PMF (reference distribution).}

\item{q}{Numeric vector; comparison PMF.}

\item{eps}{Deprecated compatibility argument. Exact mathematical KL does
not smooth zero probabilities; supplying a non-\code{NULL} value is an
error.}
}
\value{
Numeric scalar; the KL divergence (non-negative).
}
\description{
Computes the Kullback-Leibler divergence \eqn{D_{KL}(p \| q)} between two
probability mass functions.
}
\details{
The KL divergence is defined as:
\deqn{D_{KL}(p \| q) = \sum_k p(k) \log\frac{p(k)}{q(k)}}

Terms with \eqn{p(k)=0} contribute zero. If \eqn{p(k)>0} and
\eqn{q(k)=0}, the result is \code{Inf}. Both inputs must already be valid,
normalized PMFs; this function never normalizes or smooths them. In keeping
with the package PMF contract, a sum differing from one by at most
\code{.TOL_PMF_SUM} is accepted without modifying the vector. A raw negative
value whose magnitude is fully explained by that validation tolerance and
floating-point roundoff is reported as zero.

\strong{Properties:}
\itemize{
\item \eqn{D_{KL}(p \| q) \geq 0}; for exactly normalized PMFs, equality
holds iff \eqn{p = q}. For sums accepted within the PMF tolerance,
equality is interpreted up to that same numerical tolerance.
\item Not symmetric: \eqn{D_{KL}(p \| q) \neq D_{KL}(q \| p)}
}
}
\examples{
p <- c(0.2, 0.5, 0.3)
kl_divergence_pmf(p, p)  # 0

q <- c(0.3, 0.4, 0.3)
kl_divergence_pmf(p, q)

}
\seealso{
\code{\link{kl_divergence_K}} for KL divergence with induced PMF
}

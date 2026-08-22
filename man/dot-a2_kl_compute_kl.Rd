% Please edit documentation in R/12_a2_kl.R
\name{.a2_kl_compute_kl}
\alias{.a2_kl_compute_kl}
\title{Compute KL Divergence from Target to Induced PMF}
\usage{
.a2_kl_compute_kl(target_pmf, induced_pmf)
}
\arguments{
\item{target_pmf}{Numeric vector; normalized PMF (length J).}

\item{induced_pmf}{Numeric vector; normalized induced PMF (length J).}
}
\value{
Numeric; KL divergence (non-negative).
}
\description{
Compute KL Divergence from Target to Induced PMF
}
\keyword{internal}

% Please edit documentation in R/12_a2_kl.R
\name{.a2_kl_normalize_pmf}
\alias{.a2_kl_normalize_pmf}
\title{Validate a Target PMF on \{1, ..., J\}}
\usage{
.a2_kl_normalize_pmf(target_pmf, J)
}
\arguments{
\item{target_pmf}{Numeric vector; target PMF.}

\item{J}{Integer; sample size.}
}
\value{
Numeric vector of length J (support k=1,...,J). The input is
validated but never normalized or otherwise modified.
}
\description{
Validate a Target PMF on \{1, ..., J\}
}
\keyword{internal}

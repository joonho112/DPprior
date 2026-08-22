% Please edit documentation in R/16_DPprior_fit.R
\name{vif_to_variance_fit}
\alias{vif_to_variance_fit}
\title{Convert a variance inflation factor to a K variance}
\usage{
vif_to_variance_fit(mu_K, vif)
}
\arguments{
\item{mu_K}{Target mean.}

\item{vif}{Variance inflation factor.}
}
\value{
vif times mu_K minus one.
}
\description{
Convert a variance inflation factor to a K variance
}
\keyword{internal}

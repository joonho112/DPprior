% Please edit documentation in R/17_s3_methods.R
\name{.dpprior_is_dual}
\alias{.dpprior_is_dual}
\title{Check whether a fit uses a dual-anchor mode}
\usage{
.dpprior_is_dual(fit)
}
\arguments{
\item{fit}{A DPprior_fit object.}
}
\value{
Logical; TRUE if dual-anchor fit with valid structure.
}
\description{
Schema-validates a DPprior_fit object and recognizes current hard, current
soft, and retained legacy dual modes.
}
\keyword{internal}

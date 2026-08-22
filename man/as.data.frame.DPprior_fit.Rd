% Please edit documentation in R/17_s3_methods.R
\name{as.data.frame.DPprior_fit}
\alias{as.data.frame.DPprior_fit}
\title{Coerce a canonical DPprior fit to a data frame}
\usage{
\method{as.data.frame}{DPprior_fit}(x, row.names = NULL, optional = FALSE, ...)
}
\arguments{
\item{x}{A canonical \code{DPprior_fit} object.}

\item{row.names}{Optional row names.}

\item{optional}{Retained for the data-frame S3 signature; it does not alter
the canonical column contract.}

\item{...}{Additional arguments, currently unused.}
}
\value{
A one-row data frame containing schema, mode, method, status,
\code{usable}, \code{verified}, candidate availability, Gamma parameters,
canonical target and achieved K moments, any named weight target evidence,
and status guidance. Parameter columns are \code{NA} when the canonical
result has no public candidate.
}
\description{
Returns a one-row, status-aware public view of a canonical
\code{dpprior.result/1} fit. The method validates the complete object before
extracting any value; it never repairs a result or reads unregistered flat
fields as scientific evidence.
}
\examples{
fit <- DPprior_fit(
  J = 50, mu_K = 5, var_K = 8, check_diagnostics = FALSE
)
as.data.frame(fit)

}
\seealso{
\code{\link{DPprior_fit}}, \code{\link{summary.DPprior_fit}}
}

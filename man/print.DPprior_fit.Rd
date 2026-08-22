% Please edit documentation in R/17_s3_methods.R
\name{print.DPprior_fit}
\alias{print.DPprior_fit}
\title{Print Method for DPprior_fit Objects}
\usage{
\method{print}{DPprior_fit}(x, digits = 4, ...)
}
\arguments{
\item{x}{A \code{DPprior_fit} object.}

\item{digits}{Integer; number of significant digits for display.
Default is 4.}

\item{...}{Additional arguments (currently unused).}
}
\value{
Invisibly returns \code{x}.
}
\description{
Prints the canonical schema, mode, method, decision fields, target, and
available candidate evidence. Migrated or failed objects without a public
candidate are labelled explicitly rather than exposing quarantined values.
}
\details{
The method first validates the complete \code{dpprior.result/1}
object. Printed convergence language is derived from \code{status},
\code{usable}, and \code{verified}; compatibility aliases and optimizer exit
codes are never treated as verification.
}
\examples{
# Create a fit object
fit <- DPprior_fit(J = 50, mu_K = 5, var_K = 8)
print(fit)

# With custom digits
print(fit, digits = 6)

}
\seealso{
\code{\link{summary.DPprior_fit}}, \code{\link{plot.DPprior_fit}},
\code{\link{DPprior_fit}}
}

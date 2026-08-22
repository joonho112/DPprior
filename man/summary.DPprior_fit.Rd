% Please edit documentation in R/17_s3_methods.R
\name{summary.DPprior_fit}
\alias{summary.DPprior_fit}
\title{Summary Method for DPprior_fit Objects}
\usage{
\method{summary}{DPprior_fit}(object, print_output = TRUE, ...)
}
\arguments{
\item{object}{A \code{DPprior_fit} object.}

\item{print_output}{Logical; if \code{TRUE} (default), prints the summary
to the console. If \code{FALSE}, returns the summary list silently.}

\item{...}{Additional arguments (currently unused).}
}
\value{
An object of class \code{"summary.DPprior_fit"} containing:
\itemize{
\item \code{schema}, \code{mode}, \code{method}, \code{status},
\code{usable}, and \code{verified};
\item \code{candidate_available}, \code{gamma_prior}, and derived
\code{alpha_summary};
\item canonical \code{target}, \code{achieved}, \code{weight}, and
\code{errors} views;
\item \code{message}, \code{guidance}, \code{iterations}, and attached
\code{diagnostics}.
}
}
\description{
Produces a status-aware summary from validated canonical fields, including
candidate availability, Gamma parameters, target-versus-achieved evidence,
and any attached diagnostics.
}
\details{
The complete input is schema-validated before summarization. A
migrated result without a public candidate remains visibly unavailable;
quarantined compatibility evidence is not promoted into the summary.
}
\examples{
# Create a fit object
fit <- DPprior_fit(J = 50, mu_K = 5, var_K = 8, check_diagnostics = TRUE)
summary(fit)

# Store summary without printing
summ <- summary(fit, print_output = FALSE)
str(summ)

}
\seealso{
\code{\link{print.DPprior_fit}}, \code{\link{DPprior_diagnostics}}
}

% Please edit documentation in R/DPprior-package.R
\docType{package}
\name{DPprior-package}
\alias{DPprior}
\alias{DPprior-package}
\title{DPprior: Decision-Ready Prior Elicitation for Dirichlet Process Mixtures}
\description{
\pkg{DPprior} implements Design-Conditional Elicitation (DCE) for a Gamma
shape-rate prior on the Dirichlet Process concentration parameter. Version
2 translates beliefs about the occupied-cluster count \eqn{K_J} and named
weight estimands into validated, auditable calibration results.
}
\details{
Start with \code{\link{DPprior_target_K}} when a reusable cluster-count
target is useful, then call \code{\link{DPprior_fit}}. The wrapper dispatches
to A1, A2-MN, or A2-KL and returns the numerical producer's native canonical
result. A1 is an explicitly approximate shifted-Negative-Binomial proxy;
A2-MN performs exact finite-design moment matching; A2-KL targets a complete
cluster-count distribution.

For an additional weight requirement, use
\code{\link{DPprior_dual_hard}} for a verified inequality or
\code{\link{DPprior_dual_soft}} for a fixed-scale trade-off. The hard and
soft modes have separate semantics and neither silently changes the
cluster-count target inherited from the input fit.
}
\section{Canonical object contract}{

Current producers return a \code{dpprior.result/1} object. Its common
top-level fields include \code{status}, \code{usable}, \code{verified},
\code{parameters}, \code{target}, \code{achieved}, \code{residuals},
\code{computation}, \code{verification}, and \code{provenance}. Scientific
consumers should use these canonical fields. Compatibility aliases are
explicitly non-authoritative migration views.

A result is decision-ready only when its mode-specific contract and the
public \code{status}/\code{usable}/\code{verified} fields permit that use.
Optimizer exit alone is never verification. When calibration cannot return
an approved result, a typed condition retains the complete canonical object
in \code{condition$result}; inspect that evidence before revising the target
or controls.
}

\section{Legacy retention policy}{

\code{\link{DPprior_dual}} is the retained historical equality-loss adapter.
It remains available throughout v2.x, will not be removed before v3.0, and
requires a migration review before any removal. It emits one lifecycle
warning and always labels its output legacy, approximate, and unverified;
it is not a hard-constraint certificate.
}

\section{Key function groups}{

\describe{
\item{Targets and fitting}{\code{\link{DPprior_target_K}},
\code{\link{DPprior_fit}}, \code{\link{DPprior_a1}},
\code{\link{DPprior_a2_newton}}, and \code{\link{DPprior_a2_kl}}.}
\item{Dual-anchor calibration}{\code{\link{DPprior_dual_hard}},
\code{\link{DPprior_dual_soft}}, and retained
\code{\link{DPprior_dual}}.}
\item{Migration}{\code{\link{upgrade_DPprior_object}} validates current
objects and conservatively upgrades exactly recognized 1.1 objects.}
\item{Diagnostics and visualization}{\code{\link{DPprior_diagnostics}},
\code{\link{DPprior_error_bounds}}, \code{\link{plot.DPprior_fit}}, and
\code{\link{summary.DPprior_fit}}.}
}
}

\section{Acknowledgments}{

This research was supported by the Institute of Education Sciences, U.S.
Department of Education, through Grant R305D240078 to the University of
Alabama.

The opinions expressed are those of the author and do not represent views of
the Institute or the U.S. Department of Education.
}

\references{
Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
\emph{arXiv preprint} arXiv:2602.06301.

Lee, J., Che, J., Rabe-Hesketh, S., Feller, A., and Miratrix, L. (2025).
Improving the estimation of site-specific effects and their distribution
in multisite trials.
\emph{Journal of Educational and Behavioral Statistics}, 50(5), 731--764.
\doi{10.3102/10769986241254286}
}
\seealso{
Useful links:
\itemize{
  \item \url{https://joonho112.github.io/DPprior/}
  \item \url{https://github.com/joonho112/DPprior}
  \item Report bugs at \url{https://github.com/joonho112/DPprior/issues}
}

}
\author{
\strong{Maintainer}: JoonHo Lee \email{jlee296@ua.edu} (\href{https://orcid.org/0009-0006-4019-8703}{ORCID})

Authors:
\itemize{
  \item JoonHo Lee \email{jlee296@ua.edu} (\href{https://orcid.org/0009-0006-4019-8703}{ORCID})
}

}
\keyword{internal}

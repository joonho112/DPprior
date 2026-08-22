% Please edit documentation in R/09_weights_rho.R
\name{summary_rho}
\alias{summary_rho}
\title{Summary Statistics for rho Distribution}
\usage{
summary_rho(
  a,
  b,
  M = .QUAD_NODES_DEFAULT,
  M_verify = NULL,
  abs_tol = 1e-10,
  rel_tol = 1e-08,
  strict = FALSE
)
}
\arguments{
\item{a}{Numeric; shape parameter of the Gamma prior on alpha (a > 0).}

\item{b}{Numeric; rate parameter of the Gamma prior on alpha (b > 0).}

\item{M}{Integer; selected quadrature order. Default is 80.}

\item{M_verify}{Optional independent verification order. When omitted and
available, the package-wide minimum refinement order is used.}

\item{abs_tol, rel_tol}{Non-negative selected-versus-verification tolerances.}

\item{strict}{Logical; if \code{TRUE}, require verified rho moments.}
}
\value{
A list of class "rho_summary" containing:
\describe{
\item{mean}{E(rho | a, b)}
\item{var}{Var(rho | a, b)}
\item{sd}{SD(rho | a, b) = sqrt(Var)}
\item{cv}{Coefficient of variation SD/mean}
\item{status, usable, verified, message}{Numerical verification state}
\item{estimand, label, units}{Explicit estimand metadata}
\item{params}{List of input parameters (a, b)}
\item{alpha_prior}{Summary of the alpha prior (mean, sd, cv)}
\item{conditional_at_alpha_mean}{Conditional moments evaluated at E(alpha)}
\item{verification, provenance}{Quadrature comparison and method metadata}
}
}
\description{
Computes comprehensive summary statistics for the co-clustering probability
rho under the hierarchical prior alpha ~ Gamma(a, b).
}
\details{
The co-clustering probability rho indicates how likely two randomly
chosen observations are to belong to the same cluster a priori.

The \code{conditional_at_alpha_mean} component provides a "plug-in" estimate
for comparison: what the moments would be if alpha were fixed at its prior
mean. The function reports numerical values and verification status only;
scientific action thresholds must be supplied explicitly by an analysis
policy rather than inferred from uncalibrated labels.
}
\examples{
summary_rho(a = 2, b = 1)
summary_rho(a = 1.6, b = 1.22)

}
\seealso{
\code{\link{mean_rho}}, \code{\link{var_rho}}, \code{\link{summary_w1}}
}

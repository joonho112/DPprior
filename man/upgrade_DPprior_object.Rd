% Please edit documentation in R/24_schema_upgrade.R
\name{upgrade_DPprior_object}
\alias{upgrade_DPprior_object}
\title{Upgrade a serialized DPprior object conservatively}
\usage{
upgrade_DPprior_object(x, verify = TRUE, M_verify = NULL, allow_legacy = FALSE)
}
\arguments{
\item{x}{A canonical object, or an exactly recognized DPprior 1.1 fit or
diagnostics object.}

\item{verify}{Logical. For fits, record a fresh fixed-candidate numerical
audit without changing the candidate or scientific status. Diagnostics
migration requires \code{TRUE} and recomputes all metrics from \code{a}, \code{b}, and
\code{J}.}

\item{M_verify}{Optional current verifier order used only for the fresh
migration audit or diagnostics recomputation.}

\item{allow_legacy}{Logical. When \code{TRUE}, only migrated A1 and retained
legacy Dual-Anchor fits may be usable approximations. It never makes an
A2, hard, or soft result decision-ready.}
}
\value{
A validated canonical object. Current \code{dpprior.result/1},
\code{dpprior.target/1}, and \code{dpprior.weight-target/1} inputs are
returned after validation. An exactly recognized 1.1 object is converted
conservatively and retains explicit migration and quarantine evidence;
it is never promoted to a verified current calibration merely because a
fixed-candidate audit succeeds. The input is not modified.
}
\description{
Converts an exactly recognized DPprior 1.1 fit or diagnostics object into
the canonical result schema. A migrated fit never inherits the former
\code{success} or \code{converged} claim: it remains approximate and unverified because
v1.1 did not serialize candidate-selection and independent-verifier lineage.
}
\details{
Objects that are not exactly recognized fail closed with a typed
serialization condition. Refit with the current API when decision-ready
calibration is required; compatibility views are not scientific fields.
}
\seealso{
\code{\link{DPprior_fit}}, \code{\link{DPprior_target_K}}
}

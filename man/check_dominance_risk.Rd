% Please edit documentation in R/14_diagnostics.R
\name{check_dominance_risk}
\alias{check_dominance_risk}
\title{Deprecated Ambiguous Weight-Risk Check}
\usage{
check_dominance_risk(a, b, threshold = 0.5, risk_level = 0.3)
}
\arguments{
\item{a, b}{Numeric Gamma hyperparameters.}

\item{threshold}{W_SB mass threshold in (0,1).}

\item{risk_level}{Action threshold on the probability scale.}
}
\value{
Logical policy result, with a typed deprecation warning.
}
\description{
This compatibility wrapper evaluates the explicit policy
\code{P(W_SB > threshold) > risk_level}; it does not evaluate W_max or
assign a category.
}

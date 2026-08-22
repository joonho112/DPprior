% Please edit documentation in R/22_dual_anchor_soft.R
\name{compute_tradeoff_curve}
\alias{compute_tradeoff_curve}
\title{Compute a status-aware soft trade-off path}
\usage{
compute_tradeoff_curve(
  J,
  K_target,
  w1_target = NULL,
  lambda_seq = seq(0.1, 1, by = 0.1),
  max_iter = 100L,
  M = .QUAD_NODES_DEFAULT,
  verbose = FALSE,
  loss_type = NULL,
  target = NULL,
  M_verify = NULL,
  log_bounds = .LOG_BOUNDS_DEFAULT,
  control = list(),
  allow_approximate = TRUE
)
}
\arguments{
\item{J}{Design size.}

\item{K_target}{List with \code{mu_K} and \code{var_K}.}

\item{w1_target}{Optional legacy W_SB target, translated to an explicit
equality target.}

\item{lambda_seq}{Unique values satisfying \code{0 < lambda <= 1}.}

\item{max_iter, M, M_verify, log_bounds, control}{Soft solver controls.}

\item{verbose}{Show path progress.}

\item{loss_type}{Deprecated curve argument. Only \code{NULL} is accepted;
the v2 curve always uses fixed input-derived scales.}

\item{target}{Preferred explicit named weight target.}

\item{allow_approximate}{Retain approximate points without making them
usable warm starts.}
}
\value{
A deterministic \code{dpprior_tradeoff_curve} data frame with
retained fit and condition lists as attributes. No lambda is selected as
best.
}
\description{
Compute a status-aware soft trade-off path
}

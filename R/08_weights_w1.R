# =============================================================================
# Module 08: First Size-Biased Weight and Largest-Weight Diagnostics
# =============================================================================
#
# This module implements the closed-form distribution of W_SB, the first
# stick-breaking weight in the size-biased GEM(α) representation, and a
# separate diagnostic for W_max, the largest ranked population weight.
#
# Key Results (Lee, 2026, Section 4; Vicentini & Jermyn, 2025):
#   Conditional: w₁ | α ~ Beta(1, α)
#   Unconditional CDF: F(x | a, b) = 1 - (b / (b - log(1-x)))^a
#   Quantile: Q(u | a, b) = 1 - exp(b × [1 - (1-u)^{-1/a}])
#   Density: p(x | a, b) = a × b^a / [(1-x) × (b - log(1-x))^{a+1}]
#
# Numerical Notes:
#   - Uses log1p/expm1 for stability near x ~ 0 and x ~ 1
#   - Boundary behavior handled explicitly:
#       CDF: x<=0 -> 0, x>=1 -> 1
#       Survival: t<=0 -> 1, t>=1 -> 0
#       Density: outside (0,1) -> 0 (or -Inf on log-scale)
#
# Author: JoonHo Lee (jlee296@ua.edu)
# Date: December 2025
# Part of: DPprior R Package
# Dependencies: Module 00 (constants.R), Module 02 (quadrature.R)
# =============================================================================


# =============================================================================
# CDF Function
# =============================================================================

# Stable log(1 + numerator / denominator) for non-negative numerator and a
# positive denominator. Forming the ratio first can overflow for valid
# subnormal Gamma rates even though the logarithm and final probability are
# representable.
.w1_log1p_positive_ratio <- function(numerator, denominator) {
  log_ratio <- log(numerator) - log(denominator)
  out <- numeric(length(log_ratio))
  numerator_dominates <- log_ratio > 0
  out[numerator_dominates] <- log_ratio[numerator_dominates] +
    log1p(exp(-log_ratio[numerator_dominates]))
  out[!numerator_dominates] <- log1p(exp(log_ratio[!numerator_dominates]))
  out
}

#' CDF of the First Size-Biased DP Weight
#'
#' Computes \eqn{P(w_1 \le x \mid a, b)} using the closed-form expression
#' derived by marginalizing over \eqn{\alpha \sim Gamma(a, b)}.
#'
#' @param x Numeric vector. Values outside the unit interval are allowed but
#'   are mapped to the boundary values of the CDF (0 for \eqn{x \le 0}, 1 for
#'   \eqn{x \ge 1}).
#' @param a Numeric; shape parameter of the Gamma prior on \eqn{\alpha}
#'   (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior on \eqn{\alpha}
#'   (b > 0).
#'
#' @return Numeric vector of CDF values F(x | a, b) with same length as x.
#'
#' @details
#' The unconditional CDF is given by:
#' \deqn{F_{w_1}(x | a, b) = 1 - \left(\frac{b}{b - \log(1-x)}\right)^a}
#'
#' The implementation uses \code{log1p} and \code{expm1} for numerical
#' stability, particularly when the CDF is close to 0 (small x).
#'
#' @section Interpretation:
#' The weight \eqn{w_1} is in **GEM (size-biased) order**, not ranked by size.
#' It represents the asymptotic cluster share of a randomly chosen unit,
#' **not** the largest cluster proportion. See Lee (2026, Section 4) for details.
#'
#' @examples
#' # P(w1 <= 0.3) under standard prior
#' cdf_w1(0.3, a = 2, b = 1)
#'
#' # Vectorized computation
#' cdf_w1(c(0.1, 0.3, 0.5, 0.7), a = 1.6, b = 1.22)
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' Vicentini, C. and Jermyn, I. H. (2025). Prior selection for the precision
#'   parameter of Dirichlet Process Mixtures. arXiv:2502.00864.
#'
#' @seealso \code{\link{quantile_w1}}, \code{\link{density_w1}},
#'   \code{\link{prob_w1_exceeds}}
#'
#' @family weights_w1
#'
#' @export
cdf_w1 <- function(x, a, b) {
  # Input validation
  a <- .dpprior_validate_scalar(
    a, "a", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  b <- .dpprior_validate_scalar(
    b, "b", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )

  # Initialize output with NA
  out <- rep(NA_real_, length(x))

  # Handle boundaries explicitly
  out[x <= 0] <- 0
  out[x >= 1] <- 1

  # Identify valid interior points

  idx <- which(is.finite(x) & x > 0 & x < 1)
  if (length(idx) == 0L) {
    return(out)
  }

  xx <- x[idx]

  # Numerically stable computation using log1p and expm1
  # log(1-x) via log1p for stability when x is small
  log1m <- log1p(-xx)        # log(1-x) <= 0
  # Stable log survival without subtracting nearly equal logarithms.
  log_surv <- -a * .w1_log1p_positive_ratio(-log1m, b)

  # CDF = 1 - exp(log_surv) via expm1 for stability when CDF is near 0
  out[idx] <- -expm1(log_surv)

  out
}


# =============================================================================
# Quantile Function
# =============================================================================

#' Quantile Function of the First Size-Biased DP Weight
#'
#' Computes the inverse CDF: \eqn{Q(u \mid a, b) = F^{-1}(u)}.
#'
#' @param u Numeric vector of probability levels in the unit interval.
#'   Values \eqn{u \le 0} return 0 and \eqn{u \ge 1} return 1.
#' @param a Numeric; shape parameter of the Gamma prior on \eqn{\alpha}
#'   (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior on \eqn{\alpha}
#'   (b > 0).
#'
#' @return Numeric vector of quantile values Q(u | a, b).
#'
#' @details
#' The quantile function has the closed form:
#' \deqn{Q_{w_1}(u | a, b) = 1 - \exp\left(b \left[1 - (1-u)^{-1/a}\right]\right)}
#'
#' The implementation computes (1-u)^(-1/a) in log space for stability
#' when u is close to 1.
#'
#' \strong{Numerical Note:} For small values of a (a < 1) and u close to 1,
#' the quantile approaches 1 very rapidly and may round to 1.0 in double
#' precision.
#'
#' @examples
#' # Median of w1
#' quantile_w1(0.5, a = 2, b = 1)  # ~0.339
#'
#' # 90th percentile
#' quantile_w1(0.9, a = 2, b = 1)  # ~0.732
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @seealso \code{\link{cdf_w1}}, \code{\link{summary_w1}}
#'
#' @family weights_w1
#'
#' @export
quantile_w1 <- function(u, a, b) {
  # Input validation
  a <- .dpprior_validate_scalar(
    a, "a", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  b <- .dpprior_validate_scalar(
    b, "b", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  assert_probability(u, "u")

  # Initialize output with NA
  out <- rep(NA_real_, length(u))

  # Handle boundaries explicitly
  out[u <= 0] <- 0
  out[u >= 1] <- 1

  # Identify valid interior points
  idx <- which(is.finite(u) & u > 0 & u < 1)
  if (length(idx) == 0L) {
    return(out)
  }

  uu <- u[idx]

  # Let L = -log(1-u)/a.  Evaluate log{b[exp(L)-1]} without first
  # forming exp(L): exp(L) can overflow even when multiplication by a tiny b
  # brings the final exponent back into the representable range.
  log_L <- log(-log1p(-uu)) - log(a)
  L <- exp(log_L)
  log_pow_minus_one <- numeric(length(L))
  small <- is.finite(L) & L <= log(2)
  log_pow_minus_one[small] <- ifelse(
    L[small] == 0,
    log_L[small],
    log(expm1(L[small]))
  )
  large <- is.finite(L) & !small
  log_pow_minus_one[large] <-
    L[large] + log1p(-exp(-L[large]))
  log_pow_minus_one[!is.finite(L)] <- Inf

  positive_exponent <- exp(log(b) + log_pow_minus_one)
  out[idx] <- -expm1(-positive_exponent)

  out
}


# =============================================================================
# Survival Function
# =============================================================================

#' Tail Probability of the First Size-Biased DP Weight
#'
#' Computes P(W_SB > t | a, b) = 1 - F(t), the probability that a randomly
#' selected unit belongs to a population cluster whose DP mass exceeds t.
#'
#' @param t Numeric vector of thresholds. Values outside the unit interval are
#'   allowed but are mapped to the boundary values (1 for \eqn{t \le 0}, 0 for
#'   \eqn{t \ge 1}).
#' @param a Numeric; shape parameter of the Gamma prior on \eqn{\alpha}
#'   (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior on \eqn{\alpha}
#'   (b > 0).
#'
#' @return Numeric vector of survival probabilities.
#'
#' @details
#' The survival function has the closed form:
#' \deqn{P(w_1 > t | a, b) = \left(\frac{b}{b - \log(1-t)}\right)^a}
#'
#' This is a mass-weighted cluster-size probability. It is not
#' \eqn{P(W_{max}>t)}, the probability that the largest population weight
#' exceeds \eqn{t}. Use \code{\link{prob_wmax_exceeds}} for that estimand.
#'
#' @examples
#' # Probability that a random unit's population cluster has mass above 0.5
#' prob_w1_exceeds(0.5, a = 1.6, b = 1.22)  # ~0.487 (Lee et al. DP-inform)
#' prob_w1_exceeds(0.5, a = 2, b = 1)       # ~0.349
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @seealso \code{\link{cdf_w1}}, \code{\link{summary_w1}}
#'
#' @family weights_w1
#'
#' @export
prob_w1_exceeds <- function(t, a, b) {
  # Input validation
  a <- .dpprior_validate_scalar(
    a, "a", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  b <- .dpprior_validate_scalar(
    b, "b", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )

  # Initialize output with NA
  out <- rep(NA_real_, length(t))

  # Handle boundaries explicitly
  out[t <= 0] <- 1
  out[t >= 1] <- 0

  # Identify valid interior points
  idx <- which(is.finite(t) & t > 0 & t < 1)
  if (length(idx) == 0L) {
    return(out)
  }

  tt <- t[idx]

  # Numerically stable computation
  c_t <- -log1p(-tt)
  log_surv <- -a * .w1_log1p_positive_ratio(c_t, b)
  out[idx] <- exp(log_surv)

  out
}


#' Tail Probability of the First Size-Biased DP Weight
#'
#' Explicitly named alias for \code{\link{prob_w1_exceeds}}. The legacy
#' \code{w1} name remains mathematically valid because the first GEM weight is
#' a size-biased pick, but \code{W_SB} makes the estimand unambiguous.
#'
#' @param threshold Numeric vector of thresholds.
#' @param a,b Positive scalar shape and rate of the Gamma prior on alpha.
#'
#' @return Numeric vector containing \eqn{P(W_{SB}>threshold)}.
#'
#' @examples
#' prob_wsb_exceeds(0.5, a = 2, b = 1)
#'
#' @seealso \code{\link{prob_w1_exceeds}},
#'   \code{\link{prob_wmax_exceeds}}
#'
#' @family weights_w1
#' @export
prob_wsb_exceeds <- function(threshold, a, b) {
  prob_w1_exceeds(threshold, a, b)
}


# =============================================================================
# Density Function
# =============================================================================

#' Density of the First Size-Biased DP Weight
#'
#' Computes the probability density \eqn{p(w_1 = x \mid a, b)}.
#'
#' @param x Numeric vector; evaluation points.
#' @param a Numeric; shape parameter of the Gamma prior on \eqn{\alpha}
#'   (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior on \eqn{\alpha}
#'   (b > 0).
#' @param log Logical; if \code{TRUE}, returns log-density. Default is \code{FALSE}.
#'
#' @return Numeric vector of density (or log-density if \code{log = TRUE}) values.
#'   Returns 0 (or -Inf on log scale) for x outside (0, 1).
#'
#' @details
#' The marginal density of \eqn{w_1} is:
#' \deqn{p(w_1 | a, b) = \frac{a \cdot b^a}{(1-w_1) \cdot [b - \log(1-w_1)]^{a+1}}}
#'
#' \strong{Important:} For small values of a (a < 1), the density has significant
#' mass concentrated very close to x = 1.
#'
#' @examples
#' # Density at several points
#' x <- seq(0.1, 0.9, by = 0.1)
#' density_w1(x, a = 2, b = 1)
#'
#' # Log-density for numerical stability
#' density_w1(0.5, a = 2, b = 1, log = TRUE)
#'
#' @seealso \code{\link{cdf_w1}}, \code{\link{quantile_w1}}
#'
#' @family weights_w1
#'
#' @export
density_w1 <- function(x, a, b, log = FALSE) {
  # Input validation
  a <- .dpprior_validate_scalar(
    a, "a", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  b <- .dpprior_validate_scalar(
    b, "b", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  log <- .dpprior_validate_control(log, "log", type = "logical")

  # Initialize output with NA
  out <- rep(NA_real_, length(x))

  # Handle boundaries explicitly
  if (log) {
    out[x <= 0] <- -Inf
    out[x >= 1] <- -Inf
  } else {
    out[x <= 0] <- 0
    out[x >= 1] <- 0
  }

  # Identify valid interior points
  idx <- which(is.finite(x) & x > 0 & x < 1)
  if (length(idx) == 0L) {
    return(out)
  }

  xx <- x[idx]

  # Numerically stable log-density computation
  log1m <- log1p(-xx)         # log(1-x)
  log_term <- -log1m          # -log(1-x)
  log_density <- log(a) - log(b) - log1m -
    (a + 1) * .w1_log1p_positive_ratio(log_term, b)

  if (log) {
    out[idx] <- log_density
  } else {
    out[idx] <- exp(log_density)
  }

  out
}


# =============================================================================
# Moment Functions
# =============================================================================

#' Mean of the First Size-Biased DP Weight
#'
#' Computes \eqn{E(w_1 \mid a, b)} via Gauss-Laguerre quadrature.
#'
#' @param a Numeric; shape parameter of the Gamma prior on \eqn{\alpha}
#'   (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior on \eqn{\alpha}
#'   (b > 0).
#' @param M Integer; number of quadrature nodes. Default is 80.
#'
#' @return Numeric; \eqn{E(w_1)}.
#'
#' @details
#' The expectation is computed using the identity:
#' \deqn{E[w_1 | a, b] = E\left[\frac{1}{1+\alpha}\right] = I_1(a, b)}
#'
#' where the integral is evaluated via Gauss-Laguerre quadrature.
#'
#' \strong{Key identity:} \eqn{E(w_1 \mid a, b) = E(\rho \mid a, b)}, where
#' \eqn{\rho = \sum_h w_h^2} is the co-clustering probability.
#'
#' @examples
#' mean_w1(a = 2, b = 1)       # ~0.404
#' mean_w1(a = 1.6, b = 1.22)  # ~0.508
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @seealso \code{\link{var_w1}}, \code{\link{summary_w1}}
#'
#' @family weights_w1
#'
#' @export
mean_w1 <- function(a, b, M = .QUAD_NODES_DEFAULT) {
  # Input validation
  a <- .dpprior_validate_scalar(
    a, "a", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  b <- .dpprior_validate_scalar(
    b, "b", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )

  # E[w₁] = E[1/(1+α)] where α ~ Gamma(a, b)
  integrate_gamma(function(alpha) 1 / (1 + alpha), a, b, M)
}


#' Variance of the First Size-Biased DP Weight
#'
#' Computes \eqn{Var(w_1 \mid a, b)} using the law of total variance.
#'
#' @param a Numeric; shape parameter of the Gamma prior on \eqn{\alpha}
#'   (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior on \eqn{\alpha}
#'   (b > 0).
#' @param M Integer; number of quadrature nodes. Default is 80.
#'
#' @return Numeric; \eqn{Var(w_1)}.
#'
#' @details
#' Uses the law of total variance:
#' \deqn{Var(w_1) = E[Var(w_1 | \alpha)] + Var(E[w_1 | \alpha])}
#'
#' where \eqn{w_1 \mid \alpha \sim Beta(1, \alpha)}, so:
#' \itemize{
#'   \item \eqn{E(w_1 \mid \alpha) = 1/(1+\alpha)}
#'   \item \eqn{Var(w_1 \mid \alpha) = \alpha / ((1+\alpha)^2(2+\alpha))}
#' }
#'
#' @examples
#' var_w1(a = 2, b = 1)       # ~0.090
#' var_w1(a = 1.6, b = 1.22)  # ~0.105
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @seealso \code{\link{mean_w1}}, \code{\link{summary_w1}}
#'
#' @family weights_w1
#'
#' @export
var_w1 <- function(a, b, M = .QUAD_NODES_DEFAULT) {
  # Input validation
  a <- .dpprior_validate_scalar(
    a, "a", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  b <- .dpprior_validate_scalar(
    b, "b", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )

  # Center the between-alpha component directly. Subtracting
  # E[m(alpha)^2] - E[m(alpha)]^2 can erase this small non-negative term for a
  # concentrated Gamma prior and can even produce a negative total variance.
  quad <- build_gamma_quadrature(a, b, M)
  weights <- quad$weights_normalized
  alpha <- quad$alpha_nodes
  conditional_mean <- 1 / (1 + alpha)
  marginal_mean <- .quadrature_weighted_sum(weights, conditional_mean)
  within_alpha <- .quadrature_weighted_sum(
    weights,
    (alpha / (1 + alpha)) * conditional_mean * (1 / (2 + alpha))
  )
  between_alpha <- .quadrature_weighted_sum(
    weights, (conditional_mean - marginal_mean)^2
  )

  as.numeric(within_alpha + between_alpha)
}


# =============================================================================
# Summary Function
# =============================================================================

#' Summary Statistics for the First Size-Biased DP Weight
#'
#' Computes comprehensive summary statistics for the \eqn{w_1} distribution.
#'
#' @param a Numeric; shape parameter of the Gamma prior on \eqn{\alpha}
#'   (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior on \eqn{\alpha}
#'   (b > 0).
#' @param probs Numeric vector; quantile probabilities. Default is
#'   \code{c(0.05, 0.25, 0.5, 0.75, 0.95)}.
#' @param M Integer; number of quadrature nodes for mean/variance. Default is 80.
#'
#' @return A list of class "w1_summary" containing:
#'   \describe{
#'     \item{mean}{\eqn{E(w_1)}}
#'     \item{var}{\eqn{Var(w_1)}}
#'     \item{sd}{\eqn{SD(w_1) = \sqrt{Var(w_1)}}}
#'     \item{median}{Median of \eqn{w_1}}
#'     \item{quantiles}{Named vector of quantiles}
#'     \item{prob_gt_50}{P(W_SB > 0.5), a size-biased weight tail}
#'     \item{prob_gt_90}{P(W_SB > 0.9), a size-biased weight tail}
#'     \item{params}{List of input parameters (a, b)}
#'     \item{estimand,label,conditioning,provenance}{Explicit scientific and
#'       computational metadata}
#'   }
#'
#' @examples
#' # Standard summary
#' summary_w1(a = 2, b = 1)
#'
#' # Lee et al. DP-inform prior
#' summary_w1(a = 1.6, b = 1.22)
#'
#' @seealso \code{\link{cdf_w1}}, \code{\link{quantile_w1}}, \code{\link{mean_w1}}
#'
#' @family weights_w1
#'
#' @export
summary_w1 <- function(a, b, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                       M = .QUAD_NODES_DEFAULT) {
  # Input validation
  a <- .dpprior_validate_scalar(
    a, "a", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  b <- .dpprior_validate_scalar(
    b, "b", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  assert_probability(probs, "probs")

  # Compute quantiles
  quantiles <- quantile_w1(probs, a, b)
  names(quantiles) <- paste0("q", 100 * probs)

  # Compute moments
  m <- mean_w1(a, b, M)
  v <- var_w1(a, b, M)

  result <- list(
    mean = m,
    var = v,
    sd = sqrt(v),
    median = unname(quantiles["q50"]),
    quantiles = quantiles,
    prob_gt_50 = prob_w1_exceeds(0.5, a, b),
    prob_gt_90 = prob_w1_exceeds(0.9, a, b),
    params = list(a = a, b = b),
    estimand = "W_SB",
    label = "First size-biased DP weight",
    conditioning = "gamma_mixed",
    provenance = list(
      schema_version = 1L,
      gamma_parameterization = "shape_rate",
      distribution_method = "closed_form_beta_gamma_mixture",
      moment_method = "generalized_gauss_laguerre",
      M = as.integer(M)
    )
  )

  class(result) <- "w1_summary"
  result
}


#' Print Method for w1_summary Objects
#'
#' @param x An object of class "w1_summary".
#' @param digits Integer; number of digits for printing.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.w1_summary <- function(x, digits = 4, ...) {
  cat("First size-biased DP weight (W_SB) summary\n")
  cat(strrep("=", 45), "\n\n")

  cat(sprintf("Gamma prior: alpha ~ Gamma(%.4f, %.4f)\n",
              x$params$a, x$params$b))
  cat(sprintf("E[alpha] = %.4f, CV(alpha) = %.2f%%\n\n",
              x$params$a / x$params$b, 100 / sqrt(x$params$a)))

  cat("Location and Scale:\n")
  cat(strrep("-", 30), "\n")
  cat(sprintf("  Mean:   %.*f\n", digits, x$mean))
  cat(sprintf("  Median: %.*f\n", digits, x$median))
  cat(sprintf("  SD:     %.*f\n", digits, x$sd))

  cat("\nQuantiles:\n")
  cat(strrep("-", 30), "\n")
  q_str <- paste(sprintf("  %s: %.*f", names(x$quantiles), digits, x$quantiles),
                 collapse = "\n")
  cat(q_str, "\n")

  cat("\nSize-biased cluster-mass tails:\n")
  cat(strrep("-", 30), "\n")
  cat(sprintf("  P(W_SB > 0.5): %.*f\n", digits, x$prob_gt_50))
  cat(sprintf("  P(W_SB > 0.9): %.*f\n", digits, x$prob_gt_90))

  invisible(x)
}


# =============================================================================
# Utility: Random Generation (for Monte Carlo validation)
# =============================================================================

#' Random Generation from the First Size-Biased DP Weight Distribution
#'
#' Generates random samples from the \eqn{w_1} distribution by first sampling
#' \eqn{\alpha \sim Gamma(a, b)}, then
#' \eqn{w_1 \mid \alpha \sim Beta(1, \alpha)}.
#'
#' @param n Integer; number of samples to generate.
#' @param a Numeric; shape parameter of the Gamma prior on \eqn{\alpha}
#'   (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior on \eqn{\alpha}
#'   (b > 0).
#'
#' @return Numeric vector of length n; random samples from the \eqn{w_1}
#'   distribution.
#'
#' @details
#' This function uses the hierarchical representation:
#' \enumerate{
#'   \item \eqn{\alpha \sim Gamma(a, b)}
#'   \item \eqn{w_1 \mid \alpha \sim Beta(1, \alpha)}
#' }
#'
#' Useful for Monte Carlo validation of the closed-form functions.
#'
#' @examples
#' # Generate samples
#' set.seed(42)
#' samples <- rw1(10000, a = 2, b = 1)
#'
#' # Compare empirical vs theoretical mean
#' mean(samples)          # ~0.404
#' mean_w1(a = 2, b = 1)  # 0.4037
#'
#' @export
rw1 <- function(n, a, b) {
  # Input validation
  a <- .dpprior_validate_scalar(
    a, "a", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  b <- .dpprior_validate_scalar(
    b, "b", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  n <- .dpprior_validate_count(
    n, "n", minimum = 1L,
    .subclass = "dpprior_weight_random_error"
  )

  # Hierarchical sampling: α ~ Gamma(a, b), then w₁ | α ~ Beta(1, α)
  alpha <- stats::rgamma(n, shape = a, rate = b)
  if (any(!is.finite(alpha)) || any(alpha <= 0)) {
    stop(.dpprior_new_condition(
      "rw1 generated a non-finite or non-positive alpha draw",
      classes = c(
        "dpprior_weight_rng_error", "dpprior_numerical_error",
        "dpprior_error", "error"
      ),
      alpha = alpha,
      gamma_shape = a,
      gamma_rate = b,
      code = "invalid_alpha_draw"
    ))
  }
  draws <- stats::rbeta(n, shape1 = 1, shape2 = alpha)
  if (any(!is.finite(draws)) || any(draws < 0) || any(draws > 1)) {
    stop(.dpprior_new_condition(
      "rw1 generated an invalid Beta draw",
      classes = c(
        "dpprior_weight_rng_error", "dpprior_numerical_error",
        "dpprior_error", "error"
      ),
      alpha = alpha,
      draws = draws,
      code = "invalid_beta_draw"
    ))
  }
  draws
}


# =============================================================================
# Utility: Grid Computation for Plotting
# =============================================================================

#' Compute the First Size-Biased DP Weight Distribution on a Grid
#'
#' Computes CDF, PDF, and survival function on a grid of x values.
#' Useful for visualization and comparison across different priors.
#'
#' @param a Numeric; shape parameter of the Gamma prior on \eqn{\alpha}.
#' @param b Numeric; rate parameter of the Gamma prior on \eqn{\alpha}.
#' @param x_grid Numeric vector; grid of x values in (0, 1).
#'   Default is \code{seq(0.01, 0.99, length.out = 100)}.
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{x}{Grid points}
#'     \item{cdf}{CDF values F(x)}
#'     \item{pdf}{Density values p(x)}
#'     \item{survival}{Survival function S(x) = 1 - F(x)}
#'   }
#'
#' @examples
#' # Compute on default grid
#' df <- w1_grid(a = 2, b = 1)
#'
#' # Plot all three functions
#' par(mfrow = c(1, 3))
#' plot(df$x, df$cdf, type = "l", main = "CDF", xlab = "x", ylab = "F(x)")
#' plot(df$x, df$pdf, type = "l", main = "PDF", xlab = "x", ylab = "p(x)")
#' plot(df$x, df$survival, type = "l", main = "Survival", xlab = "x", ylab = "S(x)")
#'
#' @export
w1_grid <- function(a, b, x_grid = seq(0.01, 0.99, length.out = 100)) {
  a <- .dpprior_validate_scalar(
    a, "a", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  b <- .dpprior_validate_scalar(
    b, "b", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )

  result <- data.frame(
    x = x_grid,
    cdf = cdf_w1(x_grid, a, b),
    pdf = density_w1(x_grid, a, b),
    survival = prob_w1_exceeds(x_grid, a, b)
  )
  attr(result, "estimand") <- "W_SB"
  attr(result, "label") <- "First size-biased DP weight"
  result
}


# =============================================================================
# Largest Population Weight (W_max)
# =============================================================================

.wmax_validate_specification <- function(alpha, a, b) {
  fixed <- !is.null(alpha)
  mixed_any <- !is.null(a) || !is.null(b)

  if (fixed && mixed_any) {
    .dpprior_abort_invalid(
      paste(
        "specify either fixed alpha or Gamma shape-rate parameters a and b,",
        "not both"
      ),
      "dpprior_weight_specification_error", "alpha/a/b",
      list(alpha = alpha, a = a, b = b),
      "exactly one of alpha or the pair (a,b)", "ambiguous_specification"
    )
  }
  if (!fixed && (is.null(a) || is.null(b))) {
    .dpprior_abort_invalid(
      "provide fixed alpha or both Gamma shape a and rate b",
      "dpprior_weight_specification_error", "alpha/a/b",
      list(alpha = alpha, a = a, b = b),
      "fixed alpha or complete pair (a,b)", "incomplete_specification"
    )
  }

  if (fixed) {
    alpha <- .dpprior_validate_scalar(
      alpha, "alpha", lower = 0, lower_open = TRUE,
      .subclass = "dpprior_weight_parameter_error"
    )
    return(list(
      conditioning = "fixed_alpha",
      alpha = alpha,
      a = NULL,
      b = NULL,
      params = list(alpha = alpha)
    ))
  }

  a <- .dpprior_validate_scalar(
    a, "a", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  b <- .dpprior_validate_scalar(
    b, "b", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  list(
    conditioning = "gamma_mixed",
    alpha = NULL,
    a = a,
    b = b,
    params = list(a = a, b = b)
  )
}


.wmax_validate_method <- function(method) {
  choices <- c("auto", "deterministic", "bounds", "monte_carlo")
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
      !(method %in% choices)) {
    .dpprior_abort_invalid(
      sprintf("method must be one of: %s", paste(choices, collapse = ", ")),
      c("dpprior_weight_method_error", "dpprior_type_error"),
      "method", method, paste(choices, collapse = ", "), "unsupported_method"
    )
  }
  method
}


.wmax_log_size_biased_tail <- function(threshold, specification) {
  if (threshold == 0) {
    return(0)
  }
  if (threshold == 1) {
    return(-Inf)
  }
  if (identical(specification$conditioning, "fixed_alpha")) {
    return(specification$alpha * log1p(-threshold))
  }
  c_t <- -log1p(-threshold)
  -specification$a * .w1_log1p_positive_ratio(c_t, specification$b)
}

.wmax_probability_next_down <- function(value) {
  if (!is.finite(value) || value <= 0) return(0)
  min_positive <- .Machine$double.xmin * .Machine$double.eps
  decrement <- max(abs(value) * .Machine$double.eps, min_positive)
  candidate <- value - decrement
  if (identical(candidate, value)) {
    candidate <- value * (1 - 2 * .Machine$double.eps)
  }
  max(0, candidate)
}

.wmax_probability_next_up <- function(value) {
  if (!is.finite(value) || value >= 1) return(1)
  min_positive <- .Machine$double.xmin * .Machine$double.eps
  increment <- max(abs(value) * .Machine$double.eps, min_positive)
  candidate <- value + increment
  if (identical(candidate, value)) {
    candidate <- value * (1 + 2 * .Machine$double.eps)
  }
  min(1, candidate)
}


#' Certified Bounds for the Largest DP Population Weight Tail
#'
#' Computes distribution-free bounds on
#' \eqn{P(W_{max}>threshold)}. If \eqn{q_t=P(W_{SB}>threshold)}, then
#' \deqn{q_t \le P(W_{max}>threshold) \le \min(q_t/threshold,1).}
#' The identity remains valid after mixing alpha over a Gamma hyperprior.
#'
#' @param threshold One finite probability in the closed unit interval.
#' @param alpha Optional positive scalar fixed DP concentration.
#' @param a,b Optional positive scalar Gamma shape and rate. Supply either
#'   \code{alpha} or both \code{a} and \code{b}, never both specifications.
#'
#' @return A \code{wmax_tail_bounds} object. The fields
#'   \code{lower_bound} and \code{upper_bound} are certified bounds, not point
#'   estimates. \code{size_biased_tail} is \eqn{P(W_{SB}>threshold)}.
#'
#' @examples
#' wmax_tail_bounds(0.4, alpha = 1)
#' wmax_tail_bounds(0.5, a = 2, b = 1)
#'
#' @seealso \code{\link{prob_wmax_exceeds}},
#'   \code{\link{prob_wsb_exceeds}}
#'
#' @family weights_w1
#' @export
wmax_tail_bounds <- function(threshold, alpha = NULL, a = NULL, b = NULL) {
  threshold <- .dpprior_validate_probability(
    threshold, "threshold", scalar = TRUE, open = FALSE,
    .subclass = "dpprior_weight_threshold_error"
  )
  specification <- .wmax_validate_specification(alpha, a, b)
  log_q_t <- .wmax_log_size_biased_tail(threshold, specification)
  q_t <- exp(log_q_t)
  lower <- q_t
  # exp()/log() round trips cannot reveal every upward rounding event. Give
  # every finite interior lower endpoint one representable downward slack;
  # exact logarithmic endpoints remain available for audit.
  lower_outward_rounded <- threshold > 0 && threshold < 1 && q_t > 0
  if (lower_outward_rounded) {
    lower <- .wmax_probability_next_down(q_t)
  }

  log_upper <- if (threshold == 0) {
    0
  } else {
    min(log_q_t - log(threshold), 0)
  }
  upper <- exp(log_upper)
  min_positive <- .Machine$double.xmin * .Machine$double.eps
  probability_scale_underflow <- is.finite(log_q_t) && q_t == 0
  upper_scale_underflow <- is.finite(log_upper) && upper == 0
  # Likewise, give every non-boundary finite upper endpoint upward slack.
  upper_outward_rounded <- threshold > 0 && threshold < 1 &&
    upper > 0 && upper < 1
  if (upper_scale_underflow) {
    # The exact upper bound is strictly positive but below ordinary double
    # precision. A smallest-positive ceiling remains conservative; the exact
    # logarithmic bound is retained separately.
    upper <- min_positive
  } else if (upper_outward_rounded) {
    upper <- .wmax_probability_next_up(upper)
  }

  result <- list(
    estimand = "W_max",
    label = "Largest DP population weight tail",
    threshold = threshold,
    conditioning = specification$conditioning,
    size_biased_tail = q_t,
    log_size_biased_tail = log_q_t,
    lower_bound = lower,
    upper_bound = upper,
    log_lower_bound = log_q_t,
    log_upper_bound = log_upper,
    probability_scale_underflow = probability_scale_underflow,
    upper_scale_underflow = upper_scale_underflow,
    lower_outward_rounded = lower_outward_rounded,
    upper_outward_rounded = upper_outward_rounded,
    certified = TRUE,
    method = "size_biased_mass_identity",
    units = "probability",
    params = specification$params,
    provenance = list(
      schema_version = 1L,
      bounds_identity = "q_t <= P(W_max > t) <= min(q_t/t, 1)",
      size_biased_estimand = "W_SB",
      probability_representation = if (
        probability_scale_underflow || upper_scale_underflow
      ) {
        "log_bounds_with_conservative_smallest_positive_real_ceiling"
      } else if (lower_outward_rounded || upper_outward_rounded) {
        "outward_rounded_ordinary_bounds_with_exact_log_bounds"
      } else {
        "ordinary_double_and_log"
      },
      gamma_parameterization = if (
        identical(specification$conditioning, "gamma_mixed")
      ) "shape_rate" else NULL
    )
  )
  class(result) <- "wmax_tail_bounds"
  result
}


#' @export
print.wmax_tail_bounds <- function(x, digits = 6, ...) {
  conditioning <- if (identical(x$conditioning, "fixed_alpha")) {
    sprintf("fixed alpha = %.*g", digits, x$params$alpha)
  } else {
    sprintf(
      "Gamma(shape = %.*g, rate = %.*g)",
      digits, x$params$a, digits, x$params$b
    )
  }
  cat("Certified bounds for P(W_max > threshold)\n")
  cat(sprintf("  Conditioning: %s\n", conditioning))
  cat(sprintf("  Threshold: %.*g\n", digits, x$threshold))
  cat(sprintf(
    "  Rounded P(W_SB > threshold): %.*g\n",
    digits, x$size_biased_tail
  ))
  cat(sprintf(
    "  Certified bounds: [%.*g, %.*g]\n",
    digits, x$lower_bound, digits, x$upper_bound
  ))
  if (isTRUE(x$probability_scale_underflow) ||
      isTRUE(x$upper_scale_underflow)) {
    cat("  Note: ordinary-scale underflow; exact bounds are retained on the log scale.\n")
  }
  invisible(x)
}


.wmax_fixed_series_enclosure <- function(alpha, threshold,
                                          abs_tol = 1e-14,
                                          max_terms = 10000L) {
  z <- 1 - threshold
  log_z <- log1p(-threshold)
  lead <- exp(alpha * log_z)
  if (lead == 0 || z == 0) {
    return(list(
      lower = 0, upper = 0, midpoint = 0,
      remainder_bound = 0, terms = 1L, passed = TRUE
    ))
  }

  partial <- 1
  remainder <- Inf
  j <- 0L
  while (j < max_terms && remainder > abs_tol) {
    j <- j + 1L
    partial <- partial + alpha * exp(j * log_z) / (alpha + j)
    remainder <- lead * alpha * exp((j + 1) * log_z) /
      ((alpha + j + 1) * (1 - z))
  }

  lower <- lead * partial
  roundoff <- .Machine$double.eps * max(1, abs(lower)) * (j + 4)
  lower_safe <- max(0, lower - roundoff)
  upper_safe <- min(1, lower + remainder + roundoff)
  list(
    lower = lower_safe,
    upper = upper_safe,
    midpoint = (lower_safe + upper_safe) / 2,
    remainder_bound = remainder,
    roundoff_allowance = roundoff,
    terms = as.integer(j + 1L),
    passed = is.finite(remainder) && remainder <= abs_tol
  )
}


.wmax_integrate_fixed <- function(alpha, threshold, rel_tol, abs_tol,
                                  subdivisions) {
  c_t <- -log1p(-threshold)
  integrand <- function(c_value) {
    exp(log(alpha) - alpha * c_value - log(-expm1(-c_value)))
  }
  tryCatch(
    stats::integrate(
      integrand, lower = c_t, upper = Inf,
      rel.tol = rel_tol, abs.tol = abs_tol,
      subdivisions = subdivisions, stop.on.error = FALSE
    ),
    error = function(error) error
  )
}


.wmax_integrate_gamma <- function(a, b, threshold, rel_tol, abs_tol,
                                  subdivisions) {
  c_t <- -log1p(-threshold)
  integrand <- function(c_value) {
    log_value <- log(a) - log(b + c_value) -
      a * .w1_log1p_positive_ratio(c_value, b) -
      log(-expm1(-c_value))
    exp(log_value)
  }
  tryCatch(
    stats::integrate(
      integrand, lower = c_t, upper = Inf,
      rel.tol = rel_tol, abs.tol = abs_tol,
      subdivisions = subdivisions, stop.on.error = FALSE
    ),
    error = function(error) error
  )
}


.wmax_result_template <- function(bounds) {
  list(
    estimand = "W_max",
    label = "Largest DP population weight tail",
    threshold = bounds$threshold,
    conditioning = bounds$conditioning,
    estimate = NA_real_,
    lower = NA_real_,
    upper = NA_real_,
    abs_error_bound = NA_real_,
    lower_bound = bounds$lower_bound,
    upper_bound = bounds$upper_bound,
    size_biased_tail = bounds$size_biased_tail,
    log_size_biased_tail = bounds$log_size_biased_tail,
    verified = FALSE,
    estimate_usable = FALSE,
    status = "failed",
    usable = FALSE,
    bounds_usable = TRUE,
    method = NA_character_,
    reason = NA_character_,
    units = "probability",
    numerical = NULL,
    sampling = NULL,
    provenance = list(
      schema_version = 1L,
      bounds_method = bounds$method,
      bounds_certified = TRUE,
      probability_scale_underflow = bounds$probability_scale_underflow,
      general_threshold_exact_method_deferred = FALSE,
      gamma_parameterization = bounds$provenance$gamma_parameterization
    ),
    params = bounds$params
  )
}


.wmax_bounds_only_result <- function(bounds, reason) {
  result <- .wmax_result_template(bounds)
  result$status <- "approximate"
  # No direct tail estimate is present. Keep result-level usability false even
  # for an explicit bounds request; the certified bounds remain independently
  # available through `bounds_usable` and `wmax_tail_bounds()`.
  result$usable <- FALSE
  result$method <- "certified_bounds_only"
  result$reason <- reason
  result$provenance$general_threshold_exact_method_deferred <-
    bounds$threshold > 0 && bounds$threshold < 0.5
  class(result) <- "wmax_tail_result"
  result
}


.wmax_boundary_result <- function(bounds, requested_method = "auto") {
  result <- .wmax_result_template(bounds)
  value <- if (bounds$threshold == 0) 1 else 0
  result$estimate <- value
  result$lower <- value
  result$upper <- value
  result$abs_error_bound <- 0
  result$verified <- TRUE
  result$estimate_usable <- TRUE
  result$status <- "converged"
  result$usable <- TRUE
  result$method <- "probability_boundary_identity"
  result$reason <- "exact_boundary"
  result$numerical <- list(reported_abs_error = 0)
  result$provenance$requested_method <- requested_method
  result$provenance$method_short_circuit <-
    !identical(requested_method, "auto")
  class(result) <- "wmax_tail_result"
  result
}


.wmax_deterministic_result <- function(bounds, specification, rel_tol,
                                       abs_tol, subdivisions) {
  result <- .wmax_result_template(bounds)
  verification_epsilon <- NA_real_

  if (identical(specification$conditioning, "fixed_alpha")) {
    primary <- .wmax_integrate_fixed(
      specification$alpha, bounds$threshold,
      rel_tol, abs_tol, subdivisions
    )
    series <- .wmax_fixed_series_enclosure(
      specification$alpha, bounds$threshold,
      abs_tol = max(.Machine$double.eps, min(abs_tol / 10, 1e-14))
    )
    result$method <-
      "one_dimensional_quadrature_verified_by_positive_series"

    if (inherits(primary, "error")) {
      result$reason <- "quadrature_error"
      result$numerical <- list(error = conditionMessage(primary))
      class(result) <- "wmax_tail_result"
      return(result)
    }

    reported_error <- if (is.finite(primary$abs.error)) {
      primary$abs.error
    } else {
      Inf
    }
    verification_epsilon <- max(1e-10, 10 * reported_error)
    agreement <- is.finite(primary$value) && is.finite(reported_error) &&
      series$passed &&
      primary$value >= series$lower - verification_epsilon &&
      primary$value <= series$upper + verification_epsilon
    abs_error_bound <- max(
      reported_error,
      abs(primary$value - series$lower),
      abs(series$upper - primary$value)
    )
    independent_pass <- agreement
    independent_difference <- abs(primary$value - series$midpoint)
    numerical <- list(
      primary_method = "transformed_one_dimensional_integral",
      reported_abs_error = reported_error,
      error_statement_type = paste(
        "reported_quadrature_error_plus_distance_to",
        "positive_series_enclosure",
        sep = "_"
      ),
      subdivisions_used = as.integer(primary$subdivisions),
      integration_message = primary$message,
      independent_method = "conditional_positive_series_enclosure",
      independent_lower = series$lower,
      independent_upper = series$upper,
      independent_difference = independent_difference,
      independent_tolerance = verification_epsilon,
      independent_passed = independent_pass,
      series_terms = series$terms,
      series_remainder_bound = series$remainder_bound
    )
  } else {
    primary <- .wmax_integrate_gamma(
      specification$a, specification$b, bounds$threshold,
      rel_tol, abs_tol, subdivisions
    )
    result$method <- paste(
      "one_dimensional_quadrature",
      "verified_by_gamma_quadrature_of_conditional_series",
      sep = "_"
    )
    if (inherits(primary, "error")) {
      result$reason <- "quadrature_error"
      result$numerical <- list(error = conditionMessage(primary))
      class(result) <- "wmax_tail_result"
      return(result)
    }

    conditional_series <- function(alpha_value) {
      vapply(
        alpha_value,
        function(one_alpha) .wmax_fixed_series_enclosure(
          one_alpha, bounds$threshold, abs_tol = 1e-14
        )$midpoint,
        numeric(1)
      )
    }
    independent <- tryCatch(
      c(
        selected = integrate_gamma(
          conditional_series, specification$a, specification$b, M = 256L
        ),
        verification = integrate_gamma(
          conditional_series, specification$a, specification$b, M = 512L
        )
      ),
      error = function(error) error
    )
    if (inherits(independent, "error")) {
      result$reason <- "independent_verification_error"
      result$numerical <- list(
        candidate = primary$value,
        reported_abs_error = primary$abs.error,
        integration_message = primary$message,
        independent_error = conditionMessage(independent)
      )
      class(result) <- "wmax_tail_result"
      return(result)
    }

    reported_error <- if (is.finite(primary$abs.error)) {
      primary$abs.error
    } else {
      Inf
    }
    order_difference <- abs(independent[["selected"]] -
                              independent[["verification"]])
    order_tolerance <- abs_tol + rel_tol * max(abs(independent))
    verification_epsilon <- max(1e-10, 10 * reported_error)
    independent_difference <- abs(
      primary$value - independent[["verification"]]
    )
    independent_pass <- all(is.finite(independent)) &&
      is.finite(reported_error) && is.finite(order_difference) &&
      is.finite(order_tolerance) && is.finite(independent_difference) &&
      order_difference <= order_tolerance &&
      independent_difference <= verification_epsilon
    abs_error_bound <- max(
      reported_error, order_difference, independent_difference
    )
    numerical <- list(
      primary_method = "gamma_mixed_transformed_one_dimensional_integral",
      reported_abs_error = reported_error,
      error_statement_type = paste(
        "reported_quadrature_error_plus_independent",
        "selected_refined_discrepancy",
        sep = "_"
      ),
      subdivisions_used = as.integer(primary$subdivisions),
      integration_message = primary$message,
      independent_method =
        "gamma_quadrature_of_conditional_positive_series",
      M_selected = 256L,
      M_verification = 512L,
      M_verification_required = 512L,
      independent_selected = unname(independent[["selected"]]),
      independent_verification = unname(independent[["verification"]]),
      independent_order_difference = order_difference,
      independent_order_tolerance = order_tolerance,
      independent_difference = independent_difference,
      independent_tolerance = verification_epsilon,
      independent_passed = independent_pass
    )
  }

  value_finite <- is.finite(primary$value) && is.finite(abs_error_bound)
  within_probability <- value_finite && primary$value >= 0 && primary$value <= 1
  within_certified_bounds <- value_finite &&
    primary$value >= bounds$lower_bound - verification_epsilon &&
    primary$value <= bounds$upper_bound + verification_epsilon
  error_pass <- value_finite && abs_error_bound <= 1e-8
  integration_pass <- identical(primary$message, "OK")
  representation_pass <- !isTRUE(bounds$probability_scale_underflow) &&
    !isTRUE(bounds$upper_scale_underflow)
  verified <- integration_pass && independent_pass && within_probability &&
    within_certified_bounds && error_pass && representation_pass

  candidate_lower <- if (value_finite) {
    max(bounds$lower_bound, 0, primary$value - abs_error_bound)
  } else {
    NA_real_
  }
  candidate_upper <- if (value_finite) {
    min(bounds$upper_bound, 1, primary$value + abs_error_bound)
  } else {
    NA_real_
  }
  interval_contains_estimate <- value_finite &&
    candidate_lower <= primary$value && primary$value <= candidate_upper
  verified <- verified && interval_contains_estimate
  if (verified) {
    result$estimate <- primary$value
    result$abs_error_bound <- abs_error_bound
    result$lower <- candidate_lower
    result$upper <- candidate_upper
  }
  result$verified <- verified
  result$estimate_usable <- verified
  result$status <- if (value_finite) {
    if (verified) "converged" else "approximate"
  } else {
    "failed"
  }
  result$usable <- verified
  result$reason <- if (verified) {
    "independently_verified"
  } else if (!integration_pass) {
    "quadrature_message_not_ok"
  } else if (!independent_pass) {
    "independent_verification_disagreement"
  } else if (!within_certified_bounds) {
    "certified_bounds_violation"
  } else if (!error_pass) {
    "numerical_error_budget_exceeded"
  } else if (!representation_pass) {
    "probability_scale_underflow"
  } else if (!interval_contains_estimate) {
    "direct_interval_not_nested_in_certified_bounds"
  } else {
    "nonfinite_or_out_of_probability_range"
  }
  numerical$verification_epsilon <- verification_epsilon
  numerical$within_certified_bounds <- within_certified_bounds
  numerical$direct_interval_contains_estimate <- interval_contains_estimate
  numerical$error_budget <- 1e-8
  numerical$error_budget_passed <- error_pass
  numerical$probability_representation_passed <- representation_pass
  numerical$candidate <- primary$value
  numerical$candidate_abs_error_bound <- abs_error_bound
  numerical$candidate_lower <- candidate_lower
  numerical$candidate_upper <- candidate_upper
  numerical$candidate_published <- verified
  result$numerical <- numerical
  class(result) <- "wmax_tail_result"
  result
}


.wmax_simulate_exact <- function(alpha, max_sticks) {
  n <- length(alpha)
  remainder <- rep.int(1, n)
  maximum <- numeric(n)
  active <- seq_len(n)
  iterations <- 0L

  while (length(active) > 0L && iterations < max_sticks) {
    iterations <- iterations + 1L
    breaks <- stats::rbeta(
      length(active), shape1 = 1, shape2 = alpha[active]
    )
    if (any(!is.finite(breaks)) || any(breaks < 0) || any(breaks > 1)) {
      return(list(
        maximum = maximum,
        unresolved = active,
        max_iterations = iterations,
        remainder = remainder,
        numerical_failure = TRUE
      ))
    }
    weights <- remainder[active] * breaks
    maximum[active] <- pmax(maximum[active], weights)
    remainder[active] <- remainder[active] * (1 - breaks)
    active <- which(remainder > maximum)
  }

  list(
    maximum = maximum,
    unresolved = active,
    max_iterations = iterations,
    remainder = remainder,
    numerical_failure = FALSE
  )
}


.wmax_wilson_interval <- function(successes, n, conf_level) {
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  estimate <- successes / n
  denominator <- 1 + z^2 / n
  center <- (estimate + z^2 / (2 * n)) / denominator
  half_width <- z * sqrt(
    estimate * (1 - estimate) / n + z^2 / (4 * n^2)
  ) / denominator
  c(lower = max(0, center - half_width),
    upper = min(1, center + half_width))
}


.wmax_monte_carlo_result <- function(bounds, specification, n, seed,
                                     conf_level, max_sticks,
                                     warn_low_successes) {
  result <- .wmax_result_template(bounds)
  result$method <- "seeded_gem_monte_carlo_exact_stopping"
  result$status <- "approximate"
  result$reason <- "monte_carlo_sampling_error"

  if (is.null(seed)) {
    .dpprior_abort_invalid(
      "seed is required for the Monte Carlo W_max method",
      "dpprior_weight_mc_error", "seed", seed,
      "one non-negative integer", "missing_seed"
    )
  }
  n <- .dpprior_validate_count(
    n, "n", minimum = 1L,
    .subclass = "dpprior_weight_mc_error"
  )
  seed <- .dpprior_validate_count(
    seed, "seed", minimum = 0L,
    .subclass = "dpprior_weight_mc_error"
  )
  conf_level <- .dpprior_validate_probability(
    conf_level, "conf_level", scalar = TRUE, open = TRUE,
    .subclass = "dpprior_weight_mc_error"
  )
  max_sticks <- .dpprior_validate_count(
    max_sticks, "max_sticks", minimum = 1L,
    .subclass = "dpprior_weight_mc_error"
  )
  warn_low_successes <- .dpprior_validate_control(
    warn_low_successes, "warn_low_successes", type = "logical"
  )

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    previous_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", previous_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)

  alpha <- if (identical(specification$conditioning, "fixed_alpha")) {
    rep.int(specification$alpha, n)
  } else {
    stats::rgamma(n, shape = specification$a, rate = specification$b)
  }
  if (any(!is.finite(alpha)) || any(alpha <= 0)) {
    result$status <- "failed"
    result$usable <- FALSE
    result$reason <- "nonpositive_or_nonfinite_simulated_alpha"
    result$sampling <- list(
      seed = seed, n = n,
      invalid_alpha_draws = sum(!is.finite(alpha) | alpha <= 0)
    )
    class(result) <- "wmax_tail_result"
    return(result)
  }

  draws <- .wmax_simulate_exact(alpha, max_sticks)
  if (isTRUE(draws$numerical_failure)) {
    result$status <- "failed"
    result$usable <- FALSE
    result$reason <- "invalid_beta_random_draw"
    result$sampling <- list(
      seed = seed, n = n, unresolved = as.integer(length(draws$unresolved)),
      max_sticks = max_sticks,
      max_iterations = as.integer(draws$max_iterations),
      stopping_rule = "remainder_not_greater_than_current_maximum"
    )
    class(result) <- "wmax_tail_result"
    return(result)
  }
  unresolved <- length(draws$unresolved)
  if (unresolved > 0L) {
    result$status <- "failed"
    result$usable <- FALSE
    result$reason <- "exact_stopping_limit_reached"
    result$sampling <- list(
      seed = seed, n = n, unresolved = as.integer(unresolved),
      max_sticks = max_sticks,
      max_iterations = as.integer(draws$max_iterations),
      stopping_rule = "remainder_not_greater_than_current_maximum"
    )
    class(result) <- "wmax_tail_result"
    return(result)
  }

  successes <- sum(draws$maximum > bounds$threshold)
  estimate <- successes / n
  wilson <- .wmax_wilson_interval(successes, n, conf_level)
  hoeffding <- sqrt(log(2 / (1 - conf_level)) / (2 * n))
  raw_hoeffding <- c(
    lower = max(0, estimate - hoeffding),
    upper = min(1, estimate + hoeffding)
  )
  nested_hoeffding <- c(
    lower = max(bounds$lower_bound, raw_hoeffding[["lower"]]),
    upper = min(bounds$upper_bound, raw_hoeffding[["upper"]])
  )
  interval_nonempty <-
    nested_hoeffding[["lower"]] <= nested_hoeffding[["upper"]]
  interval_contains_estimate <- interval_nonempty &&
    nested_hoeffding[["lower"]] <= estimate &&
    estimate <= nested_hoeffding[["upper"]]

  if (interval_contains_estimate) {
    result$estimate <- estimate
    result$lower <- unname(nested_hoeffding[["lower"]])
    result$upper <- unname(nested_hoeffding[["upper"]])
    result$abs_error_bound <- hoeffding
  }
  result$verified <- interval_contains_estimate
  result$estimate_usable <- FALSE
  # Monte Carlo remains an approximate result.  A valid nested sampling
  # interval is recorded by `verified`, but the Phase 2 status contract keeps
  # approximate results unusable without a separate public opt-in.
  result$usable <- FALSE
  if (!interval_contains_estimate) {
    result$reason <- "sampling_estimate_outside_certified_interval"
  }
  result$sampling <- list(
    seed = seed,
    n = n,
    raw_estimate = estimate,
    successes = as.integer(successes),
    failures = as.integer(n - successes),
    conf_level = conf_level,
    wilson_interval = wilson,
    hoeffding_interval_unintersected = raw_hoeffding,
    hoeffding_interval_nested = nested_hoeffding,
    hoeffding_half_width = hoeffding,
    error_statement_type = "distribution_free_sampling_probability_bound",
    interval_intersection_nonempty = interval_nonempty,
    interval_contains_estimate = interval_contains_estimate,
    unresolved = 0L,
    max_sticks = max_sticks,
    max_iterations = as.integer(draws$max_iterations),
    stopping_rule = "remainder_not_greater_than_current_maximum",
    stick_truncation_error = 0,
    low_success_count = successes < 25L
  )
  result$provenance$random_generation <- "GEM(1, alpha) sticks"
  result$provenance$sampling_error <-
    "binomial Wilson interval and distribution-free Hoeffding bound"

  if (warn_low_successes && successes < 25L) {
    .dpprior_warn(
      sprintf(
        paste(
          "W_max Monte Carlo observed only %d exceedances in %d draws;",
          "relative precision is poor. Increase n or use the deterministic",
          "method when threshold >= 0.5"
        ),
        successes, n
      ),
      "dpprior_mc_precision_warning", "n", n,
      "enough draws for at least 25 exceedances", "low_success_count"
    )
  }

  class(result) <- "wmax_tail_result"
  result
}


#' Largest DP Population Weight Tail Probability
#'
#' Computes or bounds \eqn{P(W_{max}>threshold)} without conflating it with
#' the first size-biased tail \eqn{P(W_{SB}>threshold)}. For thresholds at
#' least 0.5, \code{method = "auto"} uses a stable one-dimensional integral
#' and an independent calculation. Below 0.5, deterministic direct evaluation
#' is intentionally deferred and \code{"auto"} returns certified bounds only.
#' Explicit seeded Monte Carlo is available at every threshold and uses an
#' exact GEM stopping rule, so its only approximation is sampling error.
#'
#' @param threshold One finite probability in the closed unit interval.
#' @param alpha Optional fixed positive DP concentration.
#' @param a,b Optional positive Gamma shape and rate for alpha.
#' @param method One of \code{"auto"}, \code{"deterministic"},
#'   \code{"bounds"}, or \code{"monte_carlo"}.
#' @param rel_tol,abs_tol Positive numerical tolerances. Deterministic public
#'   results require \code{rel_tol <= 1e-10} and \code{abs_tol <= 1e-12}.
#' @param subdivisions Positive integration subdivision limit.
#' @param n Monte Carlo draw count; the default supports a worst-case normal
#'   95-percent half-width of about 0.002.
#' @param seed Required non-negative integer for Monte Carlo.
#' @param conf_level Monte Carlo confidence level in the open unit interval.
#' @param max_sticks Maximum GEM iterations. Reaching it produces an explicit
#'   failed result; an unresolved remainder is never discarded.
#' @param warn_low_successes Whether to warn when fewer than 25 exceedances
#'   are observed.
#'
#' @return A \code{wmax_tail_result} with the named estimand, conditioning,
#'   direct estimate (when available), direct numerical/sampling interval,
#'   absolute-error statement, universal certified bounds, status,
#'   verification flag, method, and provenance. A bounds-only result has
#'   \code{estimate = NA} and \code{usable = FALSE}; neither certified endpoint
#'   is used as an estimate. Its certified interval remains available with
#'   \code{bounds_usable = TRUE}, or directly from \code{wmax_tail_bounds()}.
#'
#' @examples
#' prob_wmax_exceeds(0.5, alpha = 1) # log(2)
#' prob_wmax_exceeds(0.9, a = 2, b = 1)
#' prob_wmax_exceeds(0.4, alpha = 1) # certified bounds only
#'
#' @seealso \code{\link{wmax_tail_bounds}},
#'   \code{\link{prob_wsb_exceeds}}
#'
#' @family weights_w1
#' @export
prob_wmax_exceeds <- function(
    threshold, alpha = NULL, a = NULL, b = NULL,
    method = c("auto", "deterministic", "bounds", "monte_carlo"),
    rel_tol = 1e-10, abs_tol = 1e-12, subdivisions = 1000L,
    n = 250000L, seed = NULL, conf_level = 0.95,
    max_sticks = 100000L, warn_low_successes = TRUE) {
  if (length(method) > 1L && identical(
    method, c("auto", "deterministic", "bounds", "monte_carlo")
  )) {
    method <- "auto"
  }
  method <- .wmax_validate_method(method)
  specification <- .wmax_validate_specification(alpha, a, b)
  bounds <- wmax_tail_bounds(
    threshold,
    alpha = specification$alpha,
    a = specification$a,
    b = specification$b
  )

  if (bounds$threshold %in% c(0, 1)) {
    return(.wmax_boundary_result(bounds, requested_method = method))
  }
  if (identical(method, "bounds")) {
    return(.wmax_bounds_only_result(bounds, "bounds_requested"))
  }
  if (identical(method, "monte_carlo")) {
    return(.wmax_monte_carlo_result(
      bounds, specification, n, seed, conf_level, max_sticks,
      warn_low_successes
    ))
  }
  if (bounds$threshold < 0.5) {
    return(.wmax_bounds_only_result(
      bounds, "deterministic_method_unsupported_below_half"
    ))
  }

  rel_tol <- .dpprior_validate_scalar(
    rel_tol, "rel_tol", lower = 0, upper = 1e-10,
    lower_open = TRUE, .subclass = "dpprior_weight_control_error"
  )
  abs_tol <- .dpprior_validate_scalar(
    abs_tol, "abs_tol", lower = 0, upper = 1e-12,
    lower_open = TRUE, .subclass = "dpprior_weight_control_error"
  )
  subdivisions <- .dpprior_validate_count(
    subdivisions, "subdivisions", minimum = 100L,
    .subclass = "dpprior_weight_control_error"
  )

  .wmax_deterministic_result(
    bounds, specification, rel_tol, abs_tol, subdivisions
  )
}


#' @export
print.wmax_tail_result <- function(x, digits = 6, ...) {
  conditioning <- if (identical(x$conditioning, "fixed_alpha")) {
    sprintf("fixed alpha = %.*g", digits, x$params$alpha)
  } else {
    sprintf(
      "Gamma(shape = %.*g, rate = %.*g)",
      digits, x$params$a, digits, x$params$b
    )
  }
  cat("Largest DP population weight tail (W_max)\n")
  cat(sprintf("  Status: %s\n", x$status))
  cat(sprintf("  Verified: %s\n", if (isTRUE(x$verified)) "yes" else "no"))
  cat(sprintf("  Method used: %s\n", x$method))
  cat(sprintf("  Conditioning: %s\n", conditioning))
  cat(sprintf("  Threshold: %.*g\n", digits, x$threshold))
  if (isTRUE(x$verified) && is.finite(x$estimate)) {
    cat(sprintf("  Estimate: %.*g\n", digits, x$estimate))
    cat(sprintf(
      "  Reported interval: [%.*g, %.*g]\n",
      digits, x$lower, digits, x$upper
    ))
  } else {
    cat("  Estimate: unavailable\n")
    if (!is.null(x$reason) && length(x$reason) == 1L &&
        !is.na(x$reason)) {
      cat(sprintf("  Reason: %s\n", x$reason))
    }
  }
  cat(sprintf(
    "  Certified W_max bounds: [%.*g, %.*g]\n",
    digits, x$lower_bound, digits, x$upper_bound
  ))
  cat(sprintf(
    "  First size-biased tail P(W_SB > threshold): %.*g\n",
    digits, x$size_biased_tail
  ))
  invisible(x)
}

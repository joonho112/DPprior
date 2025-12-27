# =============================================================================
# Module 00: Constants and Utility Functions
# =============================================================================
#
# This module provides:
# 1. Global constants used throughout the DPprior package
# 2. Numerically stable log-sum-exp and softmax operations
# 3. Input validation helper functions
#
# Author: JoonHo Lee
# Date: December 2025
# Part of: DPprior R Package
# =============================================================================

# =============================================================================
# Global Constants
# =============================================================================

#' Euler-Mascheroni Constant
#' @description The Euler-Mascheroni constant (gamma), approximately 0.5772.
#'   Used in harmonic sum approximations and asymptotic expansions.
#' @keywords internal
.EULER_GAMMA <- 0.5772156649015329

#' Maximum Supported Sample Size
#' @description Default maximum value of J (sample size) supported by the package.
#'   Pre-computed Stirling number tables are limited to this size.
#' @keywords internal
.MAX_J_DEFAULT <- 500L

#' Default Number of Quadrature Nodes
#' @description Default number of Gauss-Laguerre quadrature nodes for
#'   numerical integration against Gamma distributions.
#' @keywords internal
.QUAD_NODES_DEFAULT <- 80L

#' Newton Method Convergence Tolerance
#' @description Convergence tolerance for Newton's method in moment matching.
#' @keywords internal
.TOL_NEWTON <- 1e-8

#' Feasibility Projection Buffer
#' @description Small buffer for projecting to the feasible region of the
#'   negative binomial approximation.
#' @keywords internal
.TOL_PROJECTION_BUFFER <- 1e-6


# =============================================================================
# Numerically Stable Operations
# =============================================================================

#' Numerically Stable Log-Sum-Exp (Binary)
#'
#' Computes \code{log(exp(a) + exp(b))} in a numerically stable way,
#' avoiding overflow and underflow.
#'
#' @param a Numeric vector of log-scale values.
#' @param b Numeric vector of log-scale values (recycled to match length of \code{a}).
#'
#' @return Numeric vector of \code{log(exp(a) + exp(b))}.
#'
#' @details
#' Uses the identity:
#' \deqn{\log(\exp(a) + \exp(b)) = \max(a,b) + \log(1 + \exp(-|a-b|))}
#'
#' This formulation ensures numerical stability even for extreme values
#' (e.g., \code{a = 1000} or \code{a = -1000}).
#'
#' Special cases:
#' \itemize{
#'   \item If both inputs are \code{-Inf}, returns \code{-Inf}.
#'   \item If either input is \code{Inf}, returns \code{Inf}.
#' }
#'
#' @examples
#' # Standard case
#' logsumexp(log(2), log(3))
#'
#' # Extreme values that would overflow with naive implementation
#' logsumexp(1000, 1000)
#'
#' # Edge cases with Inf
#' logsumexp(-Inf, -Inf)
#'
#' @seealso \code{\link{logsumexp_vec}} for vector input
#'
#' @export
logsumexp <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log1p(exp(-abs(a - b)))

  # Handle edge cases for Inf values
  both_pos_inf <- is.infinite(a) & (a > 0) & is.infinite(b) & (b > 0)
  both_neg_inf <- is.infinite(a) & (a < 0) & is.infinite(b) & (b < 0)
  out[both_pos_inf] <- Inf
  out[both_neg_inf] <- -Inf

  out
}


#' Vectorized Log-Sum-Exp
#'
#' Computes \code{log(sum(exp(x)))} for a numeric vector in a numerically
#' stable way.
#'
#' @param x Numeric vector of log-scale values.
#'
#' @return Scalar value equal to \code{log(sum(exp(x)))}.
#'
#' @details
#' Subtracts the maximum before exponentiating to prevent overflow:
#' \deqn{\log\sum_i \exp(x_i) = \max_i x_i + \log\sum_i \exp(x_i - \max_i x_i)}
#'
#' Special cases:
#' \itemize{
#'   \item If all entries are \code{-Inf}, returns \code{-Inf}.
#'   \item If any entry is \code{Inf}, returns \code{Inf}.
#'   \item Empty vector throws an error.
#' }
#'
#' @examples
#' # Sum of equal values
#' logsumexp_vec(c(0, 0, 0, 0))
#'
#' # Extreme values
#' logsumexp_vec(c(1000, 1000, 1000))
#'
#' # All -Inf
#' logsumexp_vec(c(-Inf, -Inf))
#'
#' @seealso \code{\link{logsumexp}} for binary operation
#'
#' @export
logsumexp_vec <- function(x) {
  if (length(x) == 0L) {
    stop("x must have positive length", call. = FALSE)
  }

  x_max <- max(x)

  # All -Inf -> log(0) = -Inf
  if (is.infinite(x_max) && x_max < 0 && all(is.infinite(x) & x < 0)) {
    return(-Inf)
  }

  # Any +Inf -> sum exp = Inf
  if (is.infinite(x_max) && x_max > 0) {
    return(Inf)
  }

  x_max + log(sum(exp(x - x_max)))
}


#' Numerically Stable Softmax
#'
#' Computes the softmax transformation of a numeric vector, returning
#' a probability vector that sums to 1.
#'
#' @param x Numeric vector of log-odds or arbitrary real values.
#'
#' @return Numeric vector of probabilities summing to 1.
#'
#' @details
#' The softmax function is defined as:
#' \deqn{p_i = \frac{\exp(x_i)}{\sum_j \exp(x_j)}}
#'
#' This implementation subtracts the maximum value before exponentiating
#' to ensure numerical stability for extreme inputs.
#'
#' Special cases:
#' \itemize{
#'   \item If \code{x} contains \code{Inf} values, the probability mass is
#'         split uniformly across all \code{Inf} entries.
#'   \item Empty vector throws an error.
#' }
#'
#' @examples
#' softmax(c(1, 2, 3))
#' sum(softmax(c(1, 2, 3)))
#'
#' # Works with extreme values
#' softmax(c(1000, 1001, 1002))
#'
#' # Inf handling
#' softmax(c(1, Inf, Inf))
#'
#' @export
softmax <- function(x) {
  if (length(x) == 0L) {
    stop("x must have positive length", call. = FALSE)
  }

  # If any +Inf is present, split mass uniformly across them
  pos_inf <- is.infinite(x) & x > 0
  if (any(pos_inf)) {
    out <- rep(0, length(x))
    out[pos_inf] <- 1 / sum(pos_inf)
    return(out)
  }

  x_max <- max(x)
  if (!is.finite(x_max)) {
    stop("softmax() requires at least one finite value", call. = FALSE)
  }

  exp_x <- exp(x - x_max)
  exp_x / sum(exp_x)
}


# =============================================================================
# Input Validation Helpers
# =============================================================================

#' Assert Positive Values
#'
#' Validates that all elements of a numeric vector are strictly positive
#' and finite. Throws an informative error if validation fails.
#'
#' @param x Numeric vector to validate.
#' @param name Character string naming the parameter (for error messages).
#'
#' @return Invisible \code{TRUE} if validation passes.
#'
#' @examples
#' \dontrun{
#' assert_positive(c(1, 2, 3), "alpha")
#' assert_positive(c(1, -1), "alpha")
#' }
#'
#' @keywords internal
assert_positive <- function(x, name = "x") {
  if (!is.numeric(x) || any(!is.finite(x)) || any(x <= 0)) {
    stop(sprintf("%s must be finite and positive", name), call. = FALSE)
  }
  invisible(TRUE)
}


#' Assert Valid Sample Size J
#'
#' Validates that J is a positive integer within the supported range.
#'
#' @param J Sample size to validate.
#'
#' @return Invisible \code{TRUE} if validation passes.
#'
#' @examples
#' \dontrun{
#' assert_valid_J(50)
#' assert_valid_J(0)
#' assert_valid_J(1000)
#' }
#'
#' @keywords internal
assert_valid_J <- function(J) {
  if (!is.numeric(J) || length(J) != 1L || !is.finite(J) ||
      J != floor(J) || J < 1L) {
    stop("J must be a positive integer", call. = FALSE)
  }
  if (J > .MAX_J_DEFAULT) {
    stop(sprintf("J exceeds maximum supported value (%d)", .MAX_J_DEFAULT),
         call. = FALSE)
  }
  invisible(TRUE)
}


#' Assert Valid Probability
#'
#' Validates that all elements of a numeric vector are valid probabilities
#' in the range \[0, 1\] and finite.
#'
#' @param p Numeric vector of probability values to validate.
#' @param name Character string naming the parameter (for error messages).
#'
#' @return Invisible \code{TRUE} if validation passes.
#'
#' @examples
#' \dontrun{
#' assert_probability(0.5, "p")
#' assert_probability(1.5, "p")
#' }
#'
#' @keywords internal
assert_probability <- function(p, name = "p") {
  if (!is.numeric(p) || any(!is.finite(p)) || any(p < 0) || any(p > 1)) {
    stop(sprintf("%s must be finite and in [0, 1]", name), call. = FALSE)
  }
  invisible(TRUE)
}


#' Assert Non-negative Values
#'
#' Validates that all elements of a numeric vector are non-negative and finite.
#'
#' @param x Numeric vector to validate.
#' @param name Character string naming the parameter (for error messages).
#'
#' @return Invisible \code{TRUE} if validation passes.
#'
#' @keywords internal
assert_nonnegative <- function(x, name = "x") {
  if (!is.numeric(x) || any(!is.finite(x)) || any(x < 0)) {
    stop(sprintf("%s must be finite and non-negative", name), call. = FALSE)
  }
  invisible(TRUE)
}


#' Assert Valid Cluster Count k
#'
#' Validates that k is a valid cluster count for given sample size J.
#'
#' @param k Cluster count to validate.
#' @param J Sample size (k must be in 1:J).
#'
#' @return Invisible \code{TRUE} if validation passes.
#'
#' @keywords internal
assert_valid_k <- function(k, J) {
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) ||
      k != floor(k) || k < 1L || k > J) {
    stop(sprintf("k must be an integer in [1, %d]", J), call. = FALSE)
  }
  invisible(TRUE)
}

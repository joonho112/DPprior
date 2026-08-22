# =============================================================================
# Module 00: Constants and Utility Functions
# =============================================================================
#
# This module provides:
# 1. Global constants used throughout the DPprior package
# 2. Numerically stable log-sum-exp and softmax operations
# 3. Input validation helper functions
#
# Author: JoonHo Lee (jlee296@ua.edu)
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

#' Probability-Mass Normalization Tolerance
#' @description Absolute tolerance used by the shared PMF validator when
#'   checking that probability mass sums to one.
#' @keywords internal
.TOL_PMF_SUM <- 1e-10

#' Feasibility Projection Buffer
#' @description Small buffer for projecting to the feasible region of the
#'   negative binomial approximation.
#' @keywords internal
.TOL_PROJECTION_BUFFER <- 1e-6

# --- Numerical safety thresholds ---
.ALPHA_SMALL <- 1e-10          # Floor for alpha in safe moment computations
.TOL_SINGULAR_JACOBIAN <- 1e-12  # Jacobian determinant singularity threshold
.PENALTY_INF <- 1e10            # Penalty return for infeasible parameters
.GD_FALLBACK_STEP <- 0.1        # Gradient descent fallback step size
.LOG_BOUNDS_DEFAULT <- c(-15, 15)  # L-BFGS-B bounds on log(a), log(b)
.EXP_MAX <- 700                 # Maximum safe exponent for exp()


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
#' \dontrun{
#' softmax(c(1, 2, 3))
#' sum(softmax(c(1, 2, 3)))
#'
#' # Works with extreme values
#' softmax(c(1000, 1001, 1002))
#'
#' # Inf handling
#' softmax(c(1, Inf, Inf))
#'
#' }
#' @keywords internal
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

# Construct a dependency-free package condition. More specific classes precede
# the shared dpprior_condition root so callers can catch either level.
.dpprior_new_condition <- function(message, classes, call = NULL, ...) {
  structure(
    c(list(message = as.character(message), call = call), list(...)),
    class = unique(c(classes, "dpprior_condition", "condition"))
  )
}


# Signal invalid user input with stable classes and machine-readable fields.
.dpprior_abort_invalid <- function(message,
                                   subclass = character(),
                                   argument = NULL,
                                   value = NULL,
                                   expected = NULL,
                                   code = "invalid",
                                   call = NULL) {
  stop(.dpprior_new_condition(
    message = message,
    classes = c(subclass, "dpprior_invalid_input", "dpprior_error", "error"),
    call = call,
    argument = argument,
    value = value,
    expected = expected,
    code = code
  ))
}


# Signal a package warning with the same field contract as package errors.
.dpprior_warn <- function(message,
                          subclass = character(),
                          argument = NULL,
                          value = NULL,
                          expected = NULL,
                          code = "warning",
                          call = NULL) {
  warning(.dpprior_new_condition(
    message = message,
    classes = c(subclass, "dpprior_warning", "warning"),
    call = call,
    argument = argument,
    value = value,
    expected = expected,
    code = code
  ))
  invisible(NULL)
}


.dpprior_is_plain_numeric <- function(x) {
  is.numeric(x) && is.null(dim(x)) && !is.object(x)
}


.dpprior_validate_numeric_input <- function(x, name, scalar, subclass) {
  if (!is.numeric(x)) {
    .dpprior_abort_invalid(
      sprintf("%s must be numeric", name),
      c(subclass, "dpprior_type_error"), name, x, "numeric", "type"
    )
  }
  if (!is.null(dim(x)) || is.object(x)) {
    .dpprior_abort_invalid(
      sprintf(
        "%s must be an ordinary numeric vector without dimensions or a custom class",
        name
      ),
      c(subclass, "dpprior_type_error"), name, x,
      "ordinary numeric vector without dimensions or a custom class",
      "type"
    )
  }
  if ((scalar && length(x) != 1L) || (!scalar && length(x) == 0L)) {
    expected <- if (scalar) "length 1" else "positive length"
    .dpprior_abort_invalid(
      sprintf("%s must have %s", name, expected),
      c(subclass, "dpprior_length_error"), name, x, expected, "length"
    )
  }
  if (anyNA(x)) {
    .dpprior_abort_invalid(
      sprintf("%s must not contain missing values", name),
      c(subclass, "dpprior_missing_error"),
      name, x, "non-missing", "missing"
    )
  }
  if (any(!is.finite(x))) {
    .dpprior_abort_invalid(
      sprintf("%s must contain only finite values", name),
      c(subclass, "dpprior_nonfinite_error"),
      name, x, "finite", "nonfinite"
    )
  }
  x
}


.dpprior_format_bound <- function(x) {
  format(x, digits = 15L, trim = TRUE)
}


.dpprior_validate_bounds <- function(x, name, lower, upper,
                                     lower_open, upper_open, subclass) {
  below <- !is.null(lower) && if (lower_open) any(x <= lower) else any(x < lower)
  above <- !is.null(upper) && if (upper_open) any(x >= upper) else any(x > upper)
  if (!below && !above) {
    return(x)
  }

  if (!is.null(lower) && !is.null(upper)) {
    expected <- sprintf(
      "%s%s, %s%s",
      if (lower_open) "(" else "[", .dpprior_format_bound(lower),
      .dpprior_format_bound(upper), if (upper_open) ")" else "]"
    )
    message <- sprintf("%s must be in %s", name, expected)
  } else if (!is.null(lower)) {
    expected <- sprintf(
      "%s %s", if (lower_open) ">" else ">=", .dpprior_format_bound(lower)
    )
    message <- sprintf("%s must be %s", name, expected)
  } else {
    expected <- sprintf(
      "%s %s", if (upper_open) "<" else "<=", .dpprior_format_bound(upper)
    )
    message <- sprintf("%s must be %s", name, expected)
  }

  .dpprior_abort_invalid(
    message, c(subclass, "dpprior_bounds_error"),
    name, x, expected, "bounds"
  )
}


# Validate a finite numeric scalar and optional open or closed bounds.
.dpprior_validate_scalar <- function(x, name = "x",
                                     lower = NULL, upper = NULL,
                                     lower_open = FALSE,
                                     upper_open = FALSE,
                                     .subclass = "dpprior_scalar_error") {
  x <- .dpprior_validate_numeric_input(x, name, TRUE, .subclass)
  .dpprior_validate_bounds(
    x, name, lower, upper, lower_open, upper_open, .subclass
  )
}


# Validate one probability, or a non-empty vector of probabilities.
.dpprior_validate_probability <- function(x, name = "p", scalar = TRUE,
                                          open = FALSE,
                                          .subclass = "dpprior_probability_error") {
  open <- .dpprior_validate_control(
    open, paste0(name, "_open"), type = "logical"
  )
  x <- .dpprior_validate_numeric_input(x, name, scalar, .subclass)
  .dpprior_validate_bounds(x, name, 0, 1, open, open, .subclass)
}


# Validate a non-negative integer-valued scalar, optionally within bounds.
.dpprior_validate_count <- function(x, name = "n", minimum = 0L,
                                    maximum = .Machine$integer.max,
                                    minimum_open = FALSE,
                                    maximum_open = FALSE,
                                    .subclass = "dpprior_count_error") {
  x <- .dpprior_validate_numeric_input(x, name, TRUE, .subclass)
  if (x != floor(x)) {
    .dpprior_abort_invalid(
      sprintf("%s must be an integer", name),
      c(.subclass, "dpprior_integer_error"),
      name, x, "integer-valued", "integer"
    )
  }
  x <- .dpprior_validate_bounds(
    x, name, minimum, maximum, minimum_open, maximum_open, .subclass
  )
  as.integer(x)
}


# Validate a scalar numerical or logical control value.
.dpprior_validate_control <- function(x, name = "control",
                                      type = c("numeric", "count",
                                               "probability", "logical"),
                                      lower = NULL, upper = NULL,
                                      lower_open = FALSE,
                                      upper_open = FALSE) {
  type <- match.arg(type)
  subclass <- "dpprior_control_error"
  if (type == "numeric") {
    return(.dpprior_validate_scalar(
      x, name, lower, upper, lower_open, upper_open, subclass
    ))
  }
  if (type == "count") {
    return(.dpprior_validate_count(
      x, name,
      if (is.null(lower)) 0L else lower,
      if (is.null(upper)) .Machine$integer.max else upper,
      lower_open, upper_open,
      subclass
    ))
  }
  if (type == "probability") {
    x <- .dpprior_validate_numeric_input(x, name, TRUE, subclass)
    return(.dpprior_validate_bounds(
      x, name,
      if (is.null(lower)) 0 else lower,
      if (is.null(upper)) 1 else upper,
      lower_open, upper_open, subclass
    ))
  }

  if (!is.logical(x) || !is.null(dim(x)) || is.object(x)) {
    .dpprior_abort_invalid(
      sprintf("%s must be logical", name),
      c(subclass, "dpprior_type_error"), name, x, "logical", "type"
    )
  }
  if (length(x) != 1L) {
    .dpprior_abort_invalid(
      sprintf("%s must have length 1", name),
      c(subclass, "dpprior_length_error"), name, x, "length 1", "length"
    )
  }
  if (is.na(x)) {
    .dpprior_abort_invalid(
      sprintf("%s must not contain missing values", name),
      c(subclass, "dpprior_missing_error"),
      name, x, "non-missing", "missing"
    )
  }
  x
}


# Re-signal legacy assertion failures with their historical message while
# attaching the same reason-specific classes used by the shared validators.
.dpprior_abort_legacy_numeric <- function(message, x, name, expected,
                                           subclass, scalar = FALSE,
                                           integer = FALSE) {
  if (!.dpprior_is_plain_numeric(x)) {
    reason <- "dpprior_type_error"
    code <- "type"
  } else if ((scalar && length(x) != 1L) || (!scalar && length(x) == 0L)) {
    reason <- "dpprior_length_error"
    code <- "length"
  } else if (anyNA(x)) {
    reason <- "dpprior_missing_error"
    code <- "missing"
  } else if (any(!is.finite(x))) {
    reason <- "dpprior_nonfinite_error"
    code <- "nonfinite"
  } else if (integer && any(x != floor(x))) {
    reason <- "dpprior_integer_error"
    code <- "integer"
  } else {
    reason <- "dpprior_bounds_error"
    code <- "bounds"
  }

  .dpprior_abort_invalid(
    message, c(subclass, reason), name, x, expected, code
  )
}


# Reject array/matrix and classed atomic inputs where an ordinary numeric vector
# is part of the mathematical contract. Names are harmless and are retained;
# dimensional and class attributes can otherwise be silently discarded by
# as.numeric(), changing the interpretation of a probability input.
.dpprior_validate_plain_vector <- function(x, name,
                                           subclass = "dpprior_invalid_input") {
  if (!is.null(dim(x)) || is.object(x)) {
    .dpprior_abort_invalid(
      sprintf("%s must be an ordinary vector without dimensions or a custom class.",
              name),
      c(subclass, "dpprior_type_error"),
      name, x, "ordinary vector without dimensions or a custom class",
      "vector_required"
    )
  }
  x
}


# Validate a finite, non-negative PMF without silently normalizing it.
.dpprior_validate_pmf <- function(x, name = "pmf", expected_length = NULL,
                                  require_sum_one = TRUE,
                                  tolerance = .TOL_PMF_SUM) {
  tolerance <- .dpprior_validate_control(
    tolerance, "tolerance", "numeric", lower = 0, lower_open = TRUE
  )
  require_sum_one <- .dpprior_validate_control(
    require_sum_one, "require_sum_one", "logical"
  )
  x <- .dpprior_validate_plain_vector(x, name, "dpprior_pmf_error")
  x <- .dpprior_validate_numeric_input(x, name, FALSE, "dpprior_pmf_error")

  if (!is.null(expected_length)) {
    expected_length <- .dpprior_validate_count(
      expected_length, "expected_length", 1L,
      .subclass = "dpprior_control_error"
    )
    if (length(x) != expected_length) {
      .dpprior_abort_invalid(
        sprintf("%s must have length %d", name, expected_length),
        c("dpprior_pmf_error", "dpprior_length_error"),
        name, x, sprintf("length %d", expected_length), "length"
      )
    }
  }

  x <- .dpprior_validate_bounds(
    x, name, 0, NULL, FALSE, FALSE, "dpprior_pmf_error"
  )
  total <- sum(x)
  if (!is.finite(total)) {
    .dpprior_abort_invalid(
      sprintf("%s must have a finite total mass", name),
      c("dpprior_pmf_error", "dpprior_nonfinite_error"),
      name, x, "finite total mass", "sum_nonfinite"
    )
  }
  if (total <= 0) {
    .dpprior_abort_invalid(
      sprintf("%s must have positive total mass", name),
      c("dpprior_pmf_error", "dpprior_pmf_mass_error"),
      name, x, "positive total mass", "zero_mass"
    )
  }
  if (require_sum_one && abs(total - 1) > tolerance) {
    .dpprior_abort_invalid(
      sprintf(
        "%s must sum to 1 within tolerance %.3g (sum = %.17g)",
        name, tolerance, total
      ),
      c("dpprior_pmf_error", "dpprior_normalization_error"),
      name, x, sprintf("sum to 1 within %.3g", tolerance), "normalization"
    )
  }
  x
}

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
  if (!.dpprior_is_plain_numeric(x) || length(x) == 0L ||
      any(!is.finite(x)) || any(x <= 0)) {
    .dpprior_abort_legacy_numeric(
      sprintf("%s must be finite and positive", name), x, name,
      "finite and positive", "dpprior_numeric_error"
    )
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
  if (!.dpprior_is_plain_numeric(J) || length(J) != 1L || !is.finite(J) ||
      J != floor(J) || J < 1L) {
    .dpprior_abort_legacy_numeric(
      "J must be a positive integer", J, "J", "positive integer",
      "dpprior_count_error", scalar = TRUE, integer = TRUE
    )
  }
  if (J > .MAX_J_DEFAULT) {
    .dpprior_abort_invalid(
      sprintf("J exceeds maximum supported value (%d)", .MAX_J_DEFAULT),
      c("dpprior_count_error", "dpprior_bounds_error"), "J", J,
      sprintf("<= %d", .MAX_J_DEFAULT), "bounds"
    )
  }
  invisible(TRUE)
}


# Internal scalar integer validator used for numerical control arguments.
.as_integer_scalar <- function(x, name, min = NULL, max = NULL) {
  if (!.dpprior_is_plain_numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != floor(x)) {
    .dpprior_abort_legacy_numeric(
      sprintf("%s must be an integer", name), x, name, "integer scalar",
      "dpprior_count_error", scalar = TRUE, integer = TRUE
    )
  }
  if (!is.null(min) && x < min) {
    .dpprior_abort_invalid(
      sprintf("%s must be an integer >= %d", name, as.integer(min)),
      c("dpprior_count_error", "dpprior_bounds_error"), name, x,
      sprintf(">= %d", as.integer(min)), "bounds"
    )
  }
  if (!is.null(max) && x > max) {
    .dpprior_abort_invalid(
      sprintf("%s must be an integer <= %d", name, as.integer(max)),
      c("dpprior_count_error", "dpprior_bounds_error"), name, x,
      sprintf("<= %d", as.integer(max)), "bounds"
    )
  }
  as.integer(x)
}


# Maximum variance for a random variable supported on {1, ..., J} with fixed mean.
.max_var_K_fixed_mean <- function(J, mu_K) {
  (mu_K - 1) * (J - mu_K)
}


.assert_feasible_K_moments <- function(J, mu_K, var_K, tol = 1e-12) {
  var_upper <- .max_var_K_fixed_mean(J, mu_K)
  if (var_K > var_upper + tol) {
    .dpprior_abort_invalid(
      sprintf(
        "var_K = %.4g exceeds maximum possible variance %.4g for K in {1,...,%d} with mu_K = %.4g",
        var_K, var_upper, as.integer(J), mu_K
      ),
      c("dpprior_moment_feasibility_error", "dpprior_bounds_error"),
      "var_K", var_K,
      sprintf("<= %.17g", var_upper), "bounds"
    )
  }
  invisible(var_upper)
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
  if (!.dpprior_is_plain_numeric(p) || length(p) == 0L ||
      any(!is.finite(p)) || any(p < 0) || any(p > 1)) {
    .dpprior_abort_legacy_numeric(
      sprintf("%s must be finite and in [0, 1]", name), p, name,
      "finite and in [0, 1]", "dpprior_probability_error"
    )
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
  if (!.dpprior_is_plain_numeric(x) || length(x) == 0L ||
      any(!is.finite(x)) || any(x < 0)) {
    .dpprior_abort_legacy_numeric(
      sprintf("%s must be finite and non-negative", name), x, name,
      "finite and non-negative", "dpprior_numeric_error"
    )
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
  if (!.dpprior_is_plain_numeric(k) || length(k) != 1L || !is.finite(k) ||
      k != floor(k) || k < 1L || k > J) {
    .dpprior_abort_legacy_numeric(
      sprintf("k must be an integer in [1, %d]", J), k, "k",
      sprintf("integer in [1, %d]", J), "dpprior_count_error",
      scalar = TRUE, integer = TRUE
    )
  }
  invisible(TRUE)
}


# Null-coalescing operator: returns x if not NULL, otherwise y.
# Defined once here; do NOT duplicate in other modules.
# @noRd
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

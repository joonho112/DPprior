# =============================================================================
# Module 07: Score-Based Jacobian for Newton Solver
# =============================================================================
#
# This module provides functions for computing the Jacobian matrix of the
# marginal moment map F(a,b) = (M_1(a,b), V(a,b)) using score function
# identities, enabling Newton updates for moment matching.
#
# Numerical safeguards include stable finite-sum conditional moments, an exact
# zero-limit score control variate, explicit order-refinement provenance, and
# scale-aware conditioning diagnostics.
#
# Theory Reference: Lee (2026), Section 3.2
#
# Author: JoonHo Lee (jlee296@ua.edu)
# Date: December 2025
# Part of: DPprior R Package
# =============================================================================

# =============================================================================
# Constants for Numerical Stability
# =============================================================================

# .ALPHA_SMALL is now defined in R/00_constants.R (loaded first alphabetically)

#' Default Quadrature Nodes for Verification
#' @description Higher M value recommended for Jacobian verification due to
#'   slower convergence of score-weighted integrands.
#' @keywords internal
.QUAD_NODES_VERIFICATION <- 200L

# The shape-score integrand contains log(alpha).  A zero-limit control variate
# removes the spurious constant-times-score term, and a distinct higher order
# is used to quantify the remaining quadrature error.  These tolerances govern
# the refinement diagnostic and match the independent release certificate.
# Using the same mixed budget prevents an order-refinement result from being
# labelled converged under a looser rule than the release oracle.
.JACOBIAN_SCORE_MIN_NODES <- 80L
.JACOBIAN_DERIV_ABS_TOL <- 5e-5
.JACOBIAN_DERIV_REL_TOL <- 2e-3
.JACOBIAN_RCOND_SINGULAR <- 1e-12
.JACOBIAN_RCOND_ILL <- sqrt(.Machine$double.eps)


# =============================================================================
# Score Functions for Gamma(a, b)
# =============================================================================

#' Score Function with Respect to Shape Parameter a
#'
#' Computes the score function \eqn{s_a(\alpha) = \partial/\partial a \log g_{a,b}(\alpha)}
#' for the Gamma(a, b) distribution.
#'
#' @param alpha Numeric vector; points at which to evaluate.
#' @param a Numeric scalar; shape parameter of the Gamma distribution (> 0).
#' @param b Numeric scalar; rate parameter of the Gamma distribution (> 0).
#'
#' @return Numeric vector of the same length as \code{alpha}.
#'
#' @details
#' For the Gamma(shape = a, rate = b) distribution with density
#' \deqn{g_{a,b}(\alpha) = \frac{b^a}{\Gamma(a)} \alpha^{a-1} e^{-b\alpha},}
#' the score function with respect to \code{a} is:
#' \deqn{s_a(\alpha) = \log b - \psi(a) + \log \alpha,}
#' where \eqn{\psi} is the digamma function.
#'
#' A fundamental property of score functions is that their expectation is zero:
#' \deqn{E_{\alpha \sim g_{a,b}}[s_a(\alpha)] = 0.}
#'
#' \strong{Numerical Note:} The \code{log(alpha)} term causes slow quadrature
#' convergence for the raw score expectation, especially for small shape.
#' Production derivatives use a zero-limit control variate and a separately
#' recorded order-refinement check rather than treating a larger order as exact.
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @seealso \code{\link{score_b}} for the score with respect to b
#'
#' @examples
#' # Evaluate score at several points
#' alpha_vals <- c(0.5, 1.0, 2.0, 5.0)
#' score_a(alpha_vals, a = 2.0, b = 1.0)
#'
#' @export
score_a <- function(alpha, a, b) {
  assert_positive(alpha, "alpha")
  a <- .dpprior_validate_scalar(
    a, "a", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  b <- .dpprior_validate_scalar(
    b, "b", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )

  score <- suppressWarnings(log(b) - digamma(a) + log(alpha))
  if (any(!is.finite(score))) {
    .dpprior_abort_invalid(
      "shape score is non-finite for the supplied parameters",
      c("dpprior_score_error", "dpprior_numerical_error"),
      "score_a", score, "finite score values", "nonfinite_score"
    )
  }
  score
}


#' Score Function with Respect to Rate Parameter b
#'
#' Computes the score function \eqn{s_b(\alpha) = \partial/\partial b \log g_{a,b}(\alpha)}
#' for the Gamma(a, b) distribution.
#'
#' @param alpha Numeric vector; points at which to evaluate.
#' @param a Numeric scalar; shape parameter of the Gamma distribution (> 0).
#' @param b Numeric scalar; rate parameter of the Gamma distribution (> 0).
#'
#' @return Numeric vector of the same length as \code{alpha}.
#'
#' @details
#' For the Gamma(shape = a, rate = b) distribution, the score function with
#' respect to \code{b} is:
#' \deqn{s_b(\alpha) = \frac{a}{b} - \alpha.}
#'
#' A fundamental property of score functions is that their expectation is zero:
#' \deqn{E_{\alpha \sim g_{a,b}}[s_b(\alpha)] = 0.}
#'
#' Unlike \code{score_a}, this function is linear in \eqn{\alpha}, so its
#' expectation converges very quickly with quadrature.
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @seealso \code{\link{score_a}} for the score with respect to a
#'
#' @examples
#' # Evaluate score at several points
#' alpha_vals <- c(0.5, 1.0, 2.0, 5.0)
#' score_b(alpha_vals, a = 2.0, b = 1.0)
#'
#' @export
score_b <- function(alpha, a, b) {
  assert_positive(alpha, "alpha")
  a <- .dpprior_validate_scalar(
    a, "a", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  b <- .dpprior_validate_scalar(
    b, "b", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )

  score <- a / b - alpha
  if (any(!is.finite(score))) {
    .dpprior_abort_invalid(
      "rate score is non-finite for the supplied parameters",
      c("dpprior_score_error", "dpprior_numerical_error"),
      "score_b", score, "finite score values", "nonfinite_score"
    )
  }
  score
}


# =============================================================================
# Enhanced Conditional Moments with Small-Alpha Handling
# =============================================================================

#' Conditional Mean of K Given Alpha (Enhanced)
#'
#' Compatibility wrapper for the stable conditional mean implementation.
#'
#' @param J Integer; sample size.
#' @param alpha Numeric vector; concentration parameter values.
#'
#' @return Numeric vector of conditional means.
#'
#' @details
#' Delegates to \code{mean_K_given_alpha()}, whose finite-sum implementation is
#' stable without a discontinuous small-alpha substitution.
#'
#' @keywords internal
mean_K_given_alpha_safe <- function(J, alpha) {
  mean_K_given_alpha(J, alpha)
}


#' Conditional Variance of K Given Alpha (Enhanced)
#'
#' Compatibility wrapper for the stable conditional variance implementation.
#'
#' @param J Integer; sample size.
#' @param alpha Numeric vector; concentration parameter values.
#'
#' @return Numeric vector of conditional variances.
#'
#' @details
#' Delegates to \code{var_K_given_alpha()}, whose positive-term implementation
#' is stable without clipping or a discontinuous small-alpha substitution.
#'
#' @keywords internal
var_K_given_alpha_safe <- function(J, alpha) {
  var_K_given_alpha(J, alpha)
}


# =============================================================================
# Combined Moments and Jacobian Computation
# =============================================================================


# Evaluate marginal moments and the score-identity Jacobian at one quadrature
# order.  Subtracting the exact alpha -> 0 limits mu(0)=1 and
# E[K^2|alpha=0]=1 is an exact control variate because every Gamma score has
# expectation zero.  It removes the constant-times-log(alpha) term that makes
# the uncentered shape-score quadrature inaccurate, and makes J=1 exact.
.jacobian_components_at_order <- function(J, a, b, M) {
  quad <- build_gamma_quadrature(a, b, M)
  alphas <- quad$alpha_nodes
  weights <- quad$weights_normalized

  mu <- mean_K_given_alpha(J, alphas)
  conditional_var <- var_K_given_alpha(J, alphas)
  conditional_second <- conditional_var + mu^2

  marginal_mean <- .quadrature_weighted_sum(weights, mu)
  marginal_second <- .quadrature_weighted_sum(weights, conditional_second)
  within_alpha <- .quadrature_weighted_sum(weights, conditional_var)
  between_alpha <- .quadrature_weighted_sum(
    weights, (mu - marginal_mean)^2
  )
  marginal_var <- within_alpha + between_alpha
  if (!is.finite(marginal_var) || marginal_var < 0 ||
      !is.finite(within_alpha) || within_alpha < 0 ||
      !is.finite(between_alpha) || between_alpha < 0) {
    .dpprior_abort_invalid(
      "marginal variance decomposition is non-finite or negative",
      c("dpprior_jacobian_error", "dpprior_numerical_error"),
      "marginal_var",
      c(
        total = marginal_var,
        within_alpha = within_alpha,
        between_alpha = between_alpha
      ),
      "finite non-negative variance components", "variance_decomposition"
    )
  }
  if (J == 1L) {
    marginal_mean <- 1
    marginal_var <- 0
  }

  shape_score <- score_a(alphas, a, b)
  rate_score <- score_b(alphas, a, b)
  centered_mean <- mu - 1
  centered_second <- conditional_second - 1

  dmean_da <- .quadrature_weighted_sum(
    weights, centered_mean * shape_score
  )
  dmean_db <- .quadrature_weighted_sum(
    weights, centered_mean * rate_score
  )
  dsecond_da <- .quadrature_weighted_sum(
    weights, centered_second * shape_score
  )
  dsecond_db <- .quadrature_weighted_sum(
    weights, centered_second * rate_score
  )

  jacobian <- matrix(
    c(
      dmean_da, dmean_db,
      dsecond_da - 2 * marginal_mean * dmean_da,
      dsecond_db - 2 * marginal_mean * dmean_db
    ),
    nrow = 2L, ncol = 2L, byrow = TRUE,
    dimnames = list(c("dM1", "dV"), c("da", "db"))
  )

  list(
    mean = marginal_mean,
    var = marginal_var,
    jacobian = jacobian,
    score_expectation = c(
      da = .quadrature_weighted_sum(weights, shape_score),
      db = .quadrature_weighted_sum(weights, rate_score)
    ),
    quadrature_metadata = .quadrature_metadata(quad)
  )
}


# Return a scale-aware singularity diagnostic based on singular values.  A
# determinant alone is not suitable because it changes with parameter units.
.jacobian_matrix_condition <- function(jacobian, parameterization) {
  if (!is.matrix(jacobian) || !is.numeric(jacobian) ||
      !identical(dim(jacobian), c(2L, 2L)) ||
      any(!is.finite(jacobian))) {
    .dpprior_abort_invalid(
      "jacobian must be a finite 2 by 2 numeric matrix",
      c("dpprior_jacobian_error", "dpprior_numerical_error"),
      "jacobian", jacobian, "finite 2 by 2 numeric matrix",
      "invalid_jacobian"
    )
  }

  singular_values <- tryCatch(
    svd(jacobian, nu = 0L, nv = 0L)$d,
    error = function(error) {
      .dpprior_abort_invalid(
        sprintf(
          "singular-value decomposition of jacobian failed: %s",
          conditionMessage(error)
        ),
        c("dpprior_jacobian_error", "dpprior_numerical_error"),
        "jacobian", jacobian, "finite decomposable 2 by 2 matrix",
        "svd_failure"
      )
    }
  )
  if (any(!is.finite(singular_values))) {
    .dpprior_abort_invalid(
      "singular-value decomposition returned non-finite values",
      c("dpprior_jacobian_error", "dpprior_numerical_error"),
      "jacobian", jacobian, "finite singular values", "nonfinite_svd"
    )
  }
  largest <- max(singular_values)
  reciprocal_condition <- if (largest == 0) {
    0
  } else {
    min(singular_values) / largest
  }
  condition_number <- if (reciprocal_condition == 0) {
    Inf
  } else {
    1 / reciprocal_condition
  }
  rank_tolerance <- max(
    max(dim(jacobian)) * .Machine$double.eps,
    .JACOBIAN_RCOND_SINGULAR
  ) * largest
  numerical_rank <- if (largest == 0) {
    0L
  } else {
    as.integer(sum(singular_values > rank_tolerance))
  }
  status <- if (numerical_rank < 2L ||
                reciprocal_condition <= .JACOBIAN_RCOND_SINGULAR) {
    "singular"
  } else if (reciprocal_condition <= .JACOBIAN_RCOND_ILL) {
    "ill_conditioned"
  } else {
    "well_conditioned"
  }

  list(
    parameterization = parameterization,
    status = status,
    determinant = as.numeric(det(jacobian)),
    reciprocal_condition = reciprocal_condition,
    condition_number = condition_number,
    rank = numerical_rank,
    singular_values = singular_values,
    thresholds = c(
      singular = .JACOBIAN_RCOND_SINGULAR,
      ill_conditioned = .JACOBIAN_RCOND_ILL
    )
  )
}


.jacobian_conditioning <- function(jacobian, a, b) {
  raw <- .jacobian_matrix_condition(jacobian, "shape-rate")
  log_jacobian <- jacobian %*% diag(c(a, b))
  dimnames(log_jacobian) <- dimnames(jacobian)
  log_scale <- .jacobian_matrix_condition(
    log_jacobian, "log-shape-log-rate"
  )

  list(
    schema_version = 1L,
    status = log_scale$status,
    solver_parameterization = "log-shape-log-rate",
    raw = raw,
    solver = log_scale
  )
}


# Evaluate the score Jacobian at a release floor and, when possible, at a
# distinct higher order.  The higher-order result is returned; the difference
# is retained rather than silently assuming that refinement is monotone.
.jacobian_with_refinement <- function(J, a, b, M) {
  moments <- .jacobian_components_at_order(J, a, b, M)
  score_order <- max(as.integer(M), .JACOBIAN_SCORE_MIN_NODES)
  selected <- if (score_order == M) {
    moments
  } else {
    .jacobian_components_at_order(J, a, b, score_order)
  }
  verification_target <- if (score_order < .QUADRATURE_MAX_NODES) {
    max(2L * score_order, score_order + 40L)
  } else {
    NA_integer_
  }
  verification_order <- if (is.na(verification_target)) {
    NA_integer_
  } else {
    min(verification_target, .QUADRATURE_MAX_NODES)
  }
  full_refinement <- !is.na(verification_target) &&
    verification_target <= .QUADRATURE_MAX_NODES

  if (is.na(verification_order)) {
    return(list(
      moments = moments,
      jacobian = selected$jacobian,
      score_expectation = selected$score_expectation,
      diagnostics = list(
        schema_version = 1L,
        status = "approximate",
        method = "zero-limit score control variate",
        M_moments = as.integer(M),
        M_score = score_order,
        M_verification = NA_integer_,
        verification_available = FALSE,
        full_refinement = FALSE,
        reason = "quadrature_ceiling",
        max_abs_difference = NA_real_,
        max_budget_ratio = NA_real_,
        abs_tolerance = .JACOBIAN_DERIV_ABS_TOL,
        rel_tolerance = .JACOBIAN_DERIV_REL_TOL
      )
    ))
  }

  verification <- .jacobian_components_at_order(
    J, a, b, verification_order
  )
  absolute_difference <- abs(
    verification$jacobian - selected$jacobian
  )
  scale <- pmax(abs(verification$jacobian), abs(selected$jacobian))
  budget <- .JACOBIAN_DERIV_ABS_TOL + .JACOBIAN_DERIV_REL_TOL * scale
  budget_ratio <- absolute_difference / budget
  within_budget <- all(is.finite(budget_ratio)) && all(budget_ratio <= 1)
  passed <- full_refinement && within_budget
  reason <- if (!full_refinement) {
    "ceiling_limited_refinement"
  } else if (!within_budget) {
    "order_disagreement"
  } else {
    NA_character_
  }

  list(
    moments = moments,
    jacobian = verification$jacobian,
    score_expectation = verification$score_expectation,
    diagnostics = list(
      schema_version = 1L,
      status = if (passed) "converged" else "approximate",
      method = "zero-limit score control variate with order refinement",
      M_moments = as.integer(M),
      M_score = score_order,
      M_verification = verification_order,
      verification_available = TRUE,
      full_refinement = full_refinement,
      reason = reason,
      max_abs_difference = max(absolute_difference),
      max_budget_ratio = max(budget_ratio),
      abs_tolerance = .JACOBIAN_DERIV_ABS_TOL,
      rel_tolerance = .JACOBIAN_DERIV_REL_TOL,
      selected = selected$jacobian,
      verification = verification$jacobian
    )
  )
}

#' Compute Marginal Moments and Jacobian Simultaneously
#'
#' Computes quadrature approximations to the marginal moments
#' \eqn{M_1 = E[K_J]} and \eqn{V = Var(K_J)}
#' along with the Jacobian matrix of the moment map \eqn{F(a,b) = (M_1, V)}
#' using score function identities.
#'
#' @param J Integer; sample size (number of observations/sites).
#' @param a Numeric; shape parameter of the Gamma prior on \eqn{\alpha} (> 0).
#' @param b Numeric; rate parameter of the Gamma prior on \eqn{\alpha} (> 0).
#' @param M Integer; number of quadrature nodes (default: 80).
#'
#' @return A named list with components:
#' \describe{
#'   \item{\code{mean}}{Marginal mean \eqn{E[K_J]}}
#'   \item{\code{var}}{Marginal variance \eqn{Var(K_J)}}
#'   \item{\code{jacobian}}{2x2 Jacobian matrix with structure:
#'     \deqn{J_F = \begin{bmatrix} \partial M_1/\partial a & \partial M_1/\partial b \\
#'                                  \partial V/\partial a & \partial V/\partial b \end{bmatrix}}}
#'   \item{\code{derivative_diagnostics}}{Order-refinement status, tolerances,
#'     and quadrature orders used for the score derivatives.}
#'   \item{\code{conditioning}}{Scale-aware singularity diagnostics for both
#'     shape-rate and solver log-parameterizations.}
#' }
#'
#' @details
#' This function uses the score identity (Lee, 2026, Section 3.2, Corollary 1)
#' to compute derivatives without finite differences:
#' \deqn{\frac{\partial}{\partial\theta} E[f(\alpha)] = E[f(\alpha) \cdot s_\theta(\alpha)]}
#'
#' The Jacobian components are computed as:
#' \itemize{
#'   \item \eqn{\partial M_1/\partial \theta =
#'     E[(\mu_J(\alpha)-1) \cdot s_\theta(\alpha)]}
#'   \item \eqn{\partial V/\partial \theta = \partial E[v_J]/\partial \theta +
#'              \partial E[\mu_J^2]/\partial \theta - 2 M_1 \partial M_1/\partial \theta}
#' }
#' Subtracting the exact zero-concentration limits is an exact control variate
#' because \eqn{E[s_\theta(\alpha)]=0}. It removes the spurious constant times
#' \code{log(alpha)} term that otherwise dominates finite-order shape-score
#' quadrature. Derivatives are recomputed at a distinct higher order when the
#' supported ceiling permits; the returned diagnostic is \code{converged} only
#' when a full scheduled refinement is available and the two orders meet the
#' release mixed error budget
#' \code{5e-5 + 2e-3 * max(abs(selected), abs(verification))} componentwise.
#'
#' \strong{Numerical Considerations:}
#' \itemize{
#'   \item The score function \code{s_a} contains \code{log(alpha)}, which causes
#'         slower quadrature convergence compared to moment computation.
#'   \item The score calculation uses at least 80 nodes and a distinct higher
#'         order when available; both orders are retained in diagnostics.
#'   \item Conditional moments use stable finite sums without a small-alpha
#'         threshold or negative clipping.
#'   \item The returned marginal variance uses the non-negative within/between
#'         decomposition rather than cancellation-prone second-moment
#'         subtraction.
#'   \item Conditioning status is based on the solver-relevant log-parameter
#'         Jacobian, not on an unscaled determinant alone.
#' }
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @seealso
#' \code{\link{exact_K_moments}} for moments only,
#' \code{\link{score_a}}, \code{\link{score_b}} for score functions
#'
#' @examples
#' # Compute moments and Jacobian for J=50, a=2, b=1
#' result <- moments_with_jacobian(J = 50, a = 2.0, b = 1.0)
#' print(result$mean)      # E[K_J]
#' print(result$var)       # Var(K_J)
#' print(result$jacobian)  # 2x2 Jacobian matrix
#'
#' # Use in Newton iteration
#' target <- c(5.0, 8.0)  # Target (E[K], Var(K))
#' current <- c(result$mean, result$var)
#' residual <- current - target
#' delta <- solve(result$jacobian, -residual)
#'
#' @export
moments_with_jacobian <- function(J, a, b, M = .QUAD_NODES_DEFAULT) {
  assert_valid_J(J)
  a <- .dpprior_validate_scalar(
    a, "a", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  b <- .dpprior_validate_scalar(
    b, "b", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  M <- .dpprior_validate_count(
    M, "M", minimum = 10L, maximum = .QUADRATURE_MAX_NODES,
    .subclass = "dpprior_control_error"
  )

  derivative <- .jacobian_with_refinement(J, a, b, M)
  conditioning <- .jacobian_conditioning(derivative$jacobian, a, b)

  list(
    mean = derivative$moments$mean,
    var = derivative$moments$var,
    jacobian = derivative$jacobian,
    derivative_diagnostics = derivative$diagnostics,
    conditioning = conditioning
  )
}


# =============================================================================
# Verification Functions
# =============================================================================

.jacobian_five_point_reference <- function(J, a, b, M, relative_step) {
  moment_vector <- function(shape, rate) {
    result <- exact_K_moments(J, shape, rate, M)
    c(mean = result$mean, var = result$var)
  }
  step_for <- function(parameter) {
    step <- min(parameter / 4, relative_step * max(1, abs(parameter)))
    if (!is.finite(step) || step <= 0) {
      .dpprior_abort_invalid(
        "finite-difference step is outside the numerical domain",
        c("dpprior_jacobian_error", "dpprior_numerical_error"),
        "relative_step", relative_step, "positive representable step",
        "invalid_fd_step"
      )
    }
    step
  }
  five_point <- function(minus_two, minus_one, plus_one, plus_two, step) {
    (minus_two - 8 * minus_one + 8 * plus_one - plus_two) / (12 * step)
  }

  step_a <- step_for(a)
  step_b <- step_for(b)
  derivative_a <- five_point(
    moment_vector(a - 2 * step_a, b),
    moment_vector(a - step_a, b),
    moment_vector(a + step_a, b),
    moment_vector(a + 2 * step_a, b),
    step_a
  )
  derivative_b <- five_point(
    moment_vector(a, b - 2 * step_b),
    moment_vector(a, b - step_b),
    moment_vector(a, b + step_b),
    moment_vector(a, b + 2 * step_b),
    step_b
  )
  jacobian <- cbind(da = derivative_a, db = derivative_b)
  rownames(jacobian) <- c("dM1", "dV")

  list(
    jacobian = jacobian,
    step = c(a = step_a, b = step_b),
    M = as.integer(M)
  )
}

#' Verify Jacobian Against Finite Differences
#'
#' Compares the score-identity Jacobian against numerical finite differences to
#' validate the implementation.
#'
#' @param J Integer; sample size.
#' @param a Numeric; shape parameter.
#' @param b Numeric; rate parameter.
#' @param eps Positive relative step for five-point finite differences
#'   (default: 1e-6).
#' @param M Integer; number of quadrature nodes (default: 200 for verification).
#' @param verbose Logical; if TRUE, print detailed comparison.
#' @param abs_tol,rel_tol Non-negative mixed-error tolerances; at least one
#'   must be positive.
#'
#' @return A named list with components:
#' \describe{
#'   \item{\code{analytic}}{The score-identity Jacobian (legacy component
#'     name retained for compatibility)}
#'   \item{\code{numeric}}{The numerically computed Jacobian (finite differences)}
#'   \item{\code{abs_error}}{Matrix of absolute errors}
#'   \item{\code{rel_error}}{Matrix of relative errors}
#'   \item{\code{max_rel_error}}{Maximum relative error across all entries}
#'   \item{\code{status}}{One of \code{converged}, \code{approximate}, or
#'     \code{failed}.}
#'   \item{\code{component_status}}{Named statuses for order refinement and
#'     finite-difference agreement.}
#'   \item{\code{pass}}{Compatibility logical equal to
#'     \code{status == "converged"}.}
#'   \item{\code{conditioning}}{The structured conditioning result from the
#'     score-identity Jacobian.}
#' }
#'
#' @details
#' Uses the fourth-order, five-point central formula at a fixed recorded
#' quadrature order. The algebra is independent of the score identity, although
#' both computations use the package's marginal-moment integration layer.
#'
#' \strong{Important:} Release tests additionally compare both methods against
#' adaptive log-Gamma integration with separately coded conditional moments.
#'
#' @examples
#' # Verify Jacobian for a specific case
#' result <- verify_jacobian(J = 50, a = 2.0, b = 1.0, verbose = TRUE)
#'
#' @export
verify_jacobian <- function(
    J, a, b, eps = 1e-6, M = .QUAD_NODES_VERIFICATION,
    verbose = TRUE, abs_tol = .JACOBIAN_DERIV_ABS_TOL,
    rel_tol = .JACOBIAN_DERIV_REL_TOL) {
  assert_valid_J(J)
  a <- .dpprior_validate_scalar(
    a, "a", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  b <- .dpprior_validate_scalar(
    b, "b", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  eps <- .dpprior_validate_scalar(
    eps, "eps", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_control_error"
  )
  M <- .dpprior_validate_count(
    M, "M", minimum = 10L, maximum = .QUADRATURE_MAX_NODES,
    .subclass = "dpprior_control_error"
  )
  verbose <- .dpprior_validate_control(verbose, "verbose", "logical")
  abs_tol <- .dpprior_validate_scalar(
    abs_tol, "abs_tol", lower = 0,
    .subclass = "dpprior_control_error"
  )
  rel_tol <- .dpprior_validate_scalar(
    rel_tol, "rel_tol", lower = 0,
    .subclass = "dpprior_control_error"
  )
  if (abs_tol == 0 && rel_tol == 0) {
    .dpprior_abort_invalid(
      "at least one of abs_tol and rel_tol must be positive",
      c("dpprior_control_error", "dpprior_bounds_error"),
      "abs_tol/rel_tol", c(abs_tol = abs_tol, rel_tol = rel_tol),
      "at least one positive tolerance", "bounds"
    )
  }

  result <- moments_with_jacobian(J, a, b, M)
  J_analytic <- result$jacobian
  reference_order <- result$derivative_diagnostics$M_verification
  if (is.na(reference_order)) {
    reference_order <- result$derivative_diagnostics$M_score
  }
  reference <- .jacobian_five_point_reference(
    J, a, b, reference_order, eps
  )
  J_numeric <- reference$jacobian

  abs_error <- abs(J_analytic - J_numeric)
  safe_error_ratio <- function(denominator) {
    ratio <- abs_error
    positive <- denominator > 0
    ratio[positive] <- abs_error[positive] / denominator[positive]
    ratio[!positive] <- ifelse(abs_error[!positive] == 0, 0, Inf)
    ratio
  }
  relative_scale <- pmax(abs(J_numeric), abs_tol)
  rel_error <- safe_error_ratio(relative_scale)
  budget <- abs_tol + rel_tol * abs(J_numeric)
  budget_ratio <- safe_error_ratio(budget)
  finite <- all(is.finite(J_analytic)) && all(is.finite(J_numeric)) &&
    all(is.finite(budget_ratio))
  comparison_pass <- finite && all(budget_ratio <= 1)
  comparison_status <- if (!finite) {
    "failed"
  } else if (comparison_pass) {
    "converged"
  } else {
    "approximate"
  }
  refinement_status <- result$derivative_diagnostics$status
  status <- if ("failed" %in% c(comparison_status, refinement_status)) {
    "failed"
  } else if (identical(comparison_status, "converged") &&
             identical(refinement_status, "converged")) {
    "converged"
  } else {
    "approximate"
  }
  pass <- identical(status, "converged")
  max_rel_error <- max(rel_error)

  if (isTRUE(verbose)) {
    cat(sprintf("Jacobian Verification (J=%d, a=%.2f, b=%.2f, M=%d)\n", J, a, b, M))
    cat(strrep("-", 60), "\n")

    cat("\nAnalytic Jacobian (score-based):\n")
    cat(sprintf("  dM1/da = %12.8f  dM1/db = %12.8f\n",
                J_analytic[1L, 1L], J_analytic[1L, 2L]))
    cat(sprintf("  dV/da  = %12.8f  dV/db  = %12.8f\n",
                J_analytic[2L, 1L], J_analytic[2L, 2L]))

    cat("\nNumeric Jacobian (finite diff):\n")
    cat(sprintf("  dM1/da = %12.8f  dM1/db = %12.8f\n",
                J_numeric[1L, 1L], J_numeric[1L, 2L]))
    cat(sprintf("  dV/da  = %12.8f  dV/db  = %12.8f\n",
                J_numeric[2L, 1L], J_numeric[2L, 2L]))

    cat("\nRelative Errors:\n")
    cat(sprintf("  dM1/da: %.2e  dM1/db: %.2e\n",
                rel_error[1L, 1L], rel_error[1L, 2L]))
    cat(sprintf("  dV/da:  %.2e  dV/db:  %.2e\n",
                rel_error[2L, 1L], rel_error[2L, 2L]))

    cat(sprintf("\nMax Relative Error: %.2e\n", max_rel_error))
    cat(sprintf("Max Budget Ratio:  %.2e [%s]\n",
                max(budget_ratio), toupper(status)))

    if (!comparison_pass) {
      cat("\nThe score and five-point calculations disagree beyond the mixed budget.\n")
      cat("Retain approximate status and compare against the adaptive release oracle.\n")
    }
    if (!identical(refinement_status, "converged")) {
      cat(sprintf(
        "\nScore-order refinement status is %s; overall convergence is withheld.\n",
        refinement_status
      ))
    }
  }

  invisible(list(
    analytic = J_analytic,
    numeric = J_numeric,
    abs_error = abs_error,
    rel_error = rel_error,
    budget = budget,
    budget_ratio = budget_ratio,
    max_abs_error = max(abs_error),
    max_rel_error = max_rel_error,
    max_budget_ratio = max(budget_ratio),
    status = status,
    component_status = c(
      order_refinement = refinement_status,
      finite_difference = comparison_status
    ),
    pass = pass,
    reference_M = reference$M,
    step = reference$step,
    conditioning = result$conditioning
  ))
}


#' Verify Score Function Zero Expectation Property
#'
#' Verifies the fundamental property that \eqn{E[s_\theta(\alpha)] = 0}
#' for both score functions.
#'
#' @param a Numeric; shape parameter.
#' @param b Numeric; rate parameter.
#' @param M Integer; number of quadrature nodes.
#' @param verbose Logical; if TRUE, print results.
#' @param tolerance Non-negative absolute tolerance for the raw quadrature
#'   expectation of each score.
#'
#' @return A named list with components:
#' \describe{
#'   \item{\code{E_score_a}}{Expectation of \eqn{s_a}}
#'   \item{\code{E_score_b}}{Expectation of \eqn{s_b}}
#'   \item{\code{status}}{\code{converged} when both raw score expectations
#'     meet \code{tolerance}; otherwise \code{approximate}.}
#'   \item{\code{abs_error}}{Named absolute raw score-expectation errors.}
#'   \item{\code{passed}}{Compatibility logical equal to
#'     \code{status == "converged"}.}
#' }
#'
#' @details
#' This is a fundamental property of score functions. Due to quadrature
#' approximation error, the computed expectations may not be exactly zero.
#'
#' \strong{Expected behavior:}
#' \itemize{
#'   \item \code{E[s_b]} should be very close to zero (typically < 1e-14)
#'         because s_b is linear in alpha.
#'   \item \code{E[s_a]} can show material error for small shape because the
#'         \code{log(alpha)} singularity converges slowly under a finite rule.
#' }
#'
#' The production Jacobian uses an exact zero-limit control variate, so raw
#' score-expectation error is reported rather than silently added to every
#' derivative. Release tests use adaptive log-Gamma integration as an
#' independent score-identity reference.
#'
#' @examples
#' \dontrun{
#' verify_score_expectation(a = 2.0, b = 1.0, verbose = TRUE)
#'
#' }
#' @keywords internal
verify_score_expectation <- function(
    a, b, M = .QUAD_NODES_VERIFICATION, verbose = TRUE,
    tolerance = 1e-3) {
  a <- .dpprior_validate_scalar(
    a, "a", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  b <- .dpprior_validate_scalar(
    b, "b", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  M <- .dpprior_validate_count(
    M, "M", minimum = 1L, maximum = .QUADRATURE_MAX_NODES,
    .subclass = "dpprior_control_error"
  )
  verbose <- .dpprior_validate_control(verbose, "verbose", "logical")
  tolerance <- .dpprior_validate_scalar(
    tolerance, "tolerance", lower = 0,
    .subclass = "dpprior_control_error"
  )

  quad <- build_gamma_quadrature(a, b, M)
  alphas <- quad$alpha_nodes
  w <- quad$weights_normalized

  E_sa <- .quadrature_weighted_sum(w, score_a(alphas, a, b))
  E_sb <- .quadrature_weighted_sum(w, score_b(alphas, a, b))
  absolute_error <- c(da = abs(E_sa), db = abs(E_sb))
  passed <- all(is.finite(absolute_error)) && all(absolute_error <= tolerance)
  status <- if (passed) "converged" else "approximate"

  if (isTRUE(verbose)) {
    cat(sprintf("Score Expectation Verification (a=%.2f, b=%.2f, M=%d)\n", a, b, M))
    cat(sprintf("  E[s_a(alpha)] = %.6e (should be ~= 0)\n", E_sa))
    cat(sprintf("  E[s_b(alpha)] = %.6e (should be ~= 0)\n", E_sb))

    if (!passed) {
      cat("\n  Note: E[s_a] has larger error because log(alpha) slows quadrature.\n")
      cat("  Raw score integration is marked approximate; the Jacobian uses\n")
      cat("  a zero-limit control variate and independent order refinement.\n")
    }
  }

  invisible(list(
    E_score_a = E_sa,
    E_score_b = E_sb,
    abs_error = absolute_error,
    tolerance = tolerance,
    status = status,
    passed = passed,
    M = M
  ))
}


#' Run All Module 07 Verification Tests
#'
#' Comprehensive verification suite for the score-based Jacobian module.
#'
#' @param verbose Logical; if TRUE, print detailed results.
#'
#' @return Logical; TRUE if all tests pass.
#'
#' @examples
#' \dontrun{
#' verify_jacobian_all()
#'
#' }
#' @keywords internal
verify_jacobian_all <- function(verbose = TRUE) {
  if (isTRUE(verbose)) {
    cat(strrep("=", 70), "\n")
    cat("Module 07: Score-Based Jacobian - Full Verification Suite\n")
    cat(strrep("=", 70), "\n\n")
  }

  all_pass <- TRUE

  # Test 1: Score function expectations
  if (isTRUE(verbose)) {
    cat("[Test 1] Score function zero-expectation property\n")
    cat(strrep("-", 50), "\n")
    cat("Note: E[s_a] may show larger errors due to log(alpha) term.\n\n")
  }

  test_cases <- list(
    c(a = 2.0, b = 1.0),
    c(a = 1.5, b = 0.5),
    c(a = 3.0, b = 1.5)
  )

  for (tc in test_cases) {
    score_result <- verify_score_expectation(
      tc["a"], tc["b"], verbose = verbose
    )
    all_pass <- all_pass && score_result$passed
    if (isTRUE(verbose)) cat("\n")
  }

  # Test 2: Jacobian matches finite differences
  if (isTRUE(verbose)) {
    cat("[Test 2] Jacobian vs finite differences (M=200)\n")
    cat(strrep("-", 50), "\n\n")
  }

  test_cases <- list(
    list(J = 30, a = 1.5, b = 0.5),
    list(J = 50, a = 2.0, b = 1.0),
    list(J = 100, a = 2.0, b = 1.0),
    list(J = 50, a = 3.0, b = 1.5)
  )

  for (tc in test_cases) {
    result <- verify_jacobian(tc$J, tc$a, tc$b, verbose = verbose)
    all_pass <- all_pass && result$pass
    if (isTRUE(verbose)) cat("\n")
  }

  # Test 3: Moments consistency
  if (isTRUE(verbose)) {
    cat("[Test 3] Moments consistency check\n")
    cat(strrep("-", 50), "\n")
  }

  J <- 50; a <- 2.0; b <- 1.0
  result_jac <- moments_with_jacobian(J, a, b)
  result_mom <- exact_K_moments(J, a, b)

  mean_diff <- abs(result_jac$mean - result_mom$mean)
  var_diff <- abs(result_jac$var - result_mom$var)
  consistent <- mean_diff < 1e-10 && var_diff < 1e-10

  if (isTRUE(verbose)) {
    cat(sprintf("Mean (moments_with_jacobian): %.10f\n", result_jac$mean))
    cat(sprintf("Mean (exact_K_moments):       %.10f\n", result_mom$mean))
    cat(sprintf("Difference: %.2e\n", mean_diff))
    cat(sprintf("Var (moments_with_jacobian):  %.10f\n", result_jac$var))
    cat(sprintf("Var (exact_K_moments):        %.10f\n", result_mom$var))
    cat(sprintf("Difference: %.2e\n", var_diff))
    cat(sprintf("Status: %s\n\n", if (consistent) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && consistent

  # Test 4: Jacobian non-singularity
  if (isTRUE(verbose)) {
    cat("[Test 4] Jacobian non-singularity check\n")
    cat(strrep("-", 50), "\n")
  }

  test_cases <- list(
    list(J = 30, a = 1.5, b = 0.5),
    list(J = 50, a = 2.0, b = 1.0),
    list(J = 100, a = 2.0, b = 1.0)
  )

  for (tc in test_cases) {
    result <- moments_with_jacobian(tc$J, tc$a, tc$b)
    condition <- result$conditioning$solver
    nonsingular <- identical(condition$status, "well_conditioned")

    if (isTRUE(verbose)) {
      cat(sprintf("J=%d, a=%.1f, b=%.1f: det=%.4f, cond=%.2e [%s]\n",
                  tc$J, tc$a, tc$b, condition$determinant,
                  condition$condition_number,
                  if (nonsingular) "PASS" else "FAIL"))
    }
    all_pass <- all_pass && nonsingular
  }

  # Test 5: Newton convergence
  if (isTRUE(verbose)) {
    cat("\n[Test 5] Newton convergence test\n")
    cat(strrep("-", 50), "\n")
  }

  newton_pass <- test_newton_convergence(J = 50, mu_target = 5.0, var_target = 8.0,
                                         verbose = verbose)
  all_pass <- all_pass && newton_pass

  # Summary
  if (isTRUE(verbose)) {
    cat("\n", strrep("=", 70), "\n", sep = "")
    cat(sprintf("Overall Result: %s\n",
                if (all_pass) "ALL TESTS PASSED" else "SOME TESTS FAILED"))
    cat(strrep("=", 70), "\n")
  }

  invisible(all_pass)
}


#' Test Newton Convergence Using the Jacobian
#'
#' Verifies that the Jacobian enables fast Newton convergence for moment matching.
#'
#' @param J Integer; sample size.
#' @param mu_target Numeric; target mean.
#' @param var_target Numeric; target variance.
#' @param a0 Numeric; initial shape parameter.
#' @param b0 Numeric; initial rate parameter.
#' @param max_iter Integer; maximum iterations.
#' @param tol Numeric; convergence tolerance.
#' @param verbose Logical; if TRUE, print iteration history.
#'
#' @return Logical; TRUE if Newton converges.
#'
#' @keywords internal
test_newton_convergence <- function(J, mu_target, var_target,
                                    a0 = 2.0, b0 = 1.0,
                                    max_iter = 15L, tol = 1e-8,
                                    verbose = TRUE) {
  a <- a0
  b <- b0

  if (isTRUE(verbose)) {
    cat(sprintf("Target: E[K]=%.2f, Var(K)=%.2f\n", mu_target, var_target))
    cat(sprintf("Start:  a=%.4f, b=%.4f\n\n", a, b))
    cat("Iter |    a      |    b      |   E[K]    |  Var(K)   | ||F||\n")
    cat(strrep("-", 70), "\n")
  }

  for (iter in 0:max_iter) {
    result <- moments_with_jacobian(J, a, b, M = .QUAD_NODES_VERIFICATION)
    F_vec <- c(result$mean - mu_target, result$var - var_target)
    norm_F <- sqrt(sum(F_vec^2))

    if (isTRUE(verbose)) {
      cat(sprintf("%4d | %9.6f | %9.6f | %9.6f | %9.6f | %.2e\n",
                  iter, a, b, result$mean, result$var, norm_F))
    }

    if (norm_F < tol) {
      if (isTRUE(verbose)) {
        cat(sprintf("\nConverged in %d iterations!\n", iter))
      }
      return(TRUE)
    }

    delta <- tryCatch(
      solve(result$jacobian, -F_vec),
      error = function(e) NULL
    )

    if (is.null(delta)) {
      if (isTRUE(verbose)) cat("Jacobian singular!\n")
      return(FALSE)
    }

    # Damped update
    step <- 1.0
    while (step > 1e-10) {
      a_new <- a + step * delta[1]
      b_new <- b + step * delta[2]
      if (a_new > 0.1 && b_new > 0.01) break
      step <- step * 0.5
    }

    a <- a_new
    b <- b_new
  }

  if (isTRUE(verbose)) cat("Max iterations reached.\n")
  return(FALSE)
}

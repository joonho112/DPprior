# =============================================================================
# Module 09: Co-Clustering Probability (rho) Distribution
# =============================================================================
#
# Author: JoonHo Lee (jlee296@ua.edu)
# Date: December 2025
# Part of: DPprior R Package
#
# This module implements moments of the co-clustering probability
#   rho = sum_{h>=1} w_h^2,
# where (w_h) are the stick-breaking weights under GEM(alpha).
#
# Probabilistic interpretation:
#   rho = P(Z1 = Z2 | w) where Z1, Z2 are cluster labels for two random units.
#
# Conditional moments given alpha (Lee, 2026, Section 4):
#   E[rho | alpha] = 1 / (1 + alpha)
#   E[rho^2 | alpha] = (alpha + 6) / ((alpha+1)(alpha+2)(alpha+3))
#   Var(rho | alpha) = 2 alpha / ((alpha+1)^2 (alpha+2)(alpha+3))
#
# Key identity: E[rho | alpha] = E[w1 | alpha] = 1/(1+alpha)
#   Therefore E[rho | a, b] = E[w1 | a, b], but full distributions differ.
#
# Marginal moments under alpha ~ Gamma(a, b) are computed via Gauss-Laguerre
# quadrature using integrate_gamma() (Module 02).
# =============================================================================


# =============================================================================
# Conditional Moments Given alpha
# =============================================================================

#' Conditional Mean of rho Given Alpha
#'
#' Computes the conditional mean E(rho | alpha) for the co-clustering
#' probability rho = sum_h w_h^2 under a Dirichlet Process.
#'
#' @param alpha Numeric vector; concentration parameter(s) (must be positive).
#'
#' @return Numeric vector; `E(rho | alpha) = 1/(1+alpha)`.
#'
#' @details
#' The co-clustering probability rho = sum(w_h^2) over h >= 1 has conditional mean:
#' \deqn{E[\rho | \alpha] = \frac{1}{1 + \alpha}}
#'
#' This equals `E(w1 | alpha)` since w1 ~ Beta(1, alpha) has mean 1/(1+alpha).
#'
#' \strong{Interpretation:}
#' \itemize{
#'   \item alpha -> 0: E(rho|alpha) -> 1 (all observations in one cluster)
#'   \item alpha -> Inf: E(rho|alpha) -> 0 (infinitely many small clusters)
#'   \item alpha = 1: E(rho|alpha) = 0.5 (moderate clustering)
#' }
#'
#' @examples
#' mean_rho_given_alpha(1.0)
#' mean_rho_given_alpha(c(0.5, 1, 2, 5, 10))
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @seealso \code{\link{var_rho_given_alpha}}, \code{\link{mean_rho}}
#'
#' @family co_clustering
#'
#' @export
mean_rho_given_alpha <- function(alpha) {

  assert_positive(alpha, "alpha")
  1 / (1 + alpha)
}


#' Conditional Second Moment of rho Given Alpha
#'
#' Computes `E(rho^2 | alpha) = (alpha + 6) / ((alpha+1)(alpha+2)(alpha+3))`.
#'
#' @param alpha Numeric vector; concentration parameter(s) (must be positive).
#'
#' @return Numeric vector; `E(rho^2 | alpha)`.
#'
#' @details
#' Used internally for variance computation via the identity:
#' `Var(rho|alpha) = E(rho^2|alpha) - E(rho|alpha)^2`
#'
#' @keywords internal
mean_rho_sq_given_alpha <- function(alpha) {
  assert_positive(alpha, "alpha")

  # Evaluate as a product of bounded ratios.  The algebraically equivalent
  # direct product in the denominator overflows around sqrt(.Machine$xmax),
  # long before the second moment itself underflows.
  (1 / (alpha + 1)) * ((alpha + 6) / (alpha + 2)) /
    (alpha + 3)
}


#' Conditional Variance of rho Given Alpha
#'
#' Computes the conditional variance Var(rho | alpha).
#'
#' @param alpha Numeric vector; concentration parameter(s) (must be positive).
#'
#' @return Numeric vector; Var(rho | alpha).
#'
#' @details
#' The conditional variance is:
#' \deqn{Var(\rho | \alpha) = \frac{2\alpha}{(1+\alpha)^2(2+\alpha)(3+\alpha)}}
#'
#' This is derived from the GEM recursion:
#' rho = V^2 + (1-V)^2 * rho' where V ~ Beta(1, alpha) and rho' is an
#' independent copy of rho.
#'
#' \strong{Properties:}
#' \itemize{
#'   \item Var(rho|alpha) = 0 when alpha -> 0 (degenerate at rho = 1)
#'   \item Var(rho|alpha) -> 0 when alpha -> Inf (degenerate at rho = 0)
#'   \item Maximum variance occurs at intermediate alpha
#' }
#'
#' @examples
#' var_rho_given_alpha(2)
#' var_rho_given_alpha(c(0.5, 1, 2, 5, 10))
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @seealso \code{\link{mean_rho_given_alpha}}, \code{\link{var_rho}}
#'
#' @family co_clustering
#'
#' @export
var_rho_given_alpha <- function(alpha) {
  assert_positive(alpha, "alpha")

  # Keep every intermediate bounded.  In particular, avoid forming
  # (alpha + 1)^2 * (alpha + 2) * (alpha + 3), which can overflow while the
  # mathematically positive variance is still representable.
  2 * (alpha / (alpha + 1)) / (alpha + 1) /
    (alpha + 2) / (alpha + 3)
}


# =============================================================================
# Marginal Moments (alpha ~ Gamma(a, b))
# =============================================================================

#' Marginal Mean of rho
#'
#' Computes E(rho | a, b) when alpha ~ Gamma(a, b) (shape-rate).
#'
#' @param a Numeric; shape parameter of the Gamma prior on alpha (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior on alpha (b > 0).
#' @param M Integer; number of quadrature nodes. Default is 80.
#'
#' @return Numeric; `E(rho | a, b)`.
#'
#' @details
#' Uses Gauss-Laguerre quadrature via \code{integrate_gamma}. A key identity is
#' `E(rho | alpha) = E(w1 | alpha) = 1/(1+alpha)`, so `E(rho | a, b)` equals
#' `E(w1 | a, b)` (but the full distributions differ).
#'
#' The result is the marginal probability that two exchangeable units share a
#' cluster, averaged over both the random DP weights and the Gamma hyperprior.
#' No qualitative category or action threshold is imposed by this function.
#'
#' @examples
#' mean_rho(a = 2, b = 1)
#' mean_rho(a = 1.6, b = 1.22)
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @seealso \code{\link{var_rho}}, \code{\link{cv_rho}}, \code{\link{mean_w1}}
#'
#' @family co_clustering
#'
#' @export
mean_rho <- function(a, b, M = .QUAD_NODES_DEFAULT) {
  assert_positive(a, "a")
  assert_positive(b, "b")

  integrate_gamma(function(alpha) 1 / (1 + alpha), a, b, M)
}


#' Marginal Second Moment of rho
#'
#' Computes `E(rho^2 | a, b)` by mixing `E(rho^2 | alpha)` over alpha ~ Gamma(a, b).
#'
#' @param a Numeric; shape parameter of the Gamma prior on alpha (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior on alpha (b > 0).
#' @param M Integer; number of quadrature nodes. Default is 80.
#'
#' @return Numeric; `E(rho^2 | a, b)`.
#'
#' @keywords internal
mean_rho_sq <- function(a, b, M = .QUAD_NODES_DEFAULT) {
  assert_positive(a, "a")
  assert_positive(b, "b")

  integrate_gamma(mean_rho_sq_given_alpha, a, b, M)
}


#' Marginal Variance of rho
#'
#' Computes Var(rho | a, b) when alpha ~ Gamma(a, b) (shape-rate).
#'
#' @param a Numeric; shape parameter of the Gamma prior on alpha (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior on alpha (b > 0).
#' @param M Integer; number of quadrature nodes. Default is 80.
#'
#' @return Numeric; Var(rho | a, b).
#'
#' @details
#' Uses the law of total variance:
#' \deqn{Var(\rho | a, b) = E[Var(\rho | \alpha)] + Var(E[\rho | \alpha])}
#'
#' where:
#' \itemize{
#'   \item `Var(rho | alpha) = 2*alpha / ((1+alpha)^2*(2+alpha)*(3+alpha))`
#'   \item `E(rho | alpha) = 1/(1+alpha)`
#' }
#'
#' \strong{Note:} Unlike E(rho), Var(rho) != Var(w1) in general, because
#' the conditional variances differ.
#'
#' @examples
#' var_rho(a = 2, b = 1)
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @seealso \code{\link{mean_rho}}, \code{\link{cv_rho}}, \code{\link{var_w1}}
#'
#' @family co_clustering
#'
#' @export
var_rho <- function(a, b, M = .QUAD_NODES_DEFAULT) {
  assert_positive(a, "a")
  assert_positive(b, "b")

  # Use a centered law-of-total-variance calculation.  Subtracting
  # E[m(alpha)^2] - E[m(alpha)]^2 loses the between-alpha component when the
  # hyperprior is concentrated, even though every component is non-negative.
  quad <- build_gamma_quadrature(a, b, M)
  weights <- quad$weights_normalized
  conditional_mean <- mean_rho_given_alpha(quad$alpha_nodes)
  marginal_mean <- .quadrature_weighted_sum(weights, conditional_mean)
  within_alpha <- .quadrature_weighted_sum(
    weights, var_rho_given_alpha(quad$alpha_nodes)
  )
  between_alpha <- .quadrature_weighted_sum(
    weights, (conditional_mean - marginal_mean)^2
  )

  as.numeric(within_alpha + between_alpha)
}


#' Coefficient of Variation of rho
#'
#' Computes CV(rho) = SD(rho) / E(rho) under alpha ~ Gamma(a, b).
#'
#' @param a Numeric; shape parameter of the Gamma prior on alpha (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior on alpha (b > 0).
#' @param M Integer; number of quadrature nodes. Default is 80.
#'
#' @return Numeric; coefficient of variation.
#'
#' @examples
#' cv_rho(a = 2, b = 1)
#' cv_rho(a = 1.6, b = 1.22)
#'
#' @seealso \code{\link{mean_rho}}, \code{\link{var_rho}}
#'
#' @family co_clustering
#'
#' @export
cv_rho <- function(a, b, M = .QUAD_NODES_DEFAULT) {
  sqrt(var_rho(a, b, M)) / mean_rho(a, b, M)
}


# =============================================================================
# Summary Functions
# =============================================================================

#' Summary Statistics for rho Distribution
#'
#' Computes comprehensive summary statistics for the co-clustering probability
#' rho under the hierarchical prior alpha ~ Gamma(a, b).
#'
#' @param a Numeric; shape parameter of the Gamma prior on alpha (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior on alpha (b > 0).
#' @param M Integer; selected quadrature order. Default is 80.
#' @param M_verify Optional independent verification order. When omitted and
#'   available, the package-wide minimum refinement order is used.
#' @param abs_tol,rel_tol Non-negative selected-versus-verification tolerances.
#' @param strict Logical; if \code{TRUE}, require verified rho moments.
#'
#' @return A list of class "rho_summary" containing:
#'   \describe{
#'     \item{mean}{E(rho | a, b)}
#'     \item{var}{Var(rho | a, b)}
#'     \item{sd}{SD(rho | a, b) = sqrt(Var)}
#'     \item{cv}{Coefficient of variation SD/mean}
#'     \item{status, usable, verified, message}{Numerical verification state}
#'     \item{estimand, label, units}{Explicit estimand metadata}
#'     \item{params}{List of input parameters (a, b)}
#'     \item{alpha_prior}{Summary of the alpha prior (mean, sd, cv)}
#'     \item{conditional_at_alpha_mean}{Conditional moments evaluated at E(alpha)}
#'     \item{verification, provenance}{Quadrature comparison and method metadata}
#'   }
#'
#' @details
#' The co-clustering probability rho indicates how likely two randomly
#' chosen observations are to belong to the same cluster a priori.
#'
#' The \code{conditional_at_alpha_mean} component provides a "plug-in" estimate
#' for comparison: what the moments would be if alpha were fixed at its prior
#' mean. The function reports numerical values and verification status only;
#' scientific action thresholds must be supplied explicitly by an analysis
#' policy rather than inferred from uncalibrated labels.
#'
#' @examples
#' summary_rho(a = 2, b = 1)
#' summary_rho(a = 1.6, b = 1.22)
#'
#' @seealso \code{\link{mean_rho}}, \code{\link{var_rho}}, \code{\link{summary_w1}}
#'
#' @export
summary_rho <- function(
    a, b, M = .QUAD_NODES_DEFAULT, M_verify = NULL,
    abs_tol = 1e-10, rel_tol = 1e-8, strict = FALSE) {
  a <- .dpprior_validate_scalar(
    a, "a", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  b <- .dpprior_validate_scalar(
    b, "b", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  controls <- .marginal_verification_controls(
    M, M_verify, abs_tol, rel_tol, strict = FALSE
  )
  if (is.null(controls$M_verify) && controls$verification_available) {
    controls <- .marginal_verification_controls(
      controls$M,
      M_verify = controls$M_verification_required,
      abs_tol = controls$abs_tol, rel_tol = controls$rel_tol,
      strict = strict
    )
  } else if (strict) {
    controls <- .marginal_verification_controls(
      controls$M, controls$M_verify,
      controls$abs_tol, controls$rel_tol, strict = TRUE
    )
  }

  # Marginal moments
  mean_val <- mean_rho(a, b, controls$M)
  var_val <- var_rho(a, b, controls$M)
  sd_val <- sqrt(var_val)
  mean_verification <- var_verification <- NULL
  if (!is.null(controls$M_verify)) {
    mean_verification <- mean_rho(a, b, controls$M_verify)
    var_verification <- var_rho(a, b, controls$M_verify)
  }
  audit <- function(selected, verification) {
    difference <- if (is.null(verification)) NA_real_ else {
      abs(selected - verification)
    }
    tolerance <- if (is.null(verification)) NA_real_ else {
      controls$abs_tol + controls$rel_tol *
        max(abs(selected), abs(verification))
    }
    list(
      selected = as.numeric(selected),
      verification = if (is.null(verification)) NA_real_ else {
        as.numeric(verification)
      },
      difference = difference,
      tolerance = tolerance,
      passed = if (is.null(verification)) NA else {
        is.finite(difference) && difference <= tolerance
      }
    )
  }
  mean_audit <- audit(mean_val, mean_verification)
  var_audit <- audit(var_val, var_verification)
  verified <- isTRUE(mean_audit$passed) && isTRUE(var_audit$passed)
  status <- if (verified) "converged" else "approximate"
  message <- if (verified) {
    "rho mean and variance passed higher-order quadrature verification"
  } else {
    "rho moments are unverified or disagree at the required refinement order"
  }
  if (strict && !verified) {
    stop(.dpprior_new_condition(
      message,
      classes = c(
        "dpprior_rho_convergence_error", "dpprior_numerical_error",
        "dpprior_error", "error"
      ),
      result = list(mean = mean_audit, variance = var_audit),
      code = "rho_quadrature_disagreement"
    ))
  }

  # alpha prior summary
  alpha_mean <- a / b
  alpha_sd <- sqrt(a) / b
  alpha_cv <- 1 / sqrt(a)

  # Conditional moments at E[alpha] (plug-in estimate)
  cond_mean <- mean_rho_given_alpha(alpha_mean)
  cond_var <- var_rho_given_alpha(alpha_mean)

  result <- list(
    status = status,
    usable = verified,
    verified = verified,
    message = message,
    estimand = "rho",
    label = "Conditional pairwise co-clustering probability",
    units = "probability",
    mean = mean_val,
    var = var_val,
    sd = sd_val,
    cv = sd_val / mean_val,
    params = list(a = a, b = b),
    alpha_prior = list(
      mean = alpha_mean,
      sd = alpha_sd,
      cv = alpha_cv
    ),
    conditional_at_alpha_mean = list(
      alpha = alpha_mean,
      mean = cond_mean,
      var = cond_var
    ),
    method = "gauss-laguerre-marginal-moments",
    verification = list(mean = mean_audit, variance = var_audit),
    provenance = list(
      schema_version = 2L,
      selected_method = "gauss_laguerre",
      selected_order = controls$M,
      verification_order = controls$M_verify,
      verification_order_required = controls$M_verification_required,
      verification_available = controls$verification_available,
      is_fallback = FALSE
    )
  )

  class(result) <- "rho_summary"
  result
}


#' Print Method for rho_summary Objects
#'
#' @param x An object of class "rho_summary".
#' @param digits Integer; number of digits for printing.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.rho_summary <- function(x, digits = 4, ...) {
  cat("Co-Clustering Probability (rho) Summary\n")
  cat(strrep("=", 50), "\n\n")

  cat(sprintf("Status:      %s\n", toupper(x$status)))
  cat(sprintf("Verified:    %s\n", if (isTRUE(x$verified)) "yes" else "no"))
  cat(sprintf("Estimand:    %s (%s)\n", x$estimand, x$label))
  cat(sprintf("Method:      %s\n", x$method))
  if (!isTRUE(x$verified)) {
    cat(sprintf("Caveat:      %s\n", x$message))
  }
  cat("\n")

  cat(sprintf("Gamma prior: alpha ~ Gamma(%.4f, %.4f)\n",
              x$params$a, x$params$b))
  cat(sprintf("E[alpha] = %.4f, SD(alpha) = %.4f, CV(alpha) = %.1f%%\n\n",
              x$alpha_prior$mean, x$alpha_prior$sd, 100 * x$alpha_prior$cv))

  cat("Marginal distribution of rho:\n")
  cat(strrep("-", 35), "\n")
  cat(sprintf("  Mean:   %.*f\n", digits, x$mean))
  cat(sprintf("  SD:     %.*f\n", digits, x$sd))
  cat(sprintf("  CV:     %.1f%%\n", 100 * x$cv))

  cat(sprintf("\nConditional at E[alpha] = %.4f (plug-in):\n", x$conditional_at_alpha_mean$alpha))
  cat(strrep("-", 35), "\n")
  cat(sprintf("  E[rho | E[alpha]]:   %.*f\n", digits, x$conditional_at_alpha_mean$mean))
  cat(sprintf("  Var(rho | E[alpha]): %.*f\n", digits, x$conditional_at_alpha_mean$var))

  invisible(x)
}


# =============================================================================
# Verification Functions
# =============================================================================

#' Verify the Identity E(w1) = E(rho)
#'
#' Checks the mean identity `E(w1 | a, b) = E(rho | a, b)`, which follows from
#' `E(rho | alpha) = E(w1 | alpha) = 1/(1+alpha)`.
#'
#' @param a Numeric; shape parameter of the Gamma prior on alpha (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior on alpha (b > 0).
#' @param tol Numeric; absolute tolerance. Default is 1e-10.
#' @param M Integer; number of quadrature nodes. Default is 80.
#'
#' @return Logical; TRUE if the identity holds within tolerance.
#'
#' @examples
#' \dontrun{
#' verify_w1_rho_identity(2, 1)
#' verify_w1_rho_identity(1.6, 1.22)
#' }
#'
#' @keywords internal
verify_w1_rho_identity <- function(a, b, tol = 1e-10, M = .QUAD_NODES_DEFAULT) {
  assert_positive(a, "a")
  assert_positive(b, "b")

  abs(mean_w1(a, b, M) - mean_rho(a, b, M)) < tol
}


#' Verify Variance Decomposition
#'
#' Verifies the law of total variance decomposition for Var(rho).
#'
#' @param a Numeric; shape parameter of the Gamma prior on alpha (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior on alpha (b > 0).
#' @param tol Numeric; tolerance for comparison. Default is 1e-10.
#' @param M Integer; number of quadrature nodes. Default is 80.
#' @param verbose Logical; if TRUE, print detailed results.
#'
#' @return Logical; TRUE if decomposition holds within tolerance.
#'
#' @keywords internal
verify_rho_variance_decomposition <- function(a, b, tol = 1e-10,
                                              M = .QUAD_NODES_DEFAULT,
                                              verbose = FALSE) {
  # Total variance via var_rho()
  var_total <- var_rho(a, b, M)

  # E[Var(rho|alpha)]
  E_var_cond <- integrate_gamma(var_rho_given_alpha, a, b, M)

  # Var(E[rho|alpha]) = E[(1/(1+alpha))^2] - (E[1/(1+alpha)])^2
  E_mean_sq <- integrate_gamma(function(alpha) (1 / (1 + alpha))^2, a, b, M)
  E_mean <- mean_rho(a, b, M)
  var_mean_cond <- E_mean_sq - E_mean^2

  var_decomposed <- E_var_cond + var_mean_cond
  passed <- abs(var_total - var_decomposed) < tol

  if (verbose) {
    cat(sprintf("Variance decomposition (a=%.2f, b=%.2f):\n", a, b))
    cat(sprintf("  E[Var(rho|alpha)] = %.10f\n", E_var_cond))
    cat(sprintf("  Var(E[rho|alpha]) = %.10f\n", var_mean_cond))
    cat(sprintf("  Sum               = %.10f\n", var_decomposed))
    cat(sprintf("  Var(rho)          = %.10f\n", var_total))
    cat(sprintf("  Match: %s\n", if (passed) "PASS" else "FAIL"))
  }

  passed
}


#' Verify Conditional Variance Formula
#'
#' Verifies that `Var(rho|alpha) = E(rho^2|alpha) - E(rho|alpha)^2`.
#'
#' @param alpha Numeric; concentration parameter (must be positive).
#' @param tol Numeric; tolerance for comparison. Default is 1e-12.
#'
#' @return Logical; TRUE if formula holds.
#'
#' @keywords internal
verify_rho_conditional_variance <- function(alpha, tol = 1e-12) {
  E_rho <- mean_rho_given_alpha(alpha)
  E_rho_sq <- mean_rho_sq_given_alpha(alpha)
  var_direct <- var_rho_given_alpha(alpha)
  var_from_moments <- E_rho_sq - E_rho^2

  abs(var_direct - var_from_moments) < tol
}


# =============================================================================
# Comparison Functions
# =============================================================================

#' Compare rho and w1 Distributions
#'
#' Compares the marginal distributions of rho and w1 under the same
#' hyperprior alpha ~ Gamma(a, b).
#'
#' @param a Numeric; shape parameter of the Gamma prior on alpha (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior on alpha (b > 0).
#' @param M Integer; number of quadrature nodes. Default is 80.
#'
#' @return A list containing:
#'   \describe{
#'     \item{mean_rho}{E(rho)}
#'     \item{mean_w1}{E(w1)}
#'     \item{mean_equal}{Logical; whether means are equal}
#'     \item{var_rho}{Var(rho)}
#'     \item{var_w1}{Var(w1)}
#'     \item{var_ratio}{Var(rho) / Var(w1)}
#'   }
#'
#' @details
#' While E(rho) = E(w1), the variances differ because:
#' \itemize{
#'   \item `Var(rho | alpha) = 2*alpha / ((1+alpha)^2*(2+alpha)*(3+alpha))`
#'   \item `Var(w1 | alpha) = alpha / ((1+alpha)^2*(2+alpha))`
#' }
#'
#' Generally, Var(rho) < Var(w1) because rho averages over all squared weights.
#'
#' @examples
#' \dontrun{
#' compare_rho_w1(a = 2, b = 1)
#' compare_rho_w1(a = 1.6, b = 1.22)
#'
#' }
#' @keywords internal
compare_rho_w1 <- function(a, b, M = .QUAD_NODES_DEFAULT) {
  assert_positive(a, "a")
  assert_positive(b, "b")

  mean_rho_val <- mean_rho(a, b, M)
  mean_w1_val <- mean_w1(a, b, M)
  var_rho_val <- var_rho(a, b, M)
  var_w1_val <- var_w1(a, b, M)

  list(
    mean_rho = mean_rho_val,
    mean_w1 = mean_w1_val,
    mean_equal = abs(mean_rho_val - mean_w1_val) < 1e-10,
    var_rho = var_rho_val,
    var_w1 = var_w1_val,
    var_ratio = var_rho_val / var_w1_val
  )
}


# =============================================================================
# Random Generation (for Monte Carlo validation)
# =============================================================================

.RRHO_MAX_STICKS <- 10000000L
.RRHO_BATCH_SIZE <- 4096L
.RRHO_SMALLEST_POSITIVE <- .Machine$double.xmin * .Machine$double.eps

.rrho_bound_from_log <- function(log_bound) {
  if (log_bound == -Inf) return(0)
  if (log_bound < log(.RRHO_SMALLEST_POSITIVE)) {
    return(.RRHO_SMALLEST_POSITIVE)
  }
  exp(log_bound)
}

#' Random Generation from rho Distribution
#'
#' Generates random samples from the rho = sum w_h^2 distribution by
#' stick-breaking simulation.
#'
#' @param n Integer; number of samples to generate.
#' @param a Numeric; shape parameter of the Gamma prior on alpha (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior on alpha (b > 0).
#' @param n_sticks Integer; minimum number of sticks generated for each draw.
#'   The sampler continues adaptively when the remaining squared mass is above
#'   \code{remainder_tol}. Default is 500.
#' @param remainder_tol Numeric; maximum allowed deterministic upper bound on
#'   the omitted contribution to rho. Because all ungenerated weights sum to
#'   remainder \eqn{R}, their contribution satisfies
#'   \eqn{0 <= \sum_{\mathrm{tail}} w_h^2 <= R^2}. Default is \code{1e-10}.
#' @param max_sticks Integer; hard per-draw stick ceiling. Default is 100,000.
#' @param strict Logical; if \code{TRUE} (the default), reaching
#'   \code{max_sticks} before the remainder contract is met raises a typed
#'   \code{dpprior_rho_truncation_error}. If \code{FALSE}, the finite lower
#'   approximation is returned with status \code{"approximate"} and its
#'   deterministic remainder bound in the \code{"rrho_diagnostics"}
#'   attribute.
#'
#' @return Numeric vector of length n. The \code{"rrho_diagnostics"} attribute
#'   records the adaptive method, status, tolerance, sticks used, and the
#'   per-draw deterministic omitted-rho bounds on both ordinary and log scales.
#'   When a positive bound is smaller than the floating-point range, its
#'   ordinary-scale representation is rounded upward to the smallest positive
#'   double rather than silently reported as zero.
#'
#' @details
#' Uses the hierarchical representation:
#' \enumerate{
#'   \item alpha ~ Gamma(a, b)
#'   \item v_h | alpha ~ Beta(1, alpha) independently
#'   \item w_1 = v_1, w_h = v_h * prod(1 - v_l) for l < h
#'   \item Generate at least \code{n_sticks} terms and continue until the
#'     unallocated remainder \eqn{R} obeys \eqn{R^2 <= remainder_tol}
#'   \item Return the partial sum; its omitted contribution is in
#'     \eqn{[0,R^2]}
#' }
#'
#' Thus finite-stick error is controlled draw by draw rather than being
#' silently discarded. Useful for Monte Carlo validation of analytical
#' formulas.
#'
#' @examples
#' set.seed(42)
#' rho_samples <- rrho(1000, a = 2, b = 1)
#' mean(rho_samples)
#' mean_rho(a = 2, b = 1)
#'
#' @seealso \code{\link{rw1}} for w1 random generation
#'
#' @export
rrho <- function(n, a, b, n_sticks = 500L, remainder_tol = 1e-10,
                 max_sticks = 100000L, strict = TRUE) {
  assert_positive(a, "a")
  assert_positive(b, "b")
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
    .subclass = "dpprior_rho_sampling_error"
  )
  n_sticks <- .dpprior_validate_count(
    n_sticks, "n_sticks", minimum = 10L, maximum = .RRHO_MAX_STICKS,
    .subclass = "dpprior_rho_sampling_error"
  )
  remainder_tol <- .dpprior_validate_scalar(
    remainder_tol, "remainder_tol", lower = 0, upper = 1,
    lower_open = TRUE,
    .subclass = "dpprior_rho_sampling_error"
  )
  max_sticks <- .dpprior_validate_count(
    max_sticks, "max_sticks", minimum = n_sticks,
    maximum = .RRHO_MAX_STICKS,
    .subclass = "dpprior_rho_sampling_error"
  )
  strict <- .dpprior_validate_control(strict, "strict", type = "logical")

  rho_samples <- numeric(n)
  sticks_used <- integer(n)
  remainder_bounds <- numeric(n)
  log_remainder_bounds <- numeric(n)
  failed <- logical(n)
  log_remainder_tol <- log(remainder_tol)
  alpha_draws <- stats::rgamma(n, shape = a, rate = b)
  if (any(!is.finite(alpha_draws)) || any(alpha_draws <= 0)) {
    stop(.dpprior_new_condition(
      "rrho generated a non-finite or non-positive alpha draw",
      classes = c(
        "dpprior_rho_rng_error", "dpprior_numerical_error",
        "dpprior_error", "error"
      ),
      alpha = alpha_draws,
      gamma_shape = a,
      gamma_rate = b,
      code = "invalid_alpha_draw"
    ))
  }

  for (i in seq_len(n)) {
    alpha <- alpha_draws[[i]]
    log_remaining <- 0.0
    rho_value <- 0.0
    used <- 0L

    repeat {
      # The initial batch preserves the historical minimum truncation depth;
      # subsequent batches avoid allocating all max_sticks at once.
      batch_size <- min(
        if (used < n_sticks) n_sticks - used else max(32L, n_sticks),
        max_sticks - used,
        .RRHO_BATCH_SIZE
      )
      if (batch_size <= 0L) break

      v <- stats::rbeta(batch_size, shape1 = 1, shape2 = alpha)
      if (any(!is.finite(v)) || any(v < 0) || any(v > 1)) {
        stop(.dpprior_new_condition(
          sprintf("rrho generated an invalid Beta stick at draw %d", i),
          classes = c(
            "dpprior_rho_rng_error", "dpprior_numerical_error",
            "dpprior_error", "error"
          ),
          sample_index = i,
          alpha = alpha,
          sticks = v,
          code = "invalid_beta_draw"
        ))
      }

      # Work in log remainder space so high-alpha draws can require thousands
      # of sticks without accumulating a long product in an R-level loop.
      cumulative_log_remaining <- log_remaining + cumsum(log1p(-v))
      eligible <- used + seq_len(batch_size) >= n_sticks &
        2 * cumulative_log_remaining <= log_remainder_tol
      stop_index <- if (any(eligible)) which(eligible)[[1L]] else batch_size
      keep <- seq_len(stop_index)
      log_remaining_before <- c(
        log_remaining,
        cumulative_log_remaining[-batch_size]
      )[keep]
      weights <- exp(log_remaining_before) * v[keep]
      rho_value <- rho_value + sum(weights * weights)
      log_remaining <- cumulative_log_remaining[[stop_index]]
      used <- used + stop_index

      if (used >= n_sticks && 2 * log_remaining <= log_remainder_tol) break
      if (used >= max_sticks) break
    }

    rho_samples[[i]] <- rho_value
    sticks_used[[i]] <- used
    log_remainder_bounds[[i]] <- 2 * log_remaining
    remainder_bounds[[i]] <- .rrho_bound_from_log(
      log_remainder_bounds[[i]]
    )
    failed[[i]] <- log_remainder_bounds[[i]] > log_remainder_tol

    if (failed[[i]] && strict) {
      diagnostics <- list(
        schema_version = 1L,
        status = "failed",
        usable = FALSE,
        verified = FALSE,
        method = "adaptive_GEM_remainder_bound",
        sample_index = i,
        alpha = alpha,
        remainder_bound = remainder_bounds[[i]],
        log_remainder_bound = log_remainder_bounds[[i]],
        remainder_tolerance = remainder_tol,
        sticks_used = used,
        max_sticks = max_sticks
      )
      stop(.dpprior_new_condition(
        sprintf(
          paste(
            "rrho truncation failed at draw %d: remainder bound %.6g",
            "exceeds tolerance %.6g after %d sticks"
          ),
          i, remainder_bounds[[i]], remainder_tol, used
        ),
        classes = c(
          "dpprior_rho_truncation_error", "dpprior_numerical_error",
          "dpprior_error", "error"
        ),
        sample_index = i,
        alpha = alpha,
        remainder_bound = remainder_bounds[[i]],
        log_remainder_bound = log_remainder_bounds[[i]],
        tolerance = remainder_tol,
        sticks_used = used,
        max_sticks = max_sticks,
        result = diagnostics
      ))
    }
  }

  diagnostics <- list(
    schema_version = 1L,
    status = if (any(failed)) "approximate" else "converged",
    usable = !any(failed),
    verified = !any(failed),
    message = if (any(failed)) {
      sprintf(
        "%d of %d draws exceeded the deterministic remainder tolerance",
        sum(failed), n
      )
    } else {
      "Every draw met the deterministic squared-remainder tolerance"
    },
    method = "adaptive_GEM_remainder_bound",
    n = n,
    gamma_shape = a,
    gamma_rate = b,
    minimum_sticks = n_sticks,
    max_sticks = max_sticks,
    remainder_tolerance = remainder_tol,
    sticks_used = sticks_used,
    remainder_bounds = remainder_bounds,
    log_remainder_bounds = log_remainder_bounds,
    max_remainder_bound = max(remainder_bounds),
    max_log_remainder_bound = max(log_remainder_bounds),
    failed_indices = which(failed)
  )
  attr(rho_samples, "rrho_diagnostics") <- diagnostics
  rho_samples
}


# =============================================================================
# Grid Computation (for visualization)
# =============================================================================

#' Compute rho Conditional Moments on alpha Grid
#'
#' Computes conditional mean and variance of rho on a grid of alpha values.
#' Useful for visualization of how rho varies with alpha.
#'
#' @param alpha_grid Numeric vector; grid of alpha values.
#'   Default is \code{seq(0.1, 10, length.out = 100)}.
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{alpha}{Grid points}
#'     \item{mean}{E(rho | alpha)}
#'     \item{var}{Var(rho | alpha)}
#'     \item{sd}{SD(rho | alpha)}
#'   }
#'
#' @examples
#' df <- rho_conditional_grid()
#' plot(df$alpha, df$mean, type = "l",
#'      xlab = expression(alpha), ylab = expression(E(rho)))
#'
#' @export
rho_conditional_grid <- function(alpha_grid = seq(0.1, 10, length.out = 100)) {
  assert_positive(alpha_grid, "alpha_grid")

  data.frame(
    alpha = alpha_grid,
    mean = mean_rho_given_alpha(alpha_grid),
    var = var_rho_given_alpha(alpha_grid),
    sd = sqrt(var_rho_given_alpha(alpha_grid))
  )
}

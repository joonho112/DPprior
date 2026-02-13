# =============================================================================
# Module 10: A1 Closed-Form Prior Elicitation (Revised)
# =============================================================================
#
# This module provides:
# 1. Closed-form mapping from (mu_K, var_K, J) to Gamma(a, b) hyperprior
# 2. VIF (Variance Inflation Factor) utilities
# 3. Confidence level to VIF conversion
# 4. S3 class DPprior_fit for results
#
# Author: JoonHo Lee
# Date: December 2025
# Part of: DPprior R Package
# Reference: Lee (2026) arXiv:2602.06301, Section 3.1 (TSMM Stage 1)
#
# Revision Notes:
# - Added epsilon validation for robustness
# - Added converged/iterations fields for A2 compatibility
# - Fixed digamma scaling to use numeric floor
# - Strengthened VIF validation to require >= 1
# =============================================================================


# =============================================================================
# Scaling Constant Computation
# =============================================================================

#' Compute Scaling Constant for A1 Mapping
#'
#' Computes the scaling constant \eqn{c_J} used in the A1 closed-form mapping.
#' Three variants are supported based on asymptotic approximations.
#'
#' @param J Integer; number of items/sites (must be >= 2).
#' @param scaling Character; one of "log", "harmonic", or "digamma".
#' @param mu_K Numeric; target mean of K (required for "digamma" scaling).
#'
#' @return Numeric scalar; the scaling constant \eqn{c_J}.
#'
#' @details
#' The scaling constant appears in the Poisson proxy:
#' \deqn{K_J - 1 \mid \alpha \approx \text{Poisson}(\alpha \cdot c_J)}
#'
#' Available variants:
#' \describe{
#'   \item{log}{\eqn{c_J = \log(J)}, the asymptotic leading term (default)}
#'   \item{harmonic}{\eqn{c_J = H_{J-1} = \psi(J) + \gamma}, improves
#'     accuracy for small/moderate J}
#'   \item{digamma}{\eqn{c_J = \psi(\tilde{\alpha} + J) - \psi(\tilde{\alpha})}
#'     where \eqn{\tilde{\alpha} = (\mu_K - 1)/\log(J)}, a local correction}
#' }
#'
#' @seealso \code{\link{DPprior_a1}} for the main elicitation function
#'
#' @examples
#' # Default log scaling
#' compute_scaling_constant(50, "log")
#'
#' # Harmonic scaling (better for moderate J)
#' compute_scaling_constant(50, "harmonic")
#'
#' # Digamma scaling (requires mu_K)
#' compute_scaling_constant(50, "digamma", mu_K = 5)
#'
#' @keywords internal
#' @export
compute_scaling_constant <- function(J, scaling = c("log", "harmonic", "digamma"),
                                     mu_K = NULL) {
  scaling <- match.arg(scaling)

  switch(scaling,
         log = log(J),
         harmonic = digamma(J) + .EULER_GAMMA,
         digamma = {
           if (is.null(mu_K)) {
             stop("mu_K required for digamma scaling", call. = FALSE)
           }
           alpha_tilde <- (mu_K - 1) / log(J)
           # Use numeric floor instead of fallback (robust to edge cases)
           alpha_tilde <- max(alpha_tilde, .Machine$double.eps)
           digamma(alpha_tilde + J) - digamma(alpha_tilde)
         }
  )
}


# =============================================================================
# Main A1 Elicitation Function
# =============================================================================

#' A1 Closed-Form Prior Elicitation
#'
#' Maps target beliefs about the number of clusters \eqn{(\mu_K, \sigma^2_K)}
#' to Gamma hyperprior parameters \eqn{(a, b)} using the A1 closed-form
#' approximation based on Negative Binomial moment matching.
#'
#' @param J Integer; number of items/sites (must be >= 2).
#' @param mu_K Numeric; target prior mean of \eqn{K_J} (must be > 1 and <= J).
#' @param var_K Numeric; target prior variance of \eqn{K_J} (must be > 0).
#' @param scaling Character; scaling constant method: "log" (default),
#'   "harmonic", or "digamma".
#' @param epsilon Numeric; buffer for feasibility projection. Default is
#'   \code{.TOL_PROJECTION_BUFFER} (1e-6).
#'
#' @return An S3 object of class \code{DPprior_fit} with components:
#'   \describe{
#'     \item{a}{Shape parameter of the Gamma prior}
#'     \item{b}{Rate parameter of the Gamma prior}
#'     \item{J}{Sample size used}
#'     \item{target}{List with target moments and type}
#'     \item{method}{"A1" indicating closed-form method}
#'     \item{status}{"success" or "projected" if boundary adjustment needed}
#'     \item{scaling}{Scaling method used}
#'     \item{cJ}{Scaling constant value}
#'     \item{var_K_used}{Actual variance used (may differ if projected)}
#'     \item{converged}{Always TRUE for A1 (for A2 compatibility)}
#'     \item{iterations}{Always 0L for A1 (for A2 compatibility)}
#'     \item{fit}{NULL (placeholder for A2 refinement)}
#'     \item{diagnostics}{NULL (placeholder for diagnostics)}
#'     \item{trace}{NULL (placeholder for optimization trace)}
#'   }
#'
#' @details
#' ## Theory (TSMM Stage 1)
#'
#' The A1 method uses a shifted Negative Binomial approximation:
#' \deqn{K_J - 1 \mid \alpha \approx \text{Poisson}(\alpha \cdot c_J)}
#'
#' With \eqn{\alpha \sim \text{Gamma}(a, b)}, the marginal becomes:
#' \deqn{K_J - 1 \approx \text{NegBin}(a, b/(b + c_J))}
#'
#' ## Inverse Formulas (Theorem 1)
#'
#' Let \eqn{m = \mu_K - 1} (shifted mean) and \eqn{D = \sigma^2_K - m}.
#' If \eqn{D > 0} (overdispersion):
#' \deqn{a = m^2 / D, \quad b = m \cdot c_J / D}
#'
#' If \eqn{D \leq 0} (infeasible), the variance is projected to the boundary.
#'
#' ## Feasibility
#'
#' The NegBin model requires overdispersion: \eqn{\sigma^2_K > \mu_K - 1}.
#' High-confidence specifications (low variance) may violate this constraint
#' under the A1 proxy, even though they may be feasible under the exact DP.
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @seealso
#' \code{\link{vif_to_variance}} for VIF conversion,
#' \code{\link{confidence_to_vif}} for confidence mapping,
#' \code{\link{print.DPprior_fit}} for print method
#'
#' @examples
#' # Basic usage with moment targets
#' fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)
#' print(fit)
#'
#' # Using VIF specification
#' fit <- DPprior_a1(J = 50, mu_K = 5, var_K = vif_to_variance(5, 2))
#'
#' # Using confidence-based specification
#' vif <- confidence_to_vif("medium")
#' fit <- DPprior_a1(J = 50, mu_K = 5, var_K = vif_to_variance(5, vif))
#'
#' # Infeasible variance (triggers projection)
#' fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 3)  # var < mu-1 = 4
#'
#' # Compare scaling methods
#' fit_log <- DPprior_a1(50, 5, 8, scaling = "log")
#' fit_harm <- DPprior_a1(50, 5, 8, scaling = "harmonic")
#'
#' @family elicitation
#'
#' @export
DPprior_a1 <- function(J, mu_K, var_K,
                       scaling = c("log", "harmonic", "digamma"),
                       epsilon = .TOL_PROJECTION_BUFFER) {

  scaling <- match.arg(scaling)

  # Input validation
  if (!is.numeric(J) || length(J) != 1L || !is.finite(J) ||
      J != floor(J) || J < 2L) {
    stop("J must be an integer >= 2", call. = FALSE)
  }
  if (!is.numeric(mu_K) || length(mu_K) != 1L || !is.finite(mu_K)) {
    stop("mu_K must be a finite numeric scalar", call. = FALSE)
  }
  if (mu_K <= 1) {
    stop("mu_K must be > 1 (at least one cluster is always present)", call. = FALSE)
  }
  if (mu_K > J) {
    stop("mu_K must be <= J (cannot have more clusters than items)", call. = FALSE)
  }
  if (!is.numeric(var_K) || length(var_K) != 1L || !is.finite(var_K) ||
      var_K <= 0) {
    stop("var_K must be a positive finite numeric scalar", call. = FALSE)
  }
  # Epsilon validation
  if (!is.numeric(epsilon) || length(epsilon) != 1L ||
      !is.finite(epsilon) || epsilon <= 0) {
    stop("epsilon must be a finite numeric scalar > 0", call. = FALSE)
  }

  # Compute scaling constant
  cJ <- compute_scaling_constant(J, scaling, mu_K)

  # Shifted mean (under A1 with shift s = 1)
  mu_S <- mu_K - 1

  # Compute denominator D = var_K - mu_S
  # For NegBin to be valid, we need var_K > mu_S (overdispersion)
  denom <- var_K - mu_S
  var_K_used <- var_K

  # Feasibility check with scaled projection buffer
  # Use epsilon * (1 + mu_S^2) for scale-invariant buffering
  epsilon_scaled <- epsilon * (1 + mu_S^2)

  if (denom <= epsilon_scaled) {
    # Project to feasible boundary
    var_K_used <- mu_S + epsilon_scaled
    denom <- epsilon_scaled
    status <- "projected"
    warning("var_K <= mu_K - 1: projected to feasible boundary", call. = FALSE)
  } else {
    status <- "success"
  }

  # Closed-form inverse (Theorem 1)
  a0 <- mu_S^2 / denom
  b0 <- mu_S * cJ / denom

  # Additional diagnostics for edge cases
  if (a0 < 0.01) {
    status <- paste(status, "[warning: quasi-improper prior, a < 0.01]")
  }
  if (a0 > 1e6) {
    status <- paste(status, "[warning: quasi-degenerate prior, a > 1e6]")
  }

  # Construct DPprior_fit object
  structure(
    list(
      a = a0,
      b = b0,
      J = J,
      target = list(
        mu_K = mu_K,
        var_K = var_K,
        type = "moments"
      ),
      method = "A1",
      status = status,
      scaling = scaling,
      cJ = cJ,
      var_K_used = var_K_used,
      converged = TRUE,      # A2 compatibility
      iterations = 0L,       # A2 compatibility
      fit = NULL,
      diagnostics = NULL,
      trace = NULL
    ),
    class = "DPprior_fit"
  )
}


# =============================================================================
# VIF (Variance Inflation Factor) Utilities
# =============================================================================

#' Convert Variance Inflation Factor to Variance
#'
#' Converts a Variance Inflation Factor (VIF) specification to the actual
#' variance of \eqn{K_J}.
#'
#' @param mu_K Numeric; target prior mean of \eqn{K_J}.
#' @param vif Numeric; Variance Inflation Factor (must be >= 1 for A1 feasibility).
#'
#' @return Numeric; variance of \eqn{K_J} computed as \eqn{(\mu_K - 1) \times \text{VIF}}.
#'
#' @details
#' The VIF is defined as:
#' \deqn{\text{VIF} = \frac{\sigma^2_K}{\mu_K - 1}}
#'
#' Interpretation:
#' \describe{
#'   \item{VIF = 1}{Poisson variance (exact boundary for A1)}
#'   \item{VIF > 1}{Overdispersion (required for A1 feasibility)}
#'   \item{VIF < 1}{Underdispersion (infeasible for A1, not allowed)}
#' }
#'
#' @seealso
#' \code{\link{confidence_to_vif}} for mapping confidence levels to VIF,
#' \code{\link{cv_alpha_to_variance}} for CV-based specification
#'
#' @examples
#' # VIF = 2 means variance is twice the Poisson variance
#' vif_to_variance(mu_K = 5, vif = 2)  # Returns 8
#'
#' # Use with DPprior_a1
#' fit <- DPprior_a1(J = 50, mu_K = 5, var_K = vif_to_variance(5, 2))
#'
#' @export
vif_to_variance <- function(mu_K, vif) {
  if (!is.numeric(mu_K) || !is.numeric(vif)) {
    stop("mu_K and vif must be numeric", call. = FALSE)
  }
  if (any(mu_K <= 1)) {
    stop("mu_K must be > 1", call. = FALSE)
  }
  # VIF >= 1 required for A1 feasibility
  if (any(vif < 1)) {
    stop("vif must be >= 1 (values < 1 imply infeasible underdispersion)",
         call. = FALSE)
  }

  (mu_K - 1) * vif
}


#' Map a Qualitative Confidence Level to a Variance Inflation Factor (VIF)
#'
#' Maps intuitive confidence levels to VIF values for easy prior specification.
#'
#' @param confidence Character; one of "low", "medium", or "high".
#'
#' @return Numeric; VIF value (5.0 for low, 2.5 for medium, 1.5 for high).
#'
#' @details
#' The mapping is:
#' \describe{
#'   \item{low}{VIF = 5.0; high uncertainty about \eqn{K_J}}
#'   \item{medium}{VIF = 2.5; moderate uncertainty}
#'   \item{high}{VIF = 1.5; high confidence (near Poisson boundary)}
#' }
#'
#' Higher confidence implies lower variance, which corresponds to lower VIF.
#' The "high" setting (VIF = 1.5) is close to the A1 feasibility boundary.
#'
#' @seealso \code{\link{vif_to_variance}} for converting VIF to variance
#'
#' @examples
#' # Get VIF for medium confidence
#' vif <- confidence_to_vif("medium")  # Returns 2.5
#'
#' # Complete workflow
#' mu_K <- 5
#' vif <- confidence_to_vif("low")
#' var_K <- vif_to_variance(mu_K, vif)
#' fit <- DPprior_a1(J = 50, mu_K = mu_K, var_K = var_K)
#'
#' @export
confidence_to_vif <- function(confidence = c("low", "medium", "high")) {
  # Use match.arg for standard R idiom
  confidence <- match.arg(confidence)

  switch(
    confidence,
    low = 5.0,
    medium = 2.5,
    high = 1.5
  )
}


#' Convert CV(alpha) to Variance
#'
#' Converts a coefficient of variation specification for \eqn{\alpha} to
#' the implied variance of \eqn{K_J} under the A1 approximation.
#'
#' @param mu_K Numeric; target prior mean of \eqn{K_J}.
#' @param cv_alpha Numeric; target coefficient of variation for \eqn{\alpha}.
#'
#' @return Numeric; implied variance of \eqn{K_J}.
#'
#' @details
#' Under the A1 approximation:
#' \deqn{\text{CV}(\alpha) = 1/\sqrt{a} = \frac{\sqrt{\sigma^2_K - m}}{m}}
#'
#' where \eqn{m = \mu_K - 1}. Inverting:
#' \deqn{\sigma^2_K = m + (\text{CV}(\alpha) \cdot m)^2 = m(1 + \text{CV}(\alpha)^2 \cdot m)}
#'
#' @seealso \code{\link{vif_to_variance}} for VIF-based specification
#'
#' @examples
#' # CV(alpha) = 0.5 means moderate prior concentration
#' var_K <- cv_alpha_to_variance(mu_K = 5, cv_alpha = 0.5)
#'
#' # Verify round-trip
#' fit <- DPprior_a1(J = 50, mu_K = 5, var_K = var_K)
#' 1 / sqrt(fit$a)  # Should be approximately 0.5
#'
#' @export
cv_alpha_to_variance <- function(mu_K, cv_alpha) {
  if (!is.numeric(mu_K) || !is.numeric(cv_alpha)) {
    stop("mu_K and cv_alpha must be numeric", call. = FALSE)
  }
  if (any(mu_K <= 1)) {
    stop("mu_K must be > 1", call. = FALSE)
  }
  if (any(cv_alpha <= 0)) {
    stop("cv_alpha must be positive", call. = FALSE)
  }

  m <- mu_K - 1
  m + (cv_alpha * m)^2
}


# =============================================================================
# S3 Methods for DPprior_fit
# =============================================================================
# Note: summary.DPprior_fit is defined in R/17_s3_methods.R (canonical location)

#' Coerce DPprior_fit to Data Frame
#'
#' @param x A \code{DPprior_fit} object.
#' @param row.names Optional row names.
#' @param optional Logical; if TRUE, avoid setting names.
#' @param ... Additional arguments (ignored).
#'
#' @return A data frame with one row containing the fit results.
#'
#' @examples
#' fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)
#' as.data.frame(fit)
#'
#' @export
as.data.frame.DPprior_fit <- function(x, row.names = NULL, optional = FALSE, ...) {
  data.frame(
    method = x$method,
    status = x$status,
    a = x$a,
    b = x$b,
    J = x$J,
    mu_K = x$target$mu_K,
    var_K = x$target$var_K,
    mean_alpha = x$a / x$b,
    cv_alpha = 1 / sqrt(x$a),
    scaling = if (!is.null(x$scaling)) x$scaling else NA_character_,
    converged = x$converged,
    iterations = x$iterations,
    stringsAsFactors = FALSE,
    row.names = row.names
  )
}


# =============================================================================
# Verification Functions
# =============================================================================

#' Verify A1 Mapping via Round-Trip
#'
#' Tests the A1 mapping by computing the forward model (NegBin moments)
#' from the derived \eqn{(a, b)} parameters and comparing to targets.
#'
#' @param fit A \code{DPprior_fit} object from \code{DPprior_a1}.
#' @param tol Numeric; tolerance for relative error comparison.
#' @param verbose Logical; if TRUE, print verification details.
#'
#' @return Logical; TRUE if round-trip succeeds within tolerance.
#'
#' @details
#' Under the A1 NegBin approximation:
#' \deqn{K_J - 1 \sim \text{NegBin}(a, p)}
#' where \eqn{p = b/(b + c_J)}.
#'
#' @examples
#' \dontrun{
#' fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)
#' verify_a1_roundtrip(fit)
#'
#' }
#' @keywords internal
verify_a1_roundtrip <- function(fit, tol = 1e-8, verbose = TRUE) {
  if (!inherits(fit, "DPprior_fit")) {
    stop("fit must be a DPprior_fit object", call. = FALSE)
  }
  if (fit$method != "A1") {
    warning("Round-trip verification is specific to A1 method", call. = FALSE)
  }

  # Extract parameters
  a <- fit$a
  b <- fit$b
  cJ <- fit$cJ

  # Forward model: NegBin moments
  # p = b / (b + cJ)
  # Mean(K-1) = a * (1-p) / p = a * cJ / b
  # Var(K-1) = a * (1-p) / p^2 = a * cJ * (b + cJ) / b^2
  p <- b / (b + cJ)
  mean_shifted <- a * (1 - p) / p
  var_shifted <- a * (1 - p) / p^2

  mu_K_recovered <- mean_shifted + 1
  var_K_recovered <- var_shifted

  # Compare to targets (use var_K_used for projected cases)
  mu_K_target <- fit$target$mu_K
  var_K_target <- if (!is.null(fit$var_K_used)) fit$var_K_used else fit$target$var_K

  rel_err_mu <- abs(mu_K_recovered - mu_K_target) / mu_K_target
  rel_err_var <- abs(var_K_recovered - var_K_target) / var_K_target

  passed <- (rel_err_mu < tol) && (rel_err_var < tol)

  if (isTRUE(verbose)) {
    cat("A1 Round-Trip Verification\n")
    cat(paste0(rep("-", 40), collapse = ""), "\n")
    cat(sprintf("mu_K: target = %.6f, recovered = %.6f, rel_err = %.2e\n",
                mu_K_target, mu_K_recovered, rel_err_mu))
    cat(sprintf("var_K: target = %.6f, recovered = %.6f, rel_err = %.2e\n",
                var_K_target, var_K_recovered, rel_err_var))
    cat(sprintf("Result: %s\n", if (passed) "PASS" else "FAIL"))
  }

  invisible(passed)
}


#' Run All Module 10 Verification Tests
#'
#' Comprehensive verification suite for the A1 closed-form mapping module.
#'
#' @param verbose Logical; if TRUE, print detailed results.
#'
#' @return Logical; TRUE if all tests pass.
#'
#' @examples
#' \dontrun{
#' verify_a1_mapping_all()
#'
#' }
#' @keywords internal
verify_a1_mapping_all <- function(verbose = TRUE) {
  if (isTRUE(verbose)) {
    cat("=", rep("=", 69), "\n", sep = "")
    cat("Module 10: A1 Closed-Form Mapping - Full Verification Suite\n")
    cat("=", rep("=", 69), "\n\n", sep = "")
  }

  all_pass <- TRUE

  # Test 1: Basic positive parameters
  if (isTRUE(verbose)) {
    cat("[Test 1] A1 produces positive parameters\n")
  }
  fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)
  test1 <- (fit$a > 0) && (fit$b > 0) && (fit$method == "A1")
  if (isTRUE(verbose)) {
    cat(sprintf("  J=50, mu_K=5, var_K=8\n"))
    cat(sprintf("  a = %.6f, b = %.6f\n", fit$a, fit$b))
    cat(sprintf("  PASS: %s\n\n", test1))
  }
  all_pass <- all_pass && test1

  # Test 2: Infeasible variance (projection)
  if (isTRUE(verbose)) {
    cat("[Test 2] A1 handles infeasible variance via projection\n")
  }
  fit2 <- suppressWarnings(DPprior_a1(J = 50, mu_K = 5, var_K = 3))
  test2 <- grepl("projected", fit2$status)
  if (isTRUE(verbose)) {
    cat(sprintf("  J=50, mu_K=5, var_K=3 (infeasible: var < mu-1=4)\n"))
    cat(sprintf("  var_K_used = %.6f\n", fit2$var_K_used))
    cat(sprintf("  status = %s\n", fit2$status))
    cat(sprintf("  PASS: %s\n\n", test2))
  }
  all_pass <- all_pass && test2

  # Test 3: VIF conversion
  if (isTRUE(verbose)) {
    cat("[Test 3] VIF conversion is correct\n")
  }
  var_computed <- vif_to_variance(5, 2)
  test3a <- abs(var_computed - 8) < 1e-10
  vif_medium <- confidence_to_vif("medium")
  test3b <- abs(vif_medium - 2.5) < 1e-10
  if (isTRUE(verbose)) {
    cat(sprintf("  vif_to_variance(5, 2) = %.1f (expected 8)\n", var_computed))
    cat(sprintf("  confidence_to_vif('medium') = %.1f (expected 2.5)\n", vif_medium))
    cat(sprintf("  PASS: %s\n\n", test3a && test3b))
  }
  all_pass <- all_pass && test3a && test3b

  # Test 4: Round-trip verification
  if (isTRUE(verbose)) {
    cat("[Test 4] Round-trip verification\n")
  }
  test_cases <- list(
    list(J = 50, mu_K = 5.0, var_K = 8.0),
    list(J = 100, mu_K = 10.0, var_K = 20.0),
    list(J = 25, mu_K = 4.0, var_K = 30.0),
    list(J = 50, mu_K = 5.0, var_K = 6.0)
  )

  test4_pass <- TRUE
  for (tc in test_cases) {
    fit <- DPprior_a1(J = tc$J, mu_K = tc$mu_K, var_K = tc$var_K)
    passed <- verify_a1_roundtrip(fit, verbose = FALSE)
    test4_pass <- test4_pass && passed
    if (isTRUE(verbose)) {
      cat(sprintf("  J=%d, mu_K=%.1f, var_K=%.1f: %s\n",
                  tc$J, tc$mu_K, tc$var_K, if (passed) "PASS" else "FAIL"))
    }
  }
  if (isTRUE(verbose)) cat("\n")
  all_pass <- all_pass && test4_pass

  # Test 5: Scaling method comparison
  if (isTRUE(verbose)) {
    cat("[Test 5] Scaling method comparison\n")
  }
  for (scaling in c("log", "harmonic", "digamma")) {
    fit <- DPprior_a1(50, 5, 8, scaling = scaling)
    if (isTRUE(verbose)) {
      cat(sprintf("  scaling='%s': cJ=%.4f, a=%.4f, b=%.4f\n",
                  scaling, fit$cJ, fit$a, fit$b))
    }
  }
  if (isTRUE(verbose)) cat("\n")

  # Test 6: CV to variance conversion
  if (isTRUE(verbose)) {
    cat("[Test 6] CV(alpha) to variance conversion\n")
  }
  test6_pass <- TRUE
  for (cv_target in c(0.5, 1.0, 2.0)) {
    var_K <- cv_alpha_to_variance(5, cv_target)
    fit <- DPprior_a1(50, 5, var_K)
    cv_recovered <- 1 / sqrt(fit$a)
    rel_err <- abs(cv_recovered - cv_target) / cv_target
    passed <- rel_err < 0.01
    test6_pass <- test6_pass && passed
    if (isTRUE(verbose)) {
      cat(sprintf("  Target CV=%.1f -> var_K=%.2f -> Recovered CV=%.4f\n",
                  cv_target, var_K, cv_recovered))
    }
  }
  if (isTRUE(verbose)) cat("\n")
  all_pass <- all_pass && test6_pass

  # Test 7: A2 compatibility fields
  if (isTRUE(verbose)) {
    cat("[Test 7] A2 compatibility fields\n")
  }
  fit <- DPprior_a1(50, 5, 8)
  test7 <- isTRUE(fit$converged) && identical(fit$iterations, 0L)
  if (isTRUE(verbose)) {
    cat(sprintf("  converged = %s (expected TRUE)\n", fit$converged))
    cat(sprintf("  iterations = %d (expected 0)\n", fit$iterations))
    cat(sprintf("  PASS: %s\n\n", test7))
  }
  all_pass <- all_pass && test7

  # Summary
  if (isTRUE(verbose)) {
    cat("=", rep("=", 69), "\n", sep = "")
    cat(sprintf("Overall Result: %s\n",
                if (all_pass) "ALL TESTS PASSED" else "SOME TESTS FAILED"))
    cat("=", rep("=", 69), "\n", sep = "")
  }

  invisible(all_pass)
}

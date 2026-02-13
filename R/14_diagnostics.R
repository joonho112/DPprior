# =============================================================================
# Module 14: Comprehensive Diagnostics for DPprior
# =============================================================================
#
# This module implements prior validation and dominance risk assessment
# following the dual-anchor diagnostic framework (Lee, 2026, Section 4).
#
# Key insight (Lee, 2026, Section 4): matching K_J does NOT automatically calibrate weight
# behavior. A prior that achieves excellent fit to target K_J may still induce:
# - High probability of weight dominance (P(w1 > 0.5) >> expected)
# - Extreme co-clustering probability
# - Unintended "degeneracy" toward single-cluster solutions
#
# Author: JoonHo Lee
# Date: December 2025
# Part of: DPprior R Package
# References: Lee (2026), Section 4
# =============================================================================


# =============================================================================
# Alpha Distribution Diagnostics
# =============================================================================

#' Alpha Distribution Diagnostics
#'
#' Computes summary statistics for alpha distributed as Gamma(a, b).
#'
#' @param a Numeric; shape parameter of the Gamma prior (> 0).
#' @param b Numeric; rate parameter of the Gamma prior (> 0).
#'
#' @return A list with components: mean (expected value E[alpha] = a/b),
#'   sd (standard deviation), cv (coefficient of variation = 1/sqrt(a)),
#'   median, and quantiles (named vector with q5, q25, q50, q75, q95).
#'
#' @details
#' The coefficient of variation depends only on the shape parameter a,
#' making it a useful summary of prior uncertainty regardless of the mean.
#'
#' @keywords internal
compute_alpha_diagnostics <- function(a, b) {
  assert_positive(a, "a")
  assert_positive(b, "b")

  mean_alpha <- a / b
  var_alpha <- a / b^2
  sd_alpha <- sqrt(var_alpha)
  cv_alpha <- 1 / sqrt(a)

  probs <- c(0.05, 0.25, 0.50, 0.75, 0.95)
  quantiles <- qgamma(probs, shape = a, rate = b)
  names(quantiles) <- paste0("q", formatC(100 * probs, format = "f", digits = 0))

  list(
    mean = mean_alpha,
    sd = sd_alpha,
    cv = cv_alpha,
    median = unname(quantiles["q50"]),
    quantiles = quantiles
  )
}


# =============================================================================
# K_J PMF Helper (with fallback logic)
# =============================================================================

#' Get K_J PMF via pmf_K_marginal
#'
#' Internal helper that computes the marginal PMF of K_J using
#' pmf_K_marginal() and Stirling numbers.
#'
#' @param J Integer; sample size.
#' @param a,b Numeric; Gamma hyperparameters.
#' @param M Integer; quadrature nodes.
#'
#' @return List with pmf (numeric vector) and support (integer vector).
#'
#' @keywords internal
.get_K_pmf_support <- function(J, a, b, M = .QUAD_NODES_DEFAULT) {
  assert_valid_J(J)
  assert_positive(a, "a")
  assert_positive(b, "b")

  J <- as.integer(J)

  if (!exists("pmf_K_marginal", mode = "function")) {
    stop("pmf_K_marginal() is not available.", call. = FALSE)
  }
  if (!exists("compute_log_stirling", mode = "function")) {
    stop("compute_log_stirling() is required to compute K_J PMF.", call. = FALSE)
  }

  logS <- compute_log_stirling(J)
  pmf0 <- pmf_K_marginal(J, a, b, logS, M)

  # pmf_K_marginal returns length J+1 for k = 0..J
  if (length(pmf0) != J + 1L) {
    stop("pmf_K_marginal() must return length J+1 for k=0..J.", call. = FALSE)
  }

  pmf <- as.numeric(pmf0[-1L])  # Drop P(K=0)
  pmf_sum <- sum(pmf)
  if (!is.finite(pmf_sum) || pmf_sum <= 0) {
    warning("PMF normalization failed: sum is zero or non-finite. Using uniform PMF.", call. = FALSE)
    pmf <- rep(1 / length(pmf), length(pmf))
  } else {
    pmf <- pmf / pmf_sum
  }

  list(pmf = pmf, support = seq_len(J))
}


# =============================================================================
# K Distribution Diagnostics
# =============================================================================

#' K Distribution Diagnostics
#'
#' Computes summary statistics for the marginal distribution of K_J
#' under alpha distributed as Gamma(a, b).
#'
#' @param J Integer; sample size (positive integer >= 1).
#' @param a Numeric; shape parameter of the Gamma prior (> 0).
#' @param b Numeric; rate parameter of the Gamma prior (> 0).
#' @param M Integer; number of Gauss-Laguerre quadrature nodes (default: 80).
#'
#' @return A list with components: mean (expected value of K), var (variance),
#'   sd, mode, median, quantiles (named integer vector), and pmf (full PMF
#'   vector for k = 1, ..., J).
#'
#' @keywords internal
compute_K_diagnostics <- function(J, a, b, M = .QUAD_NODES_DEFAULT) {
  assert_valid_J(J)
  assert_positive(a, "a")
  assert_positive(b, "b")

  J <- as.integer(J)

  # Moments via quadrature
  moments <- exact_K_moments(J, a, b, M)

  # PMF via helper with fallback
  pmf_obj <- .get_K_pmf_support(J, a, b, M)
  pmf <- pmf_obj$pmf
  k_support <- pmf_obj$support

  cdf <- cumsum(pmf)

  # Find quantiles
  probs <- c(0.05, 0.25, 0.50, 0.75, 0.95)
  qfun <- function(p) {
    idx <- which(cdf >= p)
    if (length(idx) == 0L) return(as.integer(max(k_support)))
    as.integer(k_support[min(idx)])
  }
  quantiles <- vapply(probs, qfun, integer(1))
  names(quantiles) <- paste0("q", formatC(100 * probs, format = "f", digits = 0))

  # Mode
  mode_K <- as.integer(k_support[which.max(pmf)])

  list(
    mean = as.numeric(moments$mean),
    var = as.numeric(moments$var),
    sd = sqrt(as.numeric(moments$var)),
    mode = mode_K,
    median = as.integer(quantiles["q50"]),
    quantiles = quantiles,
    pmf = pmf
  )
}


# =============================================================================
# Weight Distribution Diagnostics (Lee, 2026, Section 4)
# =============================================================================

#' Weight Distribution Diagnostics (w1)
#'
#' Computes comprehensive diagnostics for the first stick-breaking weight.
#' This is the KEY diagnostic for concerns about unintended prior behavior (Lee, 2026, Section 4).
#'
#' @param a Numeric; shape parameter of the Gamma prior (> 0).
#' @param b Numeric; rate parameter of the Gamma prior (> 0).
#' @param thresholds Numeric vector; thresholds for computing P(w1 > t).
#'   Default is c(0.3, 0.5, 0.7, 0.9). All values must be in (0, 1).
#' @param M Integer; number of quadrature nodes for moment computation (default: 80).
#'
#' @return A list with components: mean (expected value of the first weight),
#'   median, quantiles (named vector), prob_exceeds (named vector of exceedance
#'   probabilities for each threshold), and dominance_risk (character: "low",
#'   "moderate", or "high").
#'
#' @details
#' The weight w1 is the first stick-breaking weight in GEM order
#' (size-biased permutation), representing the asymptotic cluster share of
#' a randomly chosen unit.
#'
#' Dominance risk classification:
#' \itemize{
#'   \item "low": P(w1 > 0.5) < 0.2
#'   \item "moderate": 0.2 <= P(w1 > 0.5) < 0.4
#'   \item "high": P(w1 > 0.5) >= 0.4
#' }
#'
#' @examples
#' # Lee et al. DP-inform prior (known high dominance)
#' compute_weight_diagnostics(1.60, 1.22)
#'
#' # Lower dominance case
#' compute_weight_diagnostics(5, 1)
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @family diagnostics
#'
#' @export
compute_weight_diagnostics <- function(a, b,
                                       thresholds = c(0.3, 0.5, 0.7, 0.9),
                                       M = .QUAD_NODES_DEFAULT) {
  assert_positive(a, "a")
  assert_positive(b, "b")
  assert_probability(thresholds, "thresholds")

  # Mean (via quadrature)
  mean_val <- mean_w1(a, b, M)

  # Quantiles via closed-form inverse CDF
  probs <- c(0.05, 0.25, 0.50, 0.75, 0.95)
  quantiles <- quantile_w1(probs, a, b)
  names(quantiles) <- paste0("q", formatC(100 * probs, format = "f", digits = 0))

  # P(w1 > t) for each threshold
  prob_exceeds <- vapply(thresholds, function(t) prob_w1_exceeds(t, a, b), numeric(1))
  names(prob_exceeds) <- paste0("prob_gt_", thresholds)

  # Dominance risk classification (always based on 0.5 threshold)
  p_dom <- prob_w1_exceeds(0.5, a, b)
  dominance_risk <- if (p_dom < 0.2) {
    "low"
  } else if (p_dom < 0.4) {
    "moderate"
  } else {
    "high"
  }

  list(
    mean = as.numeric(mean_val),
    median = unname(quantiles["q50"]),
    quantiles = quantiles,
    prob_exceeds = prob_exceeds,
    dominance_risk = dominance_risk
  )
}


# =============================================================================
# Co-Clustering Diagnostics
# =============================================================================

#' Co-Clustering Diagnostics (rho)
#'
#' Computes diagnostics for rho = sum of w_h squared, the prior
#' co-clustering probability.
#'
#' @param a Numeric; shape parameter of the Gamma prior (> 0).
#' @param b Numeric; rate parameter of the Gamma prior (> 0).
#' @param M Integer; number of quadrature nodes (default: 80).
#'
#' @return A list with components: mean (expected co-clustering probability),
#'   var (variance), sd, and interpretation (qualitative description of
#'   co-clustering level).
#'
#' @details
#' Key identity from Lee (2026, Section 4): the expected co-clustering probability equals the
#' expected first stick-breaking weight.
#'
#' @keywords internal
compute_coclustering_diagnostics <- function(a, b, M = .QUAD_NODES_DEFAULT) {
  assert_positive(a, "a")
  assert_positive(b, "b")

  mean_val <- mean_rho(a, b, M)
  var_val <- var_rho(a, b, M)

  interpretation <- if (mean_val > 0.5) {
    "High prior co-clustering: most unit pairs expected in same cluster"
  } else if (mean_val > 0.3) {
    "Moderate prior co-clustering"
  } else if (mean_val > 0.1) {
    "Low prior co-clustering: units likely in different clusters"
  } else {
    "Very low prior co-clustering: highly fragmented prior"
  }

  list(
    mean = as.numeric(mean_val),
    var = as.numeric(var_val),
    sd = sqrt(as.numeric(var_val)),
    interpretation = interpretation
  )
}


# =============================================================================
# Quick Dominance Risk Check
# =============================================================================

#' Quick Dominance Risk Check
#'
#' Fast check for high dominance risk without computing full diagnostics.
#'
#' @param a Numeric; shape parameter of the Gamma prior (> 0).
#' @param b Numeric; rate parameter of the Gamma prior (> 0).
#' @param threshold Numeric; dominance threshold in (0, 1) (default: 0.5).
#' @param risk_level Numeric; probability threshold for flagging risk (default: 0.3).
#'
#' @return Logical; TRUE if P(w1 > threshold) > risk_level.
#'
#' @examples
#' check_dominance_risk(1.60, 1.22)
#' check_dominance_risk(5, 1)
#'
#' @export
check_dominance_risk <- function(a, b, threshold = 0.5, risk_level = 0.3) {
  assert_positive(a, "a")
  assert_positive(b, "b")
  assert_probability(threshold, "threshold")

  prob_w1_exceeds(threshold, a, b) > risk_level
}


# =============================================================================
# Comprehensive Prior Diagnostics
# =============================================================================

#' Comprehensive Prior Diagnostics
#'
#' Computes a full diagnostic report for a fitted DPprior object, implementing
#' the "unintended prior" checks from Lee (2026, Section 4).
#'
#' @param fit A DPprior_fit object from any calibration method
#'   (e.g., DPprior_a1, DPprior_a2_newton, DPprior_dual).
#'   Must contain fields a, b, and J. Optionally can contain M for quadrature nodes.
#' @param thresholds Numeric vector; thresholds for weight dominance checks
#'   (default: c(0.5, 0.9)).
#'
#' @return An S3 object of class "DPprior_diagnostics" with components:
#'   J (sample size), a and b (Gamma parameters), alpha (alpha distribution
#'   summary), K (K distribution summary), weights (w1 distribution summary
#'   with dominance risk), coclustering (rho summary), and warnings
#'   (character vector of diagnostic warnings).
#'
#' @details
#' Warnings are issued if:
#' \itemize{
#'   \item P(w1 > 0.5) > 0.4: "HIGH DOMINANCE RISK"
#'   \item P(w1 > 0.9) > 0.15: "NEAR-DEGENERATE RISK"
#' }
#'
#' @examples
#' fit <- DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8)
#' diag <- DPprior_diagnostics(fit)
#' print(diag)
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @family diagnostics
#'
#' @export
DPprior_diagnostics <- function(fit, thresholds = c(0.5, 0.9)) {

  # Input validation
  if (is.null(fit) || !is.list(fit)) {
    stop("fit must be a list-like DPprior_fit object.", call. = FALSE)
  }
  if (is.null(fit$a) || is.null(fit$b) || is.null(fit$J)) {
    stop("fit must contain fields a, b, and J.", call. = FALSE)
  }

  # Extract and coerce parameters
  a <- as.numeric(fit$a)
  b <- as.numeric(fit$b)
  J <- as.integer(fit$J)

  # Allow optional M from fit object
  M <- if (!is.null(fit$M)) as.integer(fit$M) else .QUAD_NODES_DEFAULT

  assert_positive(a, "a")
  assert_positive(b, "b")
  assert_valid_J(J)

  # Compute all diagnostics
  alpha_diag <- compute_alpha_diagnostics(a, b)
  K_diag <- compute_K_diagnostics(J, a, b, M)
  weight_diag <- compute_weight_diagnostics(a, b, thresholds, M)
  coclust_diag <- compute_coclustering_diagnostics(a, b, M)

  # Compile warnings
  warnings <- character(0)

  # Check for high dominance risk (always check 0.5 threshold)
  p_gt_05 <- prob_w1_exceeds(0.5, a, b)
  p_gt_09 <- prob_w1_exceeds(0.9, a, b)

  if (p_gt_05 > 0.4) {
    warnings <- c(warnings,
                  sprintf("HIGH DOMINANCE RISK: P(w1 > 0.5) = %.1f%% exceeds 40%%",
                          100 * p_gt_05))
  }

  if (p_gt_09 > 0.15) {
    warnings <- c(warnings,
                  sprintf("NEAR-DEGENERATE RISK: P(w1 > 0.9) = %.1f%% exceeds 15%%",
                          100 * p_gt_09))
  }

  # Consistency check: implied alpha from K vs weights
  if (is.finite(weight_diag$mean) && weight_diag$mean > 0 &&
      weight_diag$mean < 1 && J > 1L) {
    alpha_from_weights <- (1 - weight_diag$mean) / weight_diag$mean
    if (is.finite(K_diag$mean) && log(J) > 0) {
      alpha_from_K <- (K_diag$mean - 1) / log(J)
      if (is.finite(alpha_from_K) && alpha_from_K > 0) {
        if (abs(log(alpha_from_weights / alpha_from_K)) > 0.5) {
          warnings <- c(warnings,
                        sprintf("NOTE: Implied alpha from K_J (%.2f) vs weights (%.2f) differs substantially",
                                alpha_from_K, alpha_from_weights))
        }
      }
    }
  }

  result <- list(
    J = J,
    a = a,
    b = b,
    alpha = alpha_diag,
    K = K_diag,
    weights = weight_diag,
    coclustering = coclust_diag,
    warnings = warnings
  )

  class(result) <- "DPprior_diagnostics"
  result
}


# =============================================================================
# S3 Methods for DPprior_diagnostics
# =============================================================================

#' Print Method for DPprior_diagnostics Objects
#'
#' @param x An object of class "DPprior_diagnostics".
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.DPprior_diagnostics <- function(x, ...) {
  cat("DPprior Comprehensive Diagnostics\n")
  cat(strrep("=", 60), "\n\n")

  cat(sprintf("Prior: alpha ~ Gamma(%.4f, %.4f) for J = %d\n\n", x$a, x$b, x$J))

  # Alpha summary
  cat("alpha Distribution:\n")
  cat(strrep("-", 40), "\n")
  cat(sprintf("  E[alpha] = %.3f, CV(alpha) = %.3f, Median = %.3f\n",
              x$alpha$mean, x$alpha$cv, x$alpha$median))
  cat(sprintf("  90%% CI: [%.3f, %.3f]\n\n",
              x$alpha$quantiles["q5"], x$alpha$quantiles["q95"]))

  # K summary
  cat("K_J Distribution:\n")
  cat(strrep("-", 40), "\n")
  cat(sprintf("  E[K] = %.2f, SD(K) = %.2f, Mode = %d\n",
              x$K$mean, x$K$sd, x$K$mode))
  cat(sprintf("  Median = %d, IQR = [%d, %d]\n\n",
              x$K$median, x$K$quantiles["q25"], x$K$quantiles["q75"]))

  # Weight summary (KEY DIAGNOSTIC)
  cat("w1 Distribution (Size-Biased First Weight):\n")
  cat(strrep("-", 40), "\n")

  p_gt_05 <- prob_w1_exceeds(0.5, x$a, x$b)
  p_gt_09 <- prob_w1_exceeds(0.9, x$a, x$b)

  cat(sprintf("  E[w1] = %.3f, Median = %.3f\n",
              x$weights$mean, x$weights$median))
  cat(sprintf("  P(w1 > 0.5) = %.1f%% (dominance risk: %s)\n",
              100 * p_gt_05, toupper(x$weights$dominance_risk)))
  cat(sprintf("  P(w1 > 0.9) = %.1f%%\n\n", 100 * p_gt_09))

  # Co-clustering
  cat("Co-Clustering (rho = sum w_h^2):\n")
  cat(strrep("-", 40), "\n")
  cat(sprintf("  E[rho] = %.3f (%s)\n\n",
              x$coclustering$mean, x$coclustering$interpretation))

  # Warnings
  if (length(x$warnings) > 0) {
    cat("WARNINGS:\n")
    cat(strrep("-", 40), "\n")
    for (w in x$warnings) {
      cat(sprintf("  * %s\n", w))
    }
    cat("\n  Consider using DPprior_dual() for weight-constrained elicitation.\n")
  } else {
    cat("No diagnostic warnings.\n")
  }

  invisible(x)
}


#' Summary Method for DPprior_diagnostics Objects
#'
#' @param object An object of class "DPprior_diagnostics".
#' @param ... Additional arguments (ignored).
#'
#' @return A data frame with key diagnostic metrics.
#'
#' @export
summary.DPprior_diagnostics <- function(object, ...) {
  p_gt_05 <- prob_w1_exceeds(0.5, object$a, object$b)
  p_gt_09 <- prob_w1_exceeds(0.9, object$a, object$b)

  data.frame(
    J = object$J,
    a = object$a,
    b = object$b,
    E_alpha = object$alpha$mean,
    CV_alpha = object$alpha$cv,
    E_K = object$K$mean,
    SD_K = object$K$sd,
    Mode_K = object$K$mode,
    E_w1 = object$weights$mean,
    P_w1_gt_50 = p_gt_05,
    P_w1_gt_90 = p_gt_09,
    dominance_risk = object$weights$dominance_risk,
    E_rho = object$coclustering$mean,
    n_warnings = length(object$warnings),
    row.names = NULL
  )
}


# =============================================================================
# Comparison Diagnostics
# =============================================================================

#' Compare Diagnostics Across Multiple Fits
#'
#' Creates a comparison table of diagnostic metrics for multiple DPprior fits.
#'
#' @param ... Named DPprior_fit objects to compare.
#' @param M Integer; number of quadrature nodes (default: 80).
#'
#' @return A data frame with diagnostic metrics for each fit.
#'
#' @examples
#' \dontrun{
#' fit_K <- DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8)
#' fit_dual <- DPprior_dual(fit_K,
#'                          w1_target = list(prob = list(threshold = 0.5, value = 0.3)))
#' compare_diagnostics(K_only = fit_K, Dual_anchor = fit_dual)
#'
#' }
#' @keywords internal
compare_diagnostics <- function(..., M = .QUAD_NODES_DEFAULT) {
  fits <- list(...)

  if (length(fits) == 0L) {
    stop("At least one fit object required", call. = FALSE)
  }

  fit_names <- names(fits)
  if (is.null(fit_names) || any(fit_names == "")) {
    fit_names <- paste0("Fit_", seq_along(fits))
  }

  results <- lapply(seq_along(fits), function(i) {
    fit <- fits[[i]]
    diag <- DPprior_diagnostics(fit, thresholds = c(0.5, 0.9))

    data.frame(
      method = fit_names[i],
      a = diag$a,
      b = diag$b,
      E_K = diag$K$mean,
      SD_K = diag$K$sd,
      E_w1 = diag$weights$mean,
      P_w1_gt_50 = prob_w1_exceeds(0.5, diag$a, diag$b),
      P_w1_gt_90 = prob_w1_exceeds(0.9, diag$a, diag$b),
      risk = diag$weights$dominance_risk,
      E_rho = diag$coclustering$mean,
      warnings = length(diag$warnings),
      row.names = NULL
    )
  })

  do.call(rbind, results)
}


# =============================================================================
# Verification Function
# =============================================================================

#' Verify Diagnostics Module
#'
#' Runs comprehensive tests on the diagnostics module.
#'
#' @param verbose Logical; if TRUE, print detailed output.
#'
#' @return Invisibly returns TRUE if all tests pass.
#'
#' @keywords internal
verify_diagnostics <- function(verbose = TRUE) {

  if (isTRUE(verbose)) {
    cat("=== Diagnostics Module Verification ===\n\n")
  }

  # Test 1: E[w1] = E[rho] identity
  if (isTRUE(verbose)) cat("Test 1: Verify E[w1] = E[rho] identity\n")

  test_cases <- list(
    c(a = 0.5, b = 0.5),
    c(a = 1, b = 1),
    c(a = 2, b = 1),
    c(a = 1.6, b = 1.22),
    c(a = 5, b = 2)
  )

  for (params in test_cases) {
    a <- params["a"]
    b <- params["b"]
    E_w1 <- mean_w1(a, b)
    E_rho <- mean_rho(a, b)

    stopifnot(abs(E_w1 - E_rho) < 1e-6)

    if (isTRUE(verbose)) {
      cat(sprintf("  a=%.1f, b=%.2f: E[w1]=%.6f, E[rho]=%.6f - PASS\n",
                  a, b, E_w1, E_rho))
    }
  }
  if (isTRUE(verbose)) cat("\n")

  # Test 2: High dominance detection
  if (isTRUE(verbose)) cat("Test 2: Detect high dominance (Lee et al. prior)\n")

  diag <- compute_weight_diagnostics(1.60, 1.22)
  stopifnot(diag$dominance_risk == "high")
  stopifnot(diag$prob_exceeds["prob_gt_0.5"] > 0.4)

  if (isTRUE(verbose)) {
    cat(sprintf("  P(w1 > 0.5) = %.3f, risk = %s - PASS\n\n",
                diag$prob_exceeds["prob_gt_0.5"], diag$dominance_risk))
  }

  # Test 3: Quantile-CDF consistency
  if (isTRUE(verbose)) cat("Test 3: Quantile-CDF consistency\n")

  a <- 2; b <- 1
  for (p in c(0.1, 0.25, 0.5, 0.75, 0.9)) {
    q <- quantile_w1(p, a, b)
    cdf_at_q <- cdf_w1(q, a, b)
    stopifnot(abs(cdf_at_q - p) < 1e-10)
    if (isTRUE(verbose)) {
      cat(sprintf("  p=%.2f: Q(p)=%.6f, F(Q(p))=%.10f - PASS\n", p, q, cdf_at_q))
    }
  }
  if (isTRUE(verbose)) cat("\n")

  # Test 4: Reference values (from Python/mpmath)
  if (isTRUE(verbose)) cat("Test 4: Reference values verification\n")

  # P(w1 > 0.5) for a=1.6, b=1.22 should be ~0.4868
  ref_p_05 <- prob_w1_exceeds(0.5, 1.6, 1.22)
  stopifnot(abs(ref_p_05 - 0.4868311) < 1e-5)

  if (isTRUE(verbose)) {
    cat(sprintf("  P(w1>0.5) @ (1.6,1.22): expected=0.4868, got=%.6f - PASS\n\n", ref_p_05))
  }

  if (isTRUE(verbose)) {
    cat("=== All Verification Tests Passed ===\n")
  }

  invisible(TRUE)
}

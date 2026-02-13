# =============================================================================
# Module 17: S3 Methods for DPprior_fit Class
# =============================================================================
#
# This module implements enhanced S3 methods (print, summary, plot) for the
# DPprior_fit class, providing user-friendly output and visualization.
#
# Key Features:
# - print.DPprior_fit(): Concise one-line-per-concept output with dominance risk
# - summary.DPprior_fit(): Detailed comparison of target vs achieved fit
# - plot.DPprior_fit(): Flexible visualization with multiple plot types
#
# Author: JoonHo Lee
# Date: December 2025
# Part of: DPprior R Package
# Dependencies: Modules 00 (constants), 08 (weights_w1), 14 (diagnostics),
#               15 (visualization)
# =============================================================================


# =============================================================================
# S3 Method: print.DPprior_fit()
# =============================================================================

#' Print Method for DPprior_fit Objects
#'
#' Displays a concise, informative summary of a prior elicitation result,
#' including the Gamma hyperprior specification, target vs achieved fit,
#' and dominance risk assessment.
#'
#' @param x A \code{DPprior_fit} object.
#' @param digits Integer; number of significant digits for display.
#'   Default is 4.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @details
#' The output includes:
#' \itemize{
#'   \item Gamma hyperprior parameters (a, b) with moments E[alpha] and SD[alpha]
#'   \item Target specification (J, \eqn{E[K_J]}, \eqn{Var(K_J)})
#'   \item Achieved fit with residual error
#'   \item Method used and iteration count
#'   \item Quick dominance risk summary (if diagnostics available)
#' }
#'
#' @section Dominance Risk:
#' If diagnostics are computed, the dominance risk is displayed as:
#' \itemize{
#'   \item \strong{LOW}: P(w1 > 0.5) < 20\%
#'   \item \strong{MODERATE}: 20\% <= P(w1 > 0.5) < 40\%
#'   \item \strong{HIGH}: P(w1 > 0.5) >= 40\%
#' }
#'
#' @examples
#' # Create a fit object
#' fit <- DPprior_fit(J = 50, mu_K = 5, var_K = 8)
#' print(fit)
#'
#' # With custom digits
#' print(fit, digits = 6)
#'
#' @seealso \code{\link{summary.DPprior_fit}}, \code{\link{plot.DPprior_fit}},
#'   \code{\link{DPprior_fit}}
#'
#' @method print DPprior_fit
#' @export
print.DPprior_fit <- function(x, digits = 4, ...) {

  # Header
  cat("DPprior Prior Elicitation Result\n")
  cat(strrep("=", 45), "\n\n")

  # Extract core parameters
  a <- x$a
  b <- x$b
  J <- x$J

  # Gamma hyperprior specification
  E_alpha <- a / b
  SD_alpha <- sqrt(a) / b

  cat(sprintf("Gamma Hyperprior: \u03b1 ~ Gamma(a = %.*f, b = %.*f)\n",
              digits, a, digits, b))
  cat(sprintf("  E[\u03b1] = %.3f, SD[\u03b1] = %.3f\n\n", E_alpha, SD_alpha))

  # Target specification
  cat(sprintf("Target (J = %d):\n", J))
  if (!is.null(x$target)) {
    target_mu <- x$target$mu_K
    target_var <- x$target$var_K
    target_var_used <- x$target$var_K_used

    cat(sprintf("  E[K_J]   = %.2f\n", target_mu))

    # Show projection if var_K differs from var_K_used
    if (!is.null(target_var_used) && !is.null(target_var) &&
        abs(target_var - target_var_used) > 1e-10) {
      cat(sprintf("  Var(K_J) = %.2f (requested), %.4f (used after projection)\n",
                  target_var, target_var_used))
    } else {
      cat(sprintf("  Var(K_J) = %.2f\n", target_var))
    }

    # Show confidence level if specified
    if (!is.null(x$target$confidence) && nzchar(as.character(x$target$confidence))) {
      cat(sprintf("  (from confidence = '%s')\n", x$target$confidence))
    }
  }

  # Achieved fit
  if (!is.null(x$fit)) {
    cat("\nAchieved:\n")
    cat(sprintf("  E[K_J] = %.*f, Var(K_J) = %.*f\n",
                digits + 2, x$fit$mu_K, digits + 2, x$fit$var_K))

    if (!is.null(x$fit$residual) && !is.na(x$fit$residual)) {
      cat(sprintf("  Residual = %.2e\n", x$fit$residual))
    }
  }

  # Method and iterations
  cat(sprintf("\nMethod: %s", x$method))
  if (!is.null(x$iterations) && !is.na(x$iterations)) {
    cat(sprintf(" (%d iterations)", x$iterations))
  }
  cat("\n")

  # Dominance risk summary (quick view)
  if (!is.null(x$diagnostics) && !is.null(x$diagnostics$weights)) {
    weights_diag <- x$diagnostics$weights
    dr <- weights_diag$dominance_risk
    prob_dom <- weights_diag$prob_exceeds["prob_gt_0.5"]

    # Risk symbol for visual cue
    risk_symbol <- switch(
      dr,
      "low" = "\u2713",        # checkmark
      "moderate" = "\u26A0",   # warning
      "high" = "\u2718",       # X mark
      "?"
    )

    cat(sprintf("\nDominance Risk: %s %s (P(w\u2081>0.5) = %.0f%%)\n",
                toupper(dr), risk_symbol, 100 * prob_dom))
  }

  invisible(x)
}


# =============================================================================
# S3 Method: summary.DPprior_fit()
# =============================================================================

#' Summary Method for DPprior_fit Objects
#'
#' Produces a comprehensive summary of a prior elicitation result, including
#' detailed parameter information, target vs achieved comparison, and full
#' diagnostic statistics.
#'
#' @param object A \code{DPprior_fit} object.
#' @param print_output Logical; if \code{TRUE} (default), prints the summary
#'   to the console. If \code{FALSE}, returns the summary list silently.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class \code{"summary.DPprior_fit"} containing:
#'   \itemize{
#'     \item \code{method}: Calibration method used
#'     \item \code{status}: Convergence status
#'     \item \code{gamma_prior}: List with (a, b) parameters
#'     \item \code{alpha_summary}: List with E[alpha], Var[alpha], CV[alpha]
#'     \item \code{target}: Target specification
#'     \item \code{achieved}: Achieved fit
#'     \item \code{errors}: Absolute and relative fitting errors
#'     \item \code{converged}: Logical convergence flag
#'     \item \code{iterations}: Number of iterations
#'     \item \code{diagnostics}: Full diagnostic information (if available)
#'   }
#'
#' @details
#' The printed summary includes:
#' \enumerate{
#'   \item Basic information (J, method)
#'   \item Gamma hyperprior parameters with derived statistics
#'   \item Target vs Achieved comparison table with error metrics
#'   \item Full diagnostics for alpha, K, and w1 distributions (if computed)
#'   \item Iteration trace (first 10 rows, if available)
#' }
#'
#' @examples
#' # Create a fit object
#' fit <- DPprior_fit(J = 50, mu_K = 5, var_K = 8, check_diagnostics = TRUE)
#' summary(fit)
#'
#' # Store summary without printing
#' summ <- summary(fit, print_output = FALSE)
#' str(summ)
#'
#' @seealso \code{\link{print.DPprior_fit}}, \code{\link{DPprior_diagnostics}}
#'
#' @method summary DPprior_fit
#' @export
summary.DPprior_fit <- function(object, print_output = TRUE, ...) {

  # Helper function to safely extract scalar values
  safe_scalar <- function(x, default = NA) {
    if (is.null(x) || length(x) == 0L) return(default)
    if (length(x) == 1L) return(x)
    return(x[1L])
  }

  # Extract core parameters
  a <- safe_scalar(object$a, NA_real_)
  b <- safe_scalar(object$b, NA_real_)
  J <- safe_scalar(object$J, NA_integer_)

  # Compute derived quantities
  E_alpha <- if (!is.na(a) && !is.na(b) && b > 0) a / b else NA_real_
  Var_alpha <- if (!is.na(a) && !is.na(b) && b > 0) a / b^2 else NA_real_
  CV_alpha <- if (!is.na(a) && a > 0) 1 / sqrt(a) else NA_real_

  # Target
  target_mu_K <- safe_scalar(object$target$mu_K, NA_real_)
  target_var_K <- safe_scalar(object$target$var_K, NA_real_)
  target_var_K_used <- safe_scalar(
    object$target$var_K_used %||% object$target$var_K,
    NA_real_
  )

  # Achieved fit
  achieved_mu_K <- safe_scalar(object$fit$mu_K, NA_real_)
  achieved_var_K <- safe_scalar(object$fit$var_K, NA_real_)

  # Compute errors
  error_mu_abs <- abs(target_mu_K - achieved_mu_K)
  error_var_abs <- abs(target_var_K_used - achieved_var_K)
  error_mu_rel <- if (!is.na(target_mu_K) && target_mu_K > 0) {
    100 * error_mu_abs / target_mu_K
  } else NA_real_
  error_var_rel <- if (!is.na(target_var_K_used) && target_var_K_used > 0) {
    100 * error_var_abs / target_var_K_used
  } else NA_real_

  # Build result structure
  result <- list(
    method = safe_scalar(object$method, NA_character_),
    status = safe_scalar(object$status, "success"),
    gamma_prior = list(
      a = a,
      b = b
    ),
    alpha_summary = list(
      E_alpha = E_alpha,
      Var_alpha = Var_alpha,
      SD_alpha = sqrt(Var_alpha),
      CV_alpha = CV_alpha
    ),
    target = list(
      mu_K = target_mu_K,
      var_K = target_var_K,
      var_K_used = target_var_K_used,
      confidence = safe_scalar(object$target$confidence, NA_character_)
    ),
    achieved = list(
      mu_K = achieved_mu_K,
      var_K = achieved_var_K,
      residual = safe_scalar(object$fit$residual, NA_real_)
    ),
    errors = list(
      mu_K_abs = error_mu_abs,
      var_K_abs = error_var_abs,
      mu_K_rel_pct = error_mu_rel,
      var_K_rel_pct = error_var_rel
    ),
    scaling = list(
      J = J
    ),
    converged = safe_scalar(object$converged, NA),
    iterations = safe_scalar(object$iterations, NA_integer_),
    diagnostics = object$diagnostics,
    trace = object$trace,
    dual_anchor = object$dual_anchor
  )

  class(result) <- "summary.DPprior_fit"

  # Print output if requested
  if (isTRUE(print_output)) {
    print(result)
  }

  invisible(result)
}


#' Print Method for summary.DPprior_fit
#'
#' @param x A \code{summary.DPprior_fit} object.
#' @param diagnostics Logical; if TRUE, print full diagnostics. Default is FALSE.
#' @param max_trace Integer; maximum number of trace rows to display. Default is 10.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @method print summary.DPprior_fit
#' @export
print.summary.DPprior_fit <- function(x, diagnostics = FALSE, max_trace = 10L, ...) {

  # Inline helper for number formatting
  fmt_num <- function(v, digits = 6) {
    if (!is.finite(v)) return("NA")
    formatC(v, format = "f", digits = digits)
  }

  # Header
  cat("DPprior Prior Elicitation Summary\n")
  cat(strrep("=", 60), "\n\n")

  # Basic info
  cat(sprintf("Sample size: J = %d\n", x$scaling$J))
  cat(sprintf("Method: %s\n", x$method))
  if (!is.na(x$status)) {
    cat(sprintf("Status: %s\n", x$status))
  }
  cat("\n")

  # Gamma hyperprior
  cat("Gamma Hyperprior:\n")
  cat(strrep("-", 40), "\n")
  cat(sprintf("  Shape (a) = %.6f\n", x$gamma_prior$a))
  cat(sprintf("  Rate (b)  = %.6f\n", x$gamma_prior$b))
  cat(sprintf("  E[\u03b1] = %.4f, SD[\u03b1] = %.4f, CV[\u03b1] = %.4f\n",
              x$alpha_summary$E_alpha,
              x$alpha_summary$SD_alpha,
              x$alpha_summary$CV_alpha))
  cat("\n")

  # Target vs Achieved comparison table
  cat("Target vs Achieved:\n")
  cat(strrep("-", 40), "\n")
  cat(sprintf("  %-15s %12s %12s %12s\n", "", "Target", "Achieved", "Error"))
  cat(sprintf("  %-15s %12.4f %12.4f %12.2e\n",
              "E[K_J]", x$target$mu_K, x$achieved$mu_K, x$errors$mu_K_abs))
  cat(sprintf("  %-15s %12.4f %12.4f %12.2e\n",
              "Var(K_J)", x$target$var_K_used, x$achieved$var_K, x$errors$var_K_abs))

  # Note if var_K was projected
  if (!is.na(x$target$var_K) && !is.na(x$target$var_K_used)) {
    if (abs(x$target$var_K - x$target$var_K_used) > 1e-6) {
      cat(sprintf("\n  Note: Original var_K=%.2f projected to %.2f (feasibility)\n",
                  x$target$var_K, x$target$var_K_used))
    }
  }
  cat("\n")

  # Full diagnostics if available
  if (!is.null(x$diagnostics)) {
    cat("Diagnostics:\n")
    cat(strrep("-", 40), "\n")

    # Alpha diagnostics
    if (!is.null(x$diagnostics$alpha)) {
      alpha_diag <- x$diagnostics$alpha
      q <- alpha_diag$quantiles
      cat(sprintf("  Alpha:\n"))
      cat(sprintf("    Mean = %.4f, SD = %.4f, CV = %.4f\n",
                  alpha_diag$mean, alpha_diag$sd, alpha_diag$cv))
      cat(sprintf("    90%% CI: [%.3f, %.3f]\n", q["q5"], q["q95"]))
    }

    # K diagnostics
    if (!is.null(x$diagnostics$K)) {
      K_diag <- x$diagnostics$K
      cat(sprintf("  K_J:\n"))
      cat(sprintf("    Mean = %.2f, SD = %.2f, Mode = %d\n",
                  K_diag$mean, K_diag$sd, K_diag$mode))
    }

    # Weight diagnostics (critical for dominance)
    if (!is.null(x$diagnostics$weights)) {
      w_diag <- x$diagnostics$weights
      cat(sprintf("  First weight (w\u2081):\n"))
      cat(sprintf("    Mean = %.4f\n", w_diag$mean))
      cat(sprintf("    P(w\u2081 > 0.5) = %.1f%%\n",
                  100 * w_diag$prob_exceeds["prob_gt_0.5"]))
      cat(sprintf("    P(w\u2081 > 0.9) = %.1f%%\n",
                  100 * w_diag$prob_exceeds["prob_gt_0.9"]))
      cat(sprintf("    Dominance risk: %s\n", toupper(w_diag$dominance_risk)))
    }
    cat("\n")
  }

  # Dual-anchor info if available
  if (!is.null(x$dual_anchor)) {
    cat("Dual-Anchor Calibration:\n")
    cat(strrep("-", 40), "\n")
    cat(sprintf("  Lambda: %.2f\n", x$dual_anchor$lambda))
    cat(sprintf("  Loss type: %s\n", x$dual_anchor$loss_type))
    cat(sprintf("  K loss: %.4e\n", x$dual_anchor$K_loss))
    cat(sprintf("  w1 loss: %.4e\n", x$dual_anchor$w_loss))

    if (!is.null(x$dual_anchor$w1_achieved)) {
      wa <- x$dual_anchor$w1_achieved
      cat(sprintf("  Achieved P(w\u2081 > 0.5): %.1f%%\n", 100 * wa$prob_gt_50))
    }

    if (!is.null(x$dual_anchor$init)) {
      cat(sprintf("  Initial (K-only): Gamma(%.4f, %.4f)\n",
                  x$dual_anchor$init$a, x$dual_anchor$init$b))
    }
    cat("\n")
  }

  # Iteration trace if available
  if (!is.null(x$trace) && is.data.frame(x$trace) && nrow(x$trace) > 0) {
    cat("Iteration Trace:\n")
    cat(strrep("-", 40), "\n")
    print(utils::head(x$trace, max_trace))
    if (nrow(x$trace) > max_trace) {
      cat(sprintf("... (%d more rows)\n", nrow(x$trace) - max_trace))
    }
    cat("\n")
  }

  # Full diagnostics (if requested)
  if (isTRUE(diagnostics) && !is.null(x$diagnostics)) {
    cat("Full Diagnostics:\n")
    cat(strrep("-", 40), "\n")
    print(x$diagnostics)
    cat("\n")
  } else if (!is.null(x$diagnostics) && !is.null(x$diagnostics$weights)) {
    # Quick reminder about available diagnostics
    w <- x$diagnostics$weights
    if (!is.null(w$dominance_risk)) {
      cat(sprintf("Diagnostics: dominance risk = %s (use diagnostics=TRUE for full report)\n\n",
                  toupper(as.character(w$dominance_risk))))
    }
  }

  invisible(x)
}


# =============================================================================
# S3 Method: plot.DPprior_fit()
# =============================================================================

#' Plot Method for DPprior_fit Objects
#'
#' Creates visualizations of a prior elicitation result. Multiple plot types
#' are available, including individual distribution plots and comprehensive
#' dashboards.
#'
#' @param x A \code{DPprior_fit} object.
#' @param type Character; the type of plot to create:
#'   \describe{
#'     \item{"auto"}{(Default) Automatically selects the appropriate plot type.
#'       Uses \code{"dual"} for dual-anchor fits, \code{"dashboard"} otherwise.}
#'     \item{"dashboard"}{4-panel dashboard showing alpha, K, w1, and summary.}
#'     \item{"alpha"}{Prior density of the concentration parameter alpha.}
#'     \item{"K"}{Prior PMF of the number of clusters \eqn{K_J}.}
#'     \item{"w1"}{Prior density of the first stick-breaking weight w1.}
#'     \item{"dual"}{Dual-anchor comparison dashboard (for dual-anchor fits).}
#'     \item{"comparison"}{Same as "dual".}
#'   }
#' @param engine Character; graphics engine to use:
#'   \code{"ggplot2"} (default) or \code{"base"}.
#' @param ... Additional arguments passed to the underlying plot functions.
#'   Common options include:
#'   \describe{
#'     \item{base_size}{Base font size (default: 11)}
#'     \item{ci_level}{Credible interval level for alpha plot (default: 0.95)}
#'     \item{title}{Optional title for the dashboard}
#'     \item{show}{If TRUE, display the plot; if FALSE, return silently}
#'   }
#'
#' @return Depends on the plot type and engine:
#'   \itemize{
#'     \item For ggplot2: Returns a ggplot object or gtable (for dashboards)
#'     \item For base: Returns invisible(NULL)
#'   }
#'
#' @details
#' The \code{"auto"} type is recommended for most use cases. It automatically
#' detects whether the fit object is from dual-anchor calibration and selects
#' the appropriate visualization.
#'
#' For dual-anchor fits, the comparison dashboard shows:
#' \itemize{
#'   \item Alpha prior: K-only vs Dual-anchor
#'   \item K distribution comparison
#'   \item w1 distribution comparison with dominance threshold
#'   \item Summary comparison table
#' }
#'
#' @section Plot Type Details:
#' \describe{
#'   \item{dashboard}{
#'     A 2x2 grid showing:
#'     (A) Alpha prior density with CI
#'     (B) \eqn{K_J} prior PMF with mode and mean
#'     (C) w1 prior density with dominance shading
#'     (D) Summary statistics table
#'   }
#'   \item{alpha}{
#'     Gamma(a, b) density with:
#'     - Mean line (dashed)
#'     - Credible interval (shaded region)
#'     - Annotation with moments and CI
#'   }
#'   \item{K}{
#'     Bar plot of \eqn{P(K_J = k)} with:
#'     - Target mean line
#'     - Achieved mean line
#'     - Optional CDF overlay
#'   }
#'   \item{w1}{
#'     Density plot with:
#'     - Dominance region shading (w1 > 0.5)
#'     - Threshold lines
#'     - Exceedance probabilities
#'   }
#' }
#'
#' @examples
#' # Create a fit object
#' fit <- DPprior_fit(J = 50, mu_K = 5, var_K = 8)
#'
#' # Auto-detect best plot type
#' plot(fit)
#'
#' # Specific plot types
#' plot(fit, type = "alpha")
#' plot(fit, type = "K")
#' plot(fit, type = "w1")
#' plot(fit, type = "dashboard")
#'
#' # With custom options
#' plot(fit, type = "dashboard", title = "My Prior Analysis")
#'
#' # Dual-anchor comparison
#' fit_K <- DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8)
#' fit_dual <- DPprior_dual(fit_K, w1_target = list(prob = list(threshold = 0.5, value = 0.3)))
#' plot(fit_dual)  # Auto-selects dual comparison
#' plot(fit_dual, type = "comparison")  # Explicit
#'
#' @seealso \code{\link{plot_prior_dashboard}}, \code{\link{plot_alpha_prior}},
#'   \code{\link{plot_K_prior}}, \code{\link{plot_w1_prior}},
#'   \code{\link{plot_dual_comparison}}
#'
#' @method plot DPprior_fit
#' @export
plot.DPprior_fit <- function(x, type = c("auto", "dashboard", "alpha", "K", "w1",
                                         "dual", "comparison"),
                             engine = c("ggplot2", "base"),
                             ...) {
  type <- match.arg(type)
  engine <- match.arg(engine)

  # Auto-detect: use dual comparison for dual-anchor fits
  if (type == "auto") {
    type <- if (.dpprior_is_dual(x)) "dual" else "dashboard"
  }

  # Dispatch to appropriate plot function with existence checks
  switch(
    type,
    "dashboard" = {
      if (exists("plot_prior_dashboard", mode = "function")) {
        plot_prior_dashboard(x, engine = engine, ...)
      } else {
        stop("plot_prior_dashboard() not available. Required visualization function not found. Please reinstall DPprior.", call. = FALSE)
      }
    },
    "alpha" = {
      if (exists("plot_alpha_prior", mode = "function")) {
        plot_alpha_prior(x, engine = engine, ...)
      } else {
        stop("plot_alpha_prior() not available. Required visualization function not found. Please reinstall DPprior.", call. = FALSE)
      }
    },
    "K" = {
      if (exists("plot_K_prior", mode = "function")) {
        plot_K_prior(x, engine = engine, ...)
      } else {
        stop("plot_K_prior() not available. Required visualization function not found. Please reinstall DPprior.", call. = FALSE)
      }
    },
    "w1" = {
      if (exists("plot_w1_prior", mode = "function")) {
        plot_w1_prior(x, engine = engine, ...)
      } else {
        stop("plot_w1_prior() not available. Required visualization function not found. Please reinstall DPprior.", call. = FALSE)
      }
    },
    "dual" = {
      if (.dpprior_is_dual(x) && exists("plot_dual_comparison", mode = "function")) {
        plot_dual_comparison(x, engine = engine, ...)
      } else if (exists("plot_prior_dashboard", mode = "function")) {
        if (!.dpprior_is_dual(x)) {
          warning("Not a dual-anchor fit; showing standard dashboard.", call. = FALSE)
        }
        plot_prior_dashboard(x, engine = engine, ...)
      } else {
        stop("No plotting functions available. Required visualization function not found. Please reinstall DPprior.", call. = FALSE)
      }
    },
    "comparison" = {
      if (.dpprior_is_dual(x) && exists("plot_dual_comparison", mode = "function")) {
        plot_dual_comparison(x, engine = engine, ...)
      } else {
        warning("Not a dual-anchor fit or plot_dual_comparison() not available; showing standard dashboard.",
                call. = FALSE)
        if (exists("plot_prior_dashboard", mode = "function")) {
          plot_prior_dashboard(x, engine = engine, ...)
        } else {
          stop("No plotting functions available. Required visualization function not found. Please reinstall DPprior.", call. = FALSE)
        }
      }
    }
  )
}


# =============================================================================
# Helper Functions for S3 Methods
# =============================================================================

#' Check if Fit is from Dual-Anchor Calibration
#'
#' Checks whether a DPprior_fit object originated from DPprior_dual().
#' A stricter check ensures the dual_anchor component has the expected structure.
#'
#' @param fit A DPprior_fit object.
#' @return Logical; TRUE if dual-anchor fit with valid structure.
#'
#' @keywords internal
.dpprior_is_dual <- function(fit) {
  # Stricter check: requires both dual_anchor and dual_anchor$init
  # This ensures the fit actually came from DPprior_dual() with proper initialization
  !is.null(fit$dual_anchor) && !is.null(fit$dual_anchor$init)
}


# =============================================================================
# Verification Function
# =============================================================================

#' Verify S3 Methods Module
#'
#' Runs comprehensive verification tests for the S3 methods module.
#'
#' @param verbose Logical; if TRUE, print detailed test output.
#'
#' @return Invisibly returns TRUE if all tests pass.
#'
#' @examples
#' \dontrun{
#' verify_s3_methods()
#'
#' }
#' @keywords internal
verify_s3_methods <- function(verbose = TRUE) {

  if (isTRUE(verbose)) {
    cat(strrep("=", 60), "\n")
    cat("Module 17: S3 Methods - Verification Suite\n")
    cat(strrep("=", 60), "\n\n")
  }

  all_pass <- TRUE

  # -------------------------------------------------------------------------
  # Test 1: Create mock DPprior_fit object
  # -------------------------------------------------------------------------
  if (isTRUE(verbose)) {
    cat("[Test 1] Create mock DPprior_fit object\n")
    cat(strrep("-", 50), "\n")
  }

  test_fit <- list(
    J = 50L,
    a = 1.597,
    b = 1.222,
    target = list(mu_K = 5.0, var_K = 8.0, var_K_used = 8.0),
    fit = list(mu_K = 5.000001, var_K = 7.999998, residual = 1.5e-10),
    method = "A2-MN",
    converged = TRUE,
    iterations = 5L,
    status = "success"
  )
  class(test_fit) <- "DPprior_fit"

  # Add diagnostics
  test_fit$diagnostics <- list(
    alpha = list(
      mean = 1.307,
      sd = 1.034,
      cv = 0.791,
      quantiles = c(q5 = 0.15, q25 = 0.55, q50 = 1.05, q75 = 1.75, q95 = 3.25)
    ),
    weights = list(
      mean = 0.509,
      prob_exceeds = c(
        "prob_gt_0.3" = 0.72,
        "prob_gt_0.5" = 0.49,
        "prob_gt_0.7" = 0.25,
        "prob_gt_0.9" = 0.08
      ),
      dominance_risk = "high"
    )
  )

  test1_pass <- inherits(test_fit, "DPprior_fit") &&
    !is.null(test_fit$a) &&
    !is.null(test_fit$b)

  if (isTRUE(verbose)) {
    cat(sprintf("  Mock object created: %s\n", if (test1_pass) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && test1_pass

  # -------------------------------------------------------------------------
  # Test 2: print.DPprior_fit()
  # -------------------------------------------------------------------------
  if (isTRUE(verbose)) {
    cat("\n[Test 2] print.DPprior_fit()\n")
    cat(strrep("-", 50), "\n")
  }

  test2_pass <- tryCatch({
    output <- capture.output(print(test_fit))
    has_header <- any(grepl("DPprior Prior Elicitation Result", output))
    has_gamma <- any(grepl("Gamma Hyperprior", output))
    has_target <- any(grepl("Target", output))
    has_achieved <- any(grepl("Achieved", output))
    has_dominance <- any(grepl("Dominance Risk", output))

    all(has_header, has_gamma, has_target, has_achieved, has_dominance)
  }, error = function(e) FALSE)

  if (isTRUE(verbose)) {
    cat(sprintf("  print() works: %s\n", if (test2_pass) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && test2_pass

  # -------------------------------------------------------------------------
  # Test 3: summary.DPprior_fit()
  # -------------------------------------------------------------------------
  if (isTRUE(verbose)) {
    cat("\n[Test 3] summary.DPprior_fit()\n")
    cat(strrep("-", 50), "\n")
  }

  test3_pass <- tryCatch({
    summ <- summary(test_fit, print_output = FALSE)

    has_class <- inherits(summ, "summary.DPprior_fit")
    has_method <- !is.null(summ$method)
    has_gamma <- !is.null(summ$gamma_prior)
    has_alpha <- !is.null(summ$alpha_summary)
    has_errors <- !is.null(summ$errors)

    all(has_class, has_method, has_gamma, has_alpha, has_errors)
  }, error = function(e) FALSE)

  if (isTRUE(verbose)) {
    cat(sprintf("  summary() works: %s\n", if (test3_pass) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && test3_pass

  # -------------------------------------------------------------------------
  # Test 4: print.summary.DPprior_fit()
  # -------------------------------------------------------------------------
  if (isTRUE(verbose)) {
    cat("\n[Test 4] print.summary.DPprior_fit()\n")
    cat(strrep("-", 50), "\n")
  }

  test4_pass <- tryCatch({
    output <- capture.output(summary(test_fit))
    has_header <- any(grepl("DPprior Prior Elicitation Summary", output))
    has_comparison <- any(grepl("Target vs Achieved", output))
    has_diagnostics <- any(grepl("Diagnostics", output))

    all(has_header, has_comparison, has_diagnostics)
  }, error = function(e) FALSE)

  if (isTRUE(verbose)) {
    cat(sprintf("  print(summary()) works: %s\n", if (test4_pass) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && test4_pass

  # -------------------------------------------------------------------------
  # Test 5: .dpprior_is_dual() helper (stricter check)
  # -------------------------------------------------------------------------
  if (isTRUE(verbose)) {
    cat("\n[Test 5] .dpprior_is_dual() helper (stricter check)\n")
    cat(strrep("-", 50), "\n")
  }

  test5_pass <- tryCatch({
    # K-only fit should not be dual
    not_dual <- !.dpprior_is_dual(test_fit)

    # Dual-anchor fit WITH init should be dual
    test_dual <- test_fit
    test_dual$method <- "dual-anchor"
    test_dual$dual_anchor <- list(lambda = 0.5, init = list(a = 1.6, b = 1.2))
    is_dual <- .dpprior_is_dual(test_dual)

    # Dual-anchor WITHOUT init should NOT be dual (stricter check)
    test_partial <- test_fit
    test_partial$dual_anchor <- list(lambda = 0.5)  # No init
    not_dual_partial <- !.dpprior_is_dual(test_partial)

    not_dual && is_dual && not_dual_partial
  }, error = function(e) FALSE)

  if (isTRUE(verbose)) {
    cat(sprintf("  .dpprior_is_dual() works (strict): %s\n", if (test5_pass) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && test5_pass

  # -------------------------------------------------------------------------
  # Test 6: var_K projection display
  # -------------------------------------------------------------------------
  if (isTRUE(verbose)) {
    cat("\n[Test 6] var_K projection display\n")
    cat(strrep("-", 50), "\n")
  }

  test6_pass <- tryCatch({
    test_projected <- test_fit
    test_projected$target$var_K <- 3.0  # Original was infeasible
    test_projected$target$var_K_used <- 4.0  # Projected to feasible

    output <- capture.output(print(test_projected))
    has_projection <- any(grepl("used after projection", output))

    has_projection
  }, error = function(e) FALSE)

  if (isTRUE(verbose)) {
    cat(sprintf("  var_K projection shown: %s\n", if (test6_pass) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && test6_pass

  # -------------------------------------------------------------------------
  # Test 7: summary with diagnostics parameter
  # -------------------------------------------------------------------------
  if (isTRUE(verbose)) {
    cat("\n[Test 7] summary with diagnostics parameter\n")
    cat(strrep("-", 50), "\n")
  }

  test7_pass <- tryCatch({
    summ <- summary(test_fit, print_output = FALSE)

    # Default: should show reminder, not full diagnostics
    output_default <- capture.output(print(summ))
    has_reminder <- any(grepl("diagnostics=TRUE", output_default))

    # With diagnostics = TRUE: should print full diagnostics
    output_full <- capture.output(print(summ, diagnostics = TRUE))
    has_full <- any(grepl("Full Diagnostics", output_full))

    has_reminder && has_full
  }, error = function(e) FALSE)

  if (isTRUE(verbose)) {
    cat(sprintf("  diagnostics parameter works: %s\n", if (test7_pass) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && test7_pass

  # -------------------------------------------------------------------------
  # Test 8: Output without diagnostics
  # -------------------------------------------------------------------------
  if (isTRUE(verbose)) {
    cat("\n[Test 8] Output without diagnostics\n")
    cat(strrep("-", 50), "\n")
  }

  test6_pass <- tryCatch({
    test_no_diag <- test_fit
    test_no_diag$diagnostics <- NULL

    output <- capture.output(print(test_no_diag))
    # Should still work, just no dominance risk line
    has_header <- any(grepl("DPprior Prior Elicitation Result", output))
    has_gamma <- any(grepl("Gamma Hyperprior", output))

    has_header && has_gamma
  }, error = function(e) FALSE)

  if (isTRUE(verbose)) {
    cat(sprintf("  Works without diagnostics: %s\n", if (test6_pass) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && test6_pass

  # -------------------------------------------------------------------------
  # Summary
  # -------------------------------------------------------------------------
  if (isTRUE(verbose)) {
    cat("\n", strrep("=", 60), "\n", sep = "")
    if (all_pass) {
      cat("All S3 methods tests PASSED!\n")
    } else {
      cat("Some tests FAILED. Check output above.\n")
    }
    cat(strrep("=", 60), "\n", sep = "")
  }

  invisible(all_pass)
}

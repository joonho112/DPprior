# =============================================================================
# Module 15: Dual-Anchor Framework (FINAL CORRECTED VERSION)
# =============================================================================
#
# Author: JoonHo Lee
# Date: December 2025
# Part of: DPprior R Package
#
# CRITICAL FIX (2025-12-29):
# Original implementation used ABSOLUTE squared errors for L_K, creating
# severe scale mismatch where L_K could be 100x larger than L_w.
# Result: optimizer stayed near K-only solution even with lambda=0.5.
#
# SOLUTION: Three loss types are now available:
#   1. "relative" - Normalized squared errors (RECOMMENDED)
#   2. "adaptive" - Scaled by initial loss magnitudes (AGGRESSIVE)
#   3. "absolute" - Legacy behavior (NOT RECOMMENDED)
#
# =============================================================================


# =============================================================================
# Dual-Anchor Loss Function (CORRECTED)
# =============================================================================

#' Dual-Anchor Loss Function
#'
#' Computes the combined loss for K_J and first-weight targets used in the dual-anchor
#' calibration framework.
#'
#' @param a Numeric; shape parameter of the Gamma prior (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior (b > 0).
#' @param J Integer; sample size.
#' @param K_target List with components \code{mu_K} and \code{var_K}.
#' @param w1_target List specifying the first-weight target (quantile, prob, or mean).
#' @param lambda Numeric value between 0 and 1; weight on K_J loss component.
#' @param M Integer; number of Gauss-Laguerre quadrature nodes.
#' @param loss_type Character; type of loss scaling:
#'   \describe{
#'     \item{"relative"}{(DEFAULT) Normalized squared errors - both losses
#'       dimensionless and comparable. RECOMMENDED for most cases.}
#'     \item{"adaptive"}{Scaled by initial loss magnitudes - gives more
#'       aggressive weight reduction. Use when lambda=0.5 should truly
#'       mean "equal importance".}
#'     \item{"absolute"}{Legacy behavior - NOT RECOMMENDED, causes optimizer
#'       to essentially ignore weight constraint.}
#'   }
#' @param L_K_scale,L_w_scale Numeric; scaling factors for adaptive loss.
#'   Only used when \code{loss_type = "adaptive"}. If NULL, computed internally.
#'
#' @return Numeric; the combined loss value.
#'
#' @details
#' ## Loss Scaling Comparison
#'
#' At typical K-only solution with \eqn{P(w_1 > 0.5)} approximately 0.48, target = 0.30:
#' \itemize{
#'   \item \strong{ABSOLUTE}: L_K approximately 0, L_w approximately 0.03. Any move causes L_K >> L_w.
#'     Optimizer stays at K-only. BROKEN.
#'   \item \strong{RELATIVE}: L_K_rel approximately 0, L_w approximately 0.03. Both dimensionless.
#'     lambda = 0.5 gives approximately 9 percent reduction in \eqn{P(w_1 > 0.5)}.
#'   \item \strong{ADAPTIVE}: Losses normalized by initial magnitudes.
#'     lambda = 0.5 gives approximately 18 percent reduction in \eqn{P(w_1 > 0.5)}.
#' }
#'
#' @seealso \code{\link{DPprior_dual}}, \code{\link{compute_tradeoff_curve}}
#'
#' @keywords internal
dual_anchor_loss <- function(a, b, J, K_target, w1_target, lambda,
                             M = .QUAD_NODES_DEFAULT,
                             loss_type = c("relative", "adaptive", "absolute"),
                             L_K_scale = NULL, L_w_scale = NULL) {

  loss_type <- match.arg(loss_type)

  # Boundary check
  if (a <= 0 || b <= 0 || !is.finite(a) || !is.finite(b)) {
    return(.PENALTY_INF)
  }

  # -------------------------------------------------------------------------
  # K_J Component Loss
  # -------------------------------------------------------------------------
  moments <- tryCatch(
    exact_K_moments(J, a, b, M),
    error = function(e) NULL
  )

  if (is.null(moments)) {
    return(.PENALTY_INF)
  }

  if (loss_type == "absolute") {
    # Legacy: Absolute squared errors (NOT RECOMMENDED)
    L_K <- (moments$mean - K_target$mu_K)^2 +
      (moments$var - K_target$var_K)^2
  } else {
    # CORRECTED: Relative (normalized) squared errors
    L_K <- ((moments$mean - K_target$mu_K) / K_target$mu_K)^2 +
      ((moments$var - K_target$var_K) / K_target$var_K)^2
  }

  # -------------------------------------------------------------------------
  # w1 Component Loss
  # -------------------------------------------------------------------------
  L_w <- tryCatch({
    if (!is.null(w1_target$quantile)) {
      u <- w1_target$quantile$prob
      target_val <- w1_target$quantile$value
      actual_val <- quantile_w1(u, a, b)
      (actual_val - target_val)^2

    } else if (!is.null(w1_target$prob)) {
      t <- w1_target$prob$threshold
      target_prob <- w1_target$prob$value
      actual_prob <- prob_w1_exceeds(t, a, b)
      (actual_prob - target_prob)^2

    } else if (!is.null(w1_target$mean)) {
      target_mean <- w1_target$mean
      actual_mean <- mean_w1(a, b, M)
      (actual_mean - target_mean)^2

    } else {
      stop("w1_target must specify quantile, prob, or mean")
    }
  }, error = function(e) .PENALTY_INF)

  if (!is.finite(L_w)) {
    return(.PENALTY_INF)
  }

  # -------------------------------------------------------------------------
  # Apply scaling for adaptive loss type
  # -------------------------------------------------------------------------
  if (loss_type == "adaptive") {
    if (!is.null(L_K_scale) && L_K_scale > 0) {
      L_K <- L_K / L_K_scale
    }
    if (!is.null(L_w_scale) && L_w_scale > 0) {
      L_w <- L_w / L_w_scale
    }
  }

  # -------------------------------------------------------------------------
  # Combined Loss
  # -------------------------------------------------------------------------
  lambda * L_K + (1 - lambda) * L_w
}


# =============================================================================
# Dual-Anchor Calibration (CORRECTED)
# =============================================================================

#' Dual-Anchor Prior Calibration
#'
#' Refines a K_J-calibrated prior to also satisfy weight constraints,
#' implementing the dual-anchor framework from Lee (2026, Section 4).
#'
#' @param fit A \code{DPprior_fit} object from any K-based calibration method.
#' @param w1_target List specifying the first-weight target. Options:
#'   \describe{
#'     \item{\code{list(prob = list(threshold = 0.5, value = 0.3))}}{
#'       Constrain \eqn{P(w_1 > 0.5) = 0.3}}
#'     \item{\code{list(quantile = list(prob = 0.9, value = 0.4))}}{
#'       Constrain 90th percentile of \eqn{w_1} = 0.4}
#'     \item{\code{list(mean = 0.3)}}{Constrain \eqn{E[w_1] = 0.3}}
#'   }
#' @param lambda Numeric value between 0 and 1; weight on K_J anchor.
#'   Default is 0.5 (balanced). lambda = 1 recovers K-only fit.
#' @param max_iter Integer; maximum optimization iterations.
#' @param verbose Logical; if TRUE, print progress.
#' @param M Integer; quadrature nodes.
#' @param loss_type Character; "relative" (default), "adaptive", or "absolute".
#'   See \code{\link{dual_anchor_loss}} for details.
#'
#' @return An S3 object of class \code{"DPprior_fit"}.
#'
#' @details
#' ## Critical Fix (2025-12-29)
#'
#' The original implementation used absolute squared errors for L_K, causing
#' severe scale mismatch. The optimizer would essentially ignore the weight
#' constraint. Example improvement at lambda = 0.5:
#' \itemize{
#'   \item ABSOLUTE (broken): 0.5 percent reduction in \eqn{P(w_1 > 0.5)}
#'   \item RELATIVE (fixed): 9 percent reduction
#'   \item ADAPTIVE: 18 percent reduction
#' }
#'
#' @section Choosing loss_type:
#' \itemize{
#'   \item \strong{"relative"}: Good default. Both losses dimensionless.
#'   \item \strong{"adaptive"}: More aggressive. Use when you want lambda = 0.5
#'     to give significant weight reduction.
#'   \item \strong{"absolute"}: Legacy, not recommended.
#' }
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @examples
#' # K-only fit
#' fit_K <- DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8)
#' cat("K-only P(w_1 > 0.5):", prob_w1_exceeds(0.5, fit_K$a, fit_K$b), "\n")
#'
#' # Dual-anchor with relative loss (default)
#' fit_rel <- DPprior_dual(
#'   fit_K,
#'   w1_target = list(prob = list(threshold = 0.5, value = 0.3)),
#'   lambda = 0.5,
#'   loss_type = "relative"
#' )
#' cat("Relative P(w_1 > 0.5):", fit_rel$dual_anchor$w1_achieved$prob_gt_50, "\n")
#'
#' # Dual-anchor with adaptive loss (more aggressive)
#' fit_adp <- DPprior_dual(
#'   fit_K,
#'   w1_target = list(prob = list(threshold = 0.5, value = 0.3)),
#'   lambda = 0.5,
#'   loss_type = "adaptive"
#' )
#' cat("Adaptive P(w_1 > 0.5):", fit_adp$dual_anchor$w1_achieved$prob_gt_50, "\n")
#'
#' @family elicitation
#'
#' @export
DPprior_dual <- function(fit, w1_target, lambda = 0.5,
                         max_iter = 100L, verbose = FALSE,
                         M = .QUAD_NODES_DEFAULT,
                         loss_type = c("relative", "adaptive", "absolute")) {

  loss_type <- match.arg(loss_type)

  # -------------------------------------------------------------------------
  # Input Validation
  # -------------------------------------------------------------------------
  if (!inherits(fit, "DPprior_fit")) {
    stop("fit must be a DPprior_fit object", call. = FALSE)
  }

  if (!is.numeric(lambda) || length(lambda) != 1L ||
      lambda < 0 || lambda > 1) {
    stop("lambda must be a numeric value in [0, 1]", call. = FALSE)
  }

  # Validate w1_target
  has_quantile <- !is.null(w1_target$quantile)
  has_prob <- !is.null(w1_target$prob)
  has_mean <- !is.null(w1_target$mean)

  if (sum(has_quantile, has_prob, has_mean) != 1L) {
    stop("w1_target must specify exactly one of: quantile, prob, or mean",
         call. = FALSE)
  }

  # -------------------------------------------------------------------------
  # Extract Initial Values
  # -------------------------------------------------------------------------
  J <- fit$J
  K_target <- list(mu_K = fit$target$mu_K, var_K = fit$target$var_K)
  a_init <- fit$a
  b_init <- fit$b

  # -------------------------------------------------------------------------
  # Short-circuit: lambda = 1
  # -------------------------------------------------------------------------
  if (lambda == 1) {
    result <- fit
    result$method <- "dual-anchor"
    result$dual_anchor <- list(
      w1_target = w1_target,
      lambda = lambda,
      loss_type = loss_type,
      w1_achieved = list(
        mean = mean_w1(fit$a, fit$b, M),
        prob_gt_50 = prob_w1_exceeds(0.5, fit$a, fit$b),
        prob_gt_90 = prob_w1_exceeds(0.9, fit$a, fit$b)
      ),
      K_loss = 0,
      init = list(a = a_init, b = b_init),
      note = "lambda = 1 returns the K-only solution"
    )
    return(result)
  }

  # -------------------------------------------------------------------------
  # Compute Scaling Factors for Adaptive Loss
  # -------------------------------------------------------------------------
  L_K_scale <- NULL
  L_w_scale <- NULL

  if (loss_type == "adaptive") {
    # -----------------------------------------------------------------------
    # Adaptive Scaling: Estimate scale factors from both anchor extremes
    # -----------------------------------------------------------------------
    # At K-only solution: L_K ≈ 0, L_w = L_w_init
    # At w-only solution: L_K = L_K_at_w, L_w ≈ 0
    #
    # For meaningful lambda trade-off, we normalize by these scales:
    #   L_normalized = lambda * (L_K / L_K_scale) + (1-lambda) * (L_w / L_w_scale)
    # -----------------------------------------------------------------------

    # Step 1: Get L_w at K-only point (this is our w-anchor scale)
    L_w_init <- dual_anchor_loss(a_init, b_init, J, K_target, w1_target,
                                 0.0, M, "relative")

    # Step 2: Quick optimization to estimate w-only solution and its K_loss
    # This gives us the K-anchor scale
    obj_w_only <- function(eta) {
      a_curr <- exp(eta[1L])
      b_curr <- exp(eta[2L])
      dual_anchor_loss(a_curr, b_curr, J, K_target, w1_target,
                       0.0, M, "relative")  # lambda=0 means w-only
    }

    opt_w <- optim(
      par = c(log(a_init), log(b_init)),
      fn = obj_w_only,
      method = "BFGS",
      control = list(maxit = 50)  # Quick estimate
    )

    a_w <- exp(opt_w$par[1L])
    b_w <- exp(opt_w$par[2L])

    # K_loss at the w-only solution
    L_K_at_w <- dual_anchor_loss(a_w, b_w, J, K_target, w1_target,
                                 1.0, M, "relative")

    # Set scale factors
    L_K_scale <- max(L_K_at_w, 0.01)   # Scale by K-loss at w-only extreme
    L_w_scale <- max(L_w_init, 0.001)  # Scale by w-loss at K-only extreme

    if (isTRUE(verbose)) {
      cat(sprintf("Adaptive scaling:\n"))
      cat(sprintf("  L_w at K-only: %.4e\n", L_w_init))
      cat(sprintf("  L_K at w-only: %.4e (a=%.3f, b=%.3f)\n", L_K_at_w, a_w, b_w))
      cat(sprintf("  L_K_scale = %.4e, L_w_scale = %.4e\n", L_K_scale, L_w_scale))
    }
  }

  # -------------------------------------------------------------------------
  # Initial Diagnostics
  # -------------------------------------------------------------------------
  if (isTRUE(verbose)) {
    init_L_K <- dual_anchor_loss(a_init, b_init, J, K_target, w1_target,
                                 1.0, M, loss_type, L_K_scale, L_w_scale)
    init_L_w <- dual_anchor_loss(a_init, b_init, J, K_target, w1_target,
                                 0.0, M, loss_type, L_K_scale, L_w_scale)

    cat(sprintf("Starting dual-anchor (lambda=%.2f, loss_type=%s)\n",
                lambda, loss_type))
    cat(sprintf("  Initial: a=%.4f, b=%.4f\n", a_init, b_init))
    cat(sprintf("  Initial losses: L_K=%.4e, L_w=%.4e\n", init_L_K, init_L_w))
  }

  # -------------------------------------------------------------------------
  # Optimization
  # -------------------------------------------------------------------------
  objective <- function(eta) {
    eta <- pmin(pmax(eta, -.EXP_MAX), .EXP_MAX)
    a_curr <- exp(eta[1L])
    b_curr <- exp(eta[2L])
    dual_anchor_loss(a_curr, b_curr, J, K_target, w1_target, lambda, M,
                     loss_type, L_K_scale, L_w_scale)
  }

  opt <- optim(
    par = c(log(a_init), log(b_init)),
    fn = objective,
    method = "BFGS",
    control = list(maxit = max_iter)
  )

  # Fallback: Nelder-Mead
  if (opt$convergence != 0) {
    if (isTRUE(verbose)) {
      cat("  BFGS did not converge, trying Nelder-Mead...\n")
    }
    opt_nm <- optim(
      par = c(log(a_init), log(b_init)),
      fn = objective,
      method = "Nelder-Mead",
      control = list(maxit = max_iter * 2)
    )
    if (opt_nm$value < opt$value) {
      opt <- opt_nm
    }
  }

  a_opt <- exp(opt$par[1L])
  b_opt <- exp(opt$par[2L])

  # -------------------------------------------------------------------------
  # Final Diagnostics
  # -------------------------------------------------------------------------
  moments <- exact_K_moments(J, a_opt, b_opt, M)

  final_L_K <- dual_anchor_loss(a_opt, b_opt, J, K_target, w1_target,
                                1.0, M, loss_type, L_K_scale, L_w_scale)
  final_L_w <- dual_anchor_loss(a_opt, b_opt, J, K_target, w1_target,
                                0.0, M, loss_type, L_K_scale, L_w_scale)

  w1_achieved <- list(
    mean = mean_w1(a_opt, b_opt, M),
    prob_gt_50 = prob_w1_exceeds(0.5, a_opt, b_opt),
    prob_gt_90 = prob_w1_exceeds(0.9, a_opt, b_opt)
  )

  if (isTRUE(verbose)) {
    cat(sprintf("  Final: a=%.4f, b=%.4f\n", a_opt, b_opt))
    cat(sprintf("  Achieved: mu_K=%.4f, var_K=%.4f\n", moments$mean, moments$var))
    cat(sprintf("  Final losses: L_K=%.4e, L_w=%.4e\n", final_L_K, final_L_w))
    cat(sprintf("  P(w1 > 0.5) = %.4f\n", w1_achieved$prob_gt_50))
  }

  # -------------------------------------------------------------------------
  # Build Result
  # -------------------------------------------------------------------------
  result <- list(
    a = a_opt,
    b = b_opt,
    J = J,
    target = fit$target,
    method = "dual-anchor",
    converged = (opt$convergence == 0),
    iterations = if (!is.null(opt$counts[["function"]])) opt$counts[["function"]] else NA_integer_,
    fit = list(
      mu_K = moments$mean,
      var_K = moments$var,
      residual = sqrt((moments$mean - K_target$mu_K)^2 +
                        (moments$var - K_target$var_K)^2)
    ),
    dual_anchor = list(
      w1_target = w1_target,
      lambda = lambda,
      loss_type = loss_type,
      w1_achieved = w1_achieved,
      K_loss = final_L_K,
      w_loss = final_L_w,
      total_loss = opt$value,
      init = list(a = a_init, b = b_init),
      scaling = if (loss_type == "adaptive")
        list(L_K_scale = L_K_scale, L_w_scale = L_w_scale) else NULL,
      optim = list(
        convergence = opt$convergence,
        value = opt$value,
        counts = opt$counts
      )
    )
  )

  class(result) <- "DPprior_fit"
  result
}


# =============================================================================
# Trade-off Curve
# =============================================================================

#' Compute Pareto Trade-off Curve
#'
#' Explores the trade-off between K_J fit and first-weight constraint over lambda values.
#'
#' @param J Integer; sample size.
#' @param K_target List with \code{mu_K} and \code{var_K}.
#' @param w1_target Weight target specification.
#' @param lambda_seq Numeric vector of lambda values.
#' @param max_iter Integer; max iterations per optimization.
#' @param M Integer; quadrature nodes.
#' @param verbose Logical; print progress.
#' @param loss_type Character; "relative", "adaptive", or "absolute".
#'
#' @return Data frame with columns: lambda, a, b, mu_K, var_K, K_loss, w_loss,
#'   w1_prob_gt_50, E_w1.
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @examples
#' curve <- compute_tradeoff_curve(
#'   J = 50,
#'   K_target = list(mu_K = 5, var_K = 8),
#'   w1_target = list(prob = list(threshold = 0.5, value = 0.25)),
#'   lambda_seq = seq(0, 1, by = 0.1),
#'   loss_type = "relative"
#' )
#'
#' # Visualize trade-off
#' plot(curve$lambda, curve$w1_prob_gt_50, type = "b",
#'      xlab = "lambda", ylab = "P(w_1 > 0.5)",
#'      main = "Weight Dominance vs Lambda")
#' abline(h = 0.25, lty = 2, col = "red")  # target
#'
#' @export
compute_tradeoff_curve <- function(J, K_target, w1_target,
                                   lambda_seq = seq(0, 1, by = 0.1),
                                   max_iter = 100L,
                                   M = .QUAD_NODES_DEFAULT,
                                   verbose = FALSE,
                                   loss_type = c("relative", "adaptive", "absolute")) {

  loss_type <- match.arg(loss_type)

  if (is.null(K_target$mu_K) || is.null(K_target$var_K)) {
    stop("K_target must contain mu_K and var_K", call. = FALSE)
  }

  # K-only fit using DPprior_a2_newton
  init_fit <- DPprior_a2_newton(J, K_target$mu_K, K_target$var_K,
                                verbose = FALSE, M = M)
  # Ensure required fields exist for DPprior_dual
  if (is.null(init_fit$J)) init_fit$J <- J
  if (is.null(init_fit$target)) {
    init_fit$target <- list(mu_K = K_target$mu_K, var_K = K_target$var_K, type = "moments")
  }
  if (!inherits(init_fit, "DPprior_fit")) {
    class(init_fit) <- "DPprior_fit"
  }

  results <- lapply(lambda_seq, function(lam) {
    if (isTRUE(verbose)) {
      cat(sprintf("Computing lambda = %.3f\n", lam))
    }

    if (lam == 1) {
      fit <- init_fit
    } else {
      fit <- DPprior_dual(init_fit, w1_target = w1_target,
                          lambda = lam, max_iter = max_iter,
                          M = M, verbose = FALSE, loss_type = loss_type)
    }

    moments <- exact_K_moments(J, fit$a, fit$b, M)

    # Compute standardized losses for comparison
    K_loss <- dual_anchor_loss(fit$a, fit$b, J, K_target, w1_target,
                               1.0, M, "relative")
    w_loss <- dual_anchor_loss(fit$a, fit$b, J, K_target, w1_target,
                               0.0, M, "relative")

    data.frame(
      lambda = lam,
      a = fit$a,
      b = fit$b,
      mu_K = moments$mean,
      var_K = moments$var,
      K_loss = K_loss,
      w_loss = w_loss,
      w1_prob_gt_50 = prob_w1_exceeds(0.5, fit$a, fit$b),
      E_w1 = mean_w1(fit$a, fit$b, M)
    )
  })

  do.call(rbind, results)
}


# =============================================================================
# Diagnostic Helper
# =============================================================================

#' Dual-Anchor Diagnostic Comparison
#'
#' @param fit_dual Dual-anchor DPprior_fit object.
#' @param fit_K_only Optional K-only DPprior_fit object.
#' @param M Integer; quadrature nodes.
#'
#' @return Data frame comparing metrics.
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @family diagnostics
#'
#' @export
dual_anchor_diagnostics <- function(fit_dual, fit_K_only = NULL,
                                    M = .QUAD_NODES_DEFAULT) {

  if (!inherits(fit_dual, "DPprior_fit")) {
    stop("fit_dual must be a DPprior_fit object", call. = FALSE)
  }

  # Get K-only parameters
  if (is.null(fit_K_only) && !is.null(fit_dual$dual_anchor$init)) {
    a_K <- fit_dual$dual_anchor$init$a
    b_K <- fit_dual$dual_anchor$init$b
  } else if (!is.null(fit_K_only)) {
    a_K <- fit_K_only$a
    b_K <- fit_K_only$b
  } else {
    stop("Need fit_K_only or fit_dual$dual_anchor$init", call. = FALSE)
  }

  J <- fit_dual$J

  moments_K <- exact_K_moments(J, a_K, b_K, M)
  moments_dual <- exact_K_moments(J, fit_dual$a, fit_dual$b, M)

  data.frame(
    Metric = c("a", "b", "E[K_J]", "Var[K_J]",
               "P(w1 > 0.5)", "P(w1 > 0.9)", "E[w1]"),
    K_only = c(
      a_K, b_K, moments_K$mean, moments_K$var,
      prob_w1_exceeds(0.5, a_K, b_K),
      prob_w1_exceeds(0.9, a_K, b_K),
      mean_w1(a_K, b_K, M)
    ),
    Dual_anchor = c(
      fit_dual$a, fit_dual$b, moments_dual$mean, moments_dual$var,
      prob_w1_exceeds(0.5, fit_dual$a, fit_dual$b),
      prob_w1_exceeds(0.9, fit_dual$a, fit_dual$b),
      mean_w1(fit_dual$a, fit_dual$b, M)
    )
  )
}


# =============================================================================
# Verification Tests
# =============================================================================

#' Verify Dual-Anchor Module
#'
#' @param verbose Logical; print detailed output.
#' @return Invisible TRUE if all tests pass.
#'
#' @examples
#' \dontrun{
#' verify_dual_anchor(verbose = TRUE)
#'
#' }
#' @keywords internal
verify_dual_anchor <- function(verbose = TRUE) {

  if (isTRUE(verbose)) {
    cat("=== Dual-Anchor Verification Tests ===\n\n")
  }

  J <- 50
  mu_K <- 5
  var_K <- 8

  # K-only fit using DPprior_a2_newton
  fit_K <- DPprior_a2_newton(J, mu_K, var_K, verbose = FALSE)
  # Ensure required fields exist for DPprior_dual
  if (is.null(fit_K$J)) fit_K$J <- J
  if (is.null(fit_K$target)) {
    fit_K$target <- list(mu_K = mu_K, var_K = var_K)
  }
  if (!inherits(fit_K, "DPprior_fit")) {
    class(fit_K) <- "DPprior_fit"
  }

  K_target <- list(mu_K = mu_K, var_K = var_K)
  w1_target <- list(prob = list(threshold = 0.5, value = 0.30))
  p_before <- prob_w1_exceeds(0.5, fit_K$a, fit_K$b)

  # -------------------------------------------------------------------------
  # Test 1: ABSOLUTE loss fails (baseline)
  # -------------------------------------------------------------------------
  if (isTRUE(verbose)) cat("Test 1: ABSOLUTE loss barely moves (expected)\n")

  fit_abs <- DPprior_dual(fit_K, w1_target, lambda = 0.5,
                          verbose = FALSE, loss_type = "absolute")
  p_abs <- fit_abs$dual_anchor$w1_achieved$prob_gt_50
  red_abs <- (p_before - p_abs) / p_before

  if (isTRUE(verbose)) {
    cat(sprintf("  Before: %.4f, After: %.4f, Reduction: %.1f%%\n",
                p_before, p_abs, red_abs * 100))
  }

  stopifnot(red_abs < 0.02)  # Should be < 2% reduction (broken)
  if (isTRUE(verbose)) cat("  PASS (confirmed ABSOLUTE loss is inadequate)\n\n")

  # -------------------------------------------------------------------------
  # Test 2: RELATIVE loss improves significantly
  # -------------------------------------------------------------------------
  if (isTRUE(verbose)) cat("Test 2: RELATIVE loss gives meaningful reduction\n")

  fit_rel <- DPprior_dual(fit_K, w1_target, lambda = 0.5,
                          verbose = FALSE, loss_type = "relative")
  p_rel <- fit_rel$dual_anchor$w1_achieved$prob_gt_50
  red_rel <- (p_before - p_rel) / p_before

  if (isTRUE(verbose)) {
    cat(sprintf("  Before: %.4f, After: %.4f, Reduction: %.1f%%\n",
                p_before, p_rel, red_rel * 100))
  }

  stopifnot(red_rel > 0.05)  # Should be > 5% reduction
  if (isTRUE(verbose)) cat("  PASS\n\n")

  # -------------------------------------------------------------------------
  # Test 3: ADAPTIVE loss is most aggressive
  # -------------------------------------------------------------------------
  if (isTRUE(verbose)) cat("Test 3: ADAPTIVE loss is more aggressive\n")

  fit_adp <- DPprior_dual(fit_K, w1_target, lambda = 0.5,
                          verbose = FALSE, loss_type = "adaptive")
  p_adp <- fit_adp$dual_anchor$w1_achieved$prob_gt_50
  red_adp <- (p_before - p_adp) / p_before

  if (isTRUE(verbose)) {
    cat(sprintf("  Before: %.4f, After: %.4f, Reduction: %.1f%%\n",
                p_before, p_adp, red_adp * 100))
  }

  stopifnot(red_adp > red_rel)  # Adaptive should beat relative
  if (isTRUE(verbose)) cat("  PASS\n\n")

  # -------------------------------------------------------------------------
  # Test 4: lambda = 1 recovers K-only
  # -------------------------------------------------------------------------
  if (isTRUE(verbose)) cat("Test 4: lambda = 1 recovers K-only\n")

  fit_1 <- DPprior_dual(fit_K, w1_target, lambda = 1.0, verbose = FALSE)

  stopifnot(abs(fit_1$a - fit_K$a) < 0.01)
  stopifnot(abs(fit_1$b - fit_K$b) < 0.01)

  if (isTRUE(verbose)) cat("  PASS\n\n")

  # -------------------------------------------------------------------------
  # Test 5: lambda = 0 achieves weight target
  # -------------------------------------------------------------------------
  if (isTRUE(verbose)) cat("Test 5: lambda = 0 achieves weight target\n")

  fit_0 <- DPprior_dual(fit_K, w1_target, lambda = 0.0,
                        verbose = FALSE, loss_type = "relative")
  p_0 <- fit_0$dual_anchor$w1_achieved$prob_gt_50

  if (isTRUE(verbose)) {
    cat(sprintf("  Target: 0.30, Achieved: %.4f\n", p_0))
  }

  stopifnot(abs(p_0 - 0.30) < 0.02)  # Should be very close to target
  if (isTRUE(verbose)) cat("  PASS\n\n")

  # -------------------------------------------------------------------------
  # Test 6: Trade-off curve monotonicity
  # -------------------------------------------------------------------------
  if (isTRUE(verbose)) cat("Test 6: Trade-off curve\n")

  curve <- compute_tradeoff_curve(
    J = J, K_target = K_target,
    w1_target = w1_target,
    lambda_seq = c(0, 0.5, 1.0),
    verbose = FALSE, loss_type = "relative"
  )

  if (isTRUE(verbose)) {
    print(curve[, c("lambda", "w1_prob_gt_50", "mu_K", "var_K")])
  }

  stopifnot(nrow(curve) == 3)
  stopifnot(curve$w1_prob_gt_50[1] < curve$w1_prob_gt_50[3])  # lambda=0 < lambda=1

  if (isTRUE(verbose)) cat("  PASS\n\n")

  if (isTRUE(verbose)) {
    cat("=== All Tests Passed ===\n")
    cat("\nSummary:\n")
    cat(sprintf("  ABSOLUTE (broken): %.1f%% reduction\n", red_abs * 100))
    cat(sprintf("  RELATIVE (fixed):  %.1f%% reduction\n", red_rel * 100))
    cat(sprintf("  ADAPTIVE (best):   %.1f%% reduction\n", red_adp * 100))
  }

  invisible(TRUE)
}

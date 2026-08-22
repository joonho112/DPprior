# =============================================================================
# Module 02: Gauss-Laguerre Quadrature for Gamma Distributions
# =============================================================================
#
# This module provides numerical integration against Gamma distributions
# using generalized Gauss-Laguerre quadrature.
#
# Theory Background:
# -----------------
# Standard Laguerre quadrature integrates:
#   int_0^inf f(x) x^beta e^{-x} dx = sum w_m f(x_m)
#
# For Gamma(a, b) prior expectations E[g(alpha)]:
#   E[g(alpha)] = int_0^inf g(alpha) p(alpha|a,b) d_alpha
#
# Using change of variables x = b*alpha:
#   E[g(alpha)] = (1/Gamma(a)) int_0^inf g(x/b) x^{a-1} e^{-x} dx
#
# Generalized Gauss-Laguerre quadrature with parameter alpha_param = a - 1
# provides nodes and weights for exact integration of polynomial * weight.
#
# Author: JoonHo Lee (jlee296@ua.edu)
# Date: December 2025
# Part of: DPprior R Package
# Reference: Lee (2026), Sections 2--3
# =============================================================================


# Dense symmetric eigendecomposition has O(M^2) memory and O(M^3) time.
# Keep the supported node count and cache policy local to this implementation;
# callers continue to use the existing public APIs.
.QUADRATURE_MAX_NODES <- 512L
.QUADRATURE_CACHE_MAX_ENTRIES <- 16L


# A verification computation must be materially higher order than the selected
# rule. Keep this formula in one place so generic quadrature and all marginal
# distribution callers enforce the same independence contract.
.quadrature_verification_required_order <- function(M) {
  as.integer(max(2 * M, M + 40L))
}

.dpprior_quadrature_cache <- local({
  cache <- new.env(parent = emptyenv())
  cache$entries <- new.env(parent = emptyenv())
  cache$order <- character()
  cache$hits <- 0L
  cache$misses <- 0L
  cache$evictions <- 0L
  cache
})


.quadrature_cache_key <- function(M, alpha_param) {
  if (identical(alpha_param, -0)) alpha_param <- 0
  sprintf("M=%d|alpha=%.17g", as.integer(M), alpha_param)
}


.clear_quadrature_cache <- function() {
  keys <- ls(.dpprior_quadrature_cache$entries, all.names = TRUE)
  if (length(keys)) {
    rm(list = keys, envir = .dpprior_quadrature_cache$entries)
  }
  .dpprior_quadrature_cache$order <- character()
  .dpprior_quadrature_cache$hits <- 0L
  .dpprior_quadrature_cache$misses <- 0L
  .dpprior_quadrature_cache$evictions <- 0L
  invisible(TRUE)
}


.quadrature_cache_info <- function() {
  keys <- ls(.dpprior_quadrature_cache$entries, all.names = TRUE)
  bytes <- if (!length(keys)) {
    0
  } else {
    sum(vapply(keys, function(key) {
      as.numeric(utils::object.size(get(
        key,
        envir = .dpprior_quadrature_cache$entries,
        inherits = FALSE
      )))
    }, numeric(1)))
  }

  list(
    policy = "LRU",
    max_entries = .QUADRATURE_CACHE_MAX_ENTRIES,
    entries = length(keys),
    keys = .dpprior_quadrature_cache$order,
    hits = .dpprior_quadrature_cache$hits,
    misses = .dpprior_quadrature_cache$misses,
    evictions = .dpprior_quadrature_cache$evictions,
    bytes = bytes
  )
}


.quadrature_cache_get <- function(key) {
  if (!exists(key, envir = .dpprior_quadrature_cache$entries,
              inherits = FALSE)) {
    .dpprior_quadrature_cache$misses <-
      .dpprior_quadrature_cache$misses + 1L
    return(NULL)
  }

  .dpprior_quadrature_cache$hits <- .dpprior_quadrature_cache$hits + 1L
  .dpprior_quadrature_cache$order <- c(
    setdiff(.dpprior_quadrature_cache$order, key), key
  )
  result <- get(key, envir = .dpprior_quadrature_cache$entries,
                inherits = FALSE)
  metadata <- attr(result, "quadrature_metadata", exact = TRUE)
  metadata$cache_hit <- TRUE
  attr(result, "quadrature_metadata") <- metadata
  result
}


.quadrature_cache_set <- function(key, value) {
  entries <- .dpprior_quadrature_cache$entries
  is_new <- !exists(key, envir = entries, inherits = FALSE)

  if (is_new && length(.dpprior_quadrature_cache$order) >=
      .QUADRATURE_CACHE_MAX_ENTRIES) {
    evict <- .dpprior_quadrature_cache$order[1L]
    rm(list = evict, envir = entries)
    .dpprior_quadrature_cache$order <-
      .dpprior_quadrature_cache$order[-1L]
    .dpprior_quadrature_cache$evictions <-
      .dpprior_quadrature_cache$evictions + 1L
  }

  assign(key, value, envir = entries)
  .dpprior_quadrature_cache$order <- c(
    setdiff(.dpprior_quadrature_cache$order, key), key
  )
  invisible(value)
}


# Retrieve the versioned computation metadata without changing the public
# list components returned by the quadrature constructors.
.quadrature_metadata <- function(x) {
  metadata <- attr(x, "quadrature_metadata", exact = TRUE)
  if (is.null(metadata)) {
    stop("x does not contain DPprior quadrature metadata", call. = FALSE)
  }
  metadata
}


# Evaluate a scalar integrand at quadrature nodes without allowing recycling,
# vector-valued returns, or non-finite values to enter a weighted sum.
.quadrature_evaluate_integrand <- function(f, nodes) {
  if (!is.function(f)) {
    .dpprior_abort_invalid(
      "f must be a function",
      "dpprior_integrand_error", "f", f, "function", "type"
    )
  }

  vapply(nodes, function(node) {
    .dpprior_validate_scalar(
      f(node), "f(alpha)", .subclass = "dpprior_integrand_error"
    )
  }, numeric(1))
}


# Form a convex quadrature sum after scaling by the largest absolute value.
# This avoids overflow in intermediate partial sums when every nodewise value
# is finite but close to the largest representable double.
.quadrature_weighted_sum <- function(weights, values,
                                     fail_on_nonfinite = TRUE) {
  if (!is.numeric(weights) || !is.numeric(values) ||
      length(weights) != length(values) || length(weights) == 0L ||
      any(!is.finite(weights)) || any(weights < 0) ||
      any(!is.finite(values))) {
    .dpprior_abort_invalid(
      "quadrature weights and values must be equal-length finite numeric vectors with non-negative weights",
      "dpprior_numerical_error", "weights/values", list(weights, values),
      "valid weighted-sum inputs", "invalid_weighted_sum"
    )
  }
  weight_total <- sum(weights)
  if (!is.finite(weight_total) || weight_total <= 0) {
    .dpprior_abort_invalid(
      "quadrature weights must have a finite positive total",
      "dpprior_numerical_error", "weights", weights,
      "finite positive weight total", "invalid_weight_total"
    )
  }
  normalized_weights <- weights / weight_total
  max_abs <- max(abs(values))
  result <- if (max_abs == 0) {
    0
  } else {
    scaled_result <- sum(normalized_weights * (values / max_abs))
    # Exact convexity bounds the scaled result to [-1,1]. Enforce only that
    # structural bound so a 1+epsilon normalized-weight total cannot turn a
    # finite double-maximum constant into a false overflow.
    if (scaled_result > 1 && scaled_result <= 1 + 64 * .Machine$double.eps) {
      scaled_result <- 1
    }
    if (scaled_result < -1 && scaled_result >= -1 - 64 * .Machine$double.eps) {
      scaled_result <- -1
    }
    max_abs * scaled_result
  }

  if (!is.finite(result) && fail_on_nonfinite) {
    .dpprior_abort_invalid(
      "Gamma quadrature produced a non-finite weighted integral",
      "dpprior_numerical_error", "f", values,
      "finite weighted integral", "nonfinite_integral"
    )
  }
  result
}


#' Gauss-Laguerre Quadrature Nodes and Weights
#'
#' Computes generalized Gauss-Laguerre quadrature nodes and weights via
#' eigendecomposition of the Jacobi matrix. These are used to approximate
#' integrals of the form \eqn{\int_0^\infty f(x) x^\beta e^{-x} dx}.
#'
#' @param M Integer; number of quadrature nodes (at least 1, typically 40-120).
#' @param alpha_param Numeric; Laguerre parameter (must be > -1).
#'   For standard Laguerre polynomials, use 0. For integrating against
#'   Gamma(a, b), use \code{alpha_param = a - 1}.
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{nodes}}{Numeric vector of quadrature nodes \eqn{x_m}.}
#'     \item{\code{weights}}{Numeric vector of quadrature weights \eqn{w_m}.}
#'     \item{\code{weights_log}}{Numeric vector of \eqn{\log(w_m)} for
#'       numerical stability.}
#'   }
#'   The list also carries a \code{quadrature_metadata} attribute recording
#'   node count, Laguerre parameter, algorithm, and cache provenance.
#'
#' @details
#' The algorithm constructs a tridiagonal Jacobi matrix \eqn{J} of size
#' \eqn{M \times M}:
#' \itemize{
#'   \item Diagonal: \eqn{a_k = 2k + 1 + \alpha} for \eqn{k = 0, \ldots, M-1}
#'   \item Off-diagonal: \eqn{b_k = \sqrt{k(k + \alpha)}} for \eqn{k = 1, \ldots, M-1}
#' }
#' Eigendecomposition \eqn{J = V D V^T} yields:
#' \itemize{
#'   \item Nodes = eigenvalues (diagonal of \eqn{D})
#'   \item Weights = \eqn{\Gamma(\alpha + 1) \cdot V[1,:]^2}
#' }
#'
#' The generalized Laguerre polynomials \eqn{L_n^{(\alpha)}(x)} are orthogonal
#' with respect to the weight function \eqn{w(x) = x^\alpha e^{-x}} on
#' \eqn{[0, \infty)}.
#'
#' @examples
#' # Standard Laguerre (alpha_param = 0)
#' quad <- gauss_laguerre_nodes(40)
#' sum(quad$weights)  # Should equal Gamma(1) = 1
#'
#' # For Gamma(2.5, b) integration, use alpha_param = 1.5
#' quad <- gauss_laguerre_nodes(80, alpha_param = 1.5)
#' sum(quad$weights)  # Should equal Gamma(2.5)
#'
#' @references
#' Golub, G. H., & Welsch, J. H. (1969). Calculation of Gauss Quadrature Rules.
#' \emph{Mathematics of Computation}, 23(106), 221-230.
#'
#' @seealso \code{\link{build_gamma_quadrature}} for Gamma distribution integration,
#'   \code{\link{integrate_gamma}} for high-level expectation computation
#'
#' @export
gauss_laguerre_nodes <- function(M, alpha_param = 0) {
  # Input validation
  if (!is.numeric(M) || length(M) != 1L || !is.finite(M) ||
      M != floor(M) || M < 1L) {
    .dpprior_abort_legacy_numeric(
      "M must be a positive integer", M, "M", "positive integer",
      "dpprior_control_error", scalar = TRUE, integer = TRUE
    )
  }
  if (!is.numeric(alpha_param) || length(alpha_param) != 1L ||
      !is.finite(alpha_param) || alpha_param <= -1) {
    .dpprior_abort_legacy_numeric(
      "alpha_param must be a finite number > -1", alpha_param,
      "alpha_param", "finite scalar > -1",
      "dpprior_quadrature_parameter_error", scalar = TRUE
    )
  }

  M <- as.integer(M)
  if (M > .QUADRATURE_MAX_NODES) {
    .dpprior_abort_invalid(
      sprintf(
        "M exceeds maximum supported value for dense quadrature (%d)",
        .QUADRATURE_MAX_NODES
      ),
      c("dpprior_control_error", "dpprior_bounds_error"), "M", M,
      sprintf("<= %d", .QUADRATURE_MAX_NODES), "bounds"
    )
  }

  cache_key <- .quadrature_cache_key(M, alpha_param)
  cached <- .quadrature_cache_get(cache_key)
  if (!is.null(cached)) {
    return(cached)
  }

  # Handle M = 1 special case
  if (M == 1L) {
    node <- alpha_param + 1
    weights_norm <- 1.0
    weight_log <- lgamma(alpha_param + 1)
    weight <- exp(weight_log)
    result <- list(
      nodes = node,
      weights = weight,
      weights_log = weight_log,
      weights_norm = weights_norm
    )
    attr(result, "quadrature_metadata") <- list(
      schema_version = 1L,
      engine = "golub-welsch-dense-eigen",
      M_selected = M,
      alpha_param = alpha_param,
      cache_key = cache_key,
      cache_hit = FALSE,
      zero_normalized_weights = 0L,
      normalized_weight_sum = 1
    )
    .quadrature_cache_set(cache_key, result)
    return(result)
  }

  # Construct tridiagonal Jacobi matrix
  # Diagonal elements: a_k = 2k + 1 + alpha for k = 0, ..., M-1
  n <- 0:(M - 1L)
  diag_elem <- 2 * n + alpha_param + 1

  # Off-diagonal elements: b_k = sqrt(k * (k + alpha)) for k = 1, ..., M-1
  k <- 1:(M - 1L)
  off_diag <- sqrt(k * (k + alpha_param))
  if (any(!is.finite(diag_elem)) || any(!is.finite(off_diag))) {
    .dpprior_abort_invalid(
      "Gauss-Laguerre recurrence produced non-finite matrix entries; alpha_param is outside the numerical domain",
      "dpprior_numerical_domain_error", "alpha_param", alpha_param,
      "finite Jacobi-matrix entries", "nonfinite_recurrence"
    )
  }

  # Build symmetric tridiagonal matrix
  J <- diag(diag_elem)
  for (i in seq_len(M - 1L)) {
    J[i, i + 1L] <- off_diag[i]
    J[i + 1L, i] <- off_diag[i]
  }

  # Eigendecomposition (symmetric matrix)
  eig <- tryCatch(
    eigen(J, symmetric = TRUE),
    error = function(error) {
      .dpprior_abort_invalid(
        sprintf(
          "Gauss-Laguerre eigendecomposition failed: %s",
          conditionMessage(error)
        ),
        "dpprior_numerical_error", "alpha_param", alpha_param,
        "stable finite eigendecomposition", "eigendecomposition"
      )
    }
  )

  # Sort by eigenvalues (nodes)
  ord <- order(eig$values)
  nodes <- eig$values[ord]
  V <- eig$vectors[, ord, drop = FALSE]

  # Normalized weights are V[1,]^2 (sum to 1 up to rounding)
  weights_norm <- V[1L, ]^2
  weights_norm_sum <- sum(weights_norm)
  if (!is.finite(weights_norm_sum) || weights_norm_sum <= 0) {
    stop("Gauss-Laguerre normalized weights are non-finite or sum to zero",
         call. = FALSE)
  }
  weights_norm <- weights_norm / weights_norm_sum

  # Log weights (allow -Inf if some weights_norm are numerically zero)
  weights_log <- lgamma(alpha_param + 1) + suppressWarnings(log(weights_norm))

  # Raw weights can overflow for large alpha_param; kept for completeness
  weights <- suppressWarnings(exp(weights_log))

  result <- list(
    nodes = nodes,
    weights = weights,
    weights_log = weights_log,
    weights_norm = weights_norm
  )
  attr(result, "quadrature_metadata") <- list(
    schema_version = 1L,
    engine = "golub-welsch-dense-eigen",
    M_selected = M,
    alpha_param = alpha_param,
    cache_key = cache_key,
    cache_hit = FALSE,
    zero_normalized_weights = sum(weights_norm == 0),
    normalized_weight_sum = sum(weights_norm)
  )
  .quadrature_cache_set(cache_key, result)
  result
}


#' Build Quadrature for Gamma(a, b) Integration
#'
#' Transforms standard Gauss-Laguerre quadrature for integration against
#' a Gamma(a, b) distribution with shape \code{a} and rate \code{b}.
#'
#' @param a Numeric; shape parameter of Gamma distribution (must be > 0).
#' @param b Numeric; rate parameter of Gamma distribution (must be > 0).
#' @param M Integer; number of quadrature nodes (default: 80).
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{a}}{Shape parameter.}
#'     \item{\code{b}}{Rate parameter.}
#'     \item{\code{alpha_nodes}}{Numeric vector of transformed nodes
#'       \eqn{\alpha_m = x_m / b} on the \eqn{\alpha} scale.}
#'     \item{\code{weights_normalized}}{Numeric vector of normalized weights
#'       that sum to 1.}
#'   }
#'   The list also carries a \code{quadrature_metadata} attribute with the
#'   selected \code{M}, Gamma parameterization, node-cache provenance, and
#'   normalization method. List component names are unchanged.
#'
#' @details
#' For \eqn{\alpha \sim \text{Gamma}(a, b)}:
#' \deqn{E[g(\alpha)] = \frac{1}{\Gamma(a)} \int_0^\infty g(x/b) x^{a-1} e^{-x} dx}
#'
#' Using generalized Laguerre quadrature with parameter \eqn{\alpha_{\text{param}} = a - 1}:
#' \deqn{E[g(\alpha)] \approx \sum_{m=1}^M \tilde{w}_m g(\alpha_m)}
#'
#' where \eqn{\alpha_m = x_m / b} and \eqn{\tilde{w}_m} are normalized weights
#' summing to 1.
#'
#' The weights are normalized in log-space for numerical stability. Accuracy
#' as \code{M} changes is integrand-dependent and is not generally monotone.
#'
#' @examples
#' # Build quadrature for Gamma(2.5, 1.5)
#' quad <- build_gamma_quadrature(2.5, 1.5)
#'
#' # Check weights sum to 1
#' sum(quad$weights_normalized)
#'
#' # Check mean: E[alpha] should be a/b = 2.5/1.5
#' sum(quad$weights_normalized * quad$alpha_nodes)
#'
#' @seealso \code{\link{gauss_laguerre_nodes}} for raw quadrature computation,
#'   \code{\link{integrate_gamma}} for high-level expectation computation
#'
#' @export
build_gamma_quadrature <- function(a, b, M = .QUAD_NODES_DEFAULT) {
  # Input validation
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

  # Get standard Laguerre nodes and weights with parameter a - 1
  quad <- gauss_laguerre_nodes(M, alpha_param = a - 1)

  # Transform nodes to alpha scale
  alpha_nodes <- quad$nodes / b
  if (any(!is.finite(alpha_nodes))) {
    .dpprior_abort_invalid(
      "Gamma quadrature produced non-finite transformed nodes; use less extreme a/b values",
      "dpprior_numerical_domain_error", "a/b", c(a = a, b = b),
      "finite transformed quadrature nodes", "nonfinite_transform"
    )
  }

  # Normalize weights in log-space for numerical stability. This does not
  # imply monotone convergence in M for a general integrand.
  logw <- quad$weights_log
  logw_norm <- logw - logsumexp_vec(logw)
  weights_normalized <- exp(logw_norm)
  normalized_sum <- sum(weights_normalized)
  if (any(!is.finite(weights_normalized)) || !is.finite(normalized_sum) ||
      normalized_sum <= 0) {
    .dpprior_abort_invalid(
      "Gamma quadrature produced invalid normalized weights; use less extreme a values",
      "dpprior_numerical_domain_error", "a", a,
      "finite normalized quadrature weights", "nonfinite_weights"
    )
  }

  result <- list(
    a = a,
    b = b,
    alpha_nodes = alpha_nodes,
    weights_normalized = weights_normalized
  )
  node_metadata <- .quadrature_metadata(quad)
  attr(result, "quadrature_metadata") <- list(
    schema_version = 1L,
    distribution = "Gamma(shape, rate)",
    shape = a,
    rate = b,
    M_selected = M,
    M_verification = NA_integer_,
    normalization = "log-sum-exp",
    normalized_weight_sum = normalized_sum,
    node_rule = node_metadata
  )
  result
}


#' Integrate Function Against Gamma Distribution
#'
#' High-level interface for computing \eqn{E[f(\alpha)]} where
#' \eqn{\alpha \sim \text{Gamma}(a, b)}.
#'
#' @param f Function to integrate; must accept a single numeric argument
#'   and return a single numeric value.
#' @param a Numeric; shape parameter of Gamma distribution (must be > 0).
#' @param b Numeric; rate parameter of Gamma distribution (must be > 0).
#' @param M Integer; number of quadrature nodes (default: 80).
#'
#' @return Numeric; approximation of \eqn{E[f(\alpha)]}.
#'
#' @details
#' Uses Gauss-Laguerre quadrature to approximate:
#' \deqn{E[f(\alpha)] = \int_0^\infty f(\alpha) \frac{b^a}{\Gamma(a)} \alpha^{a-1} e^{-b\alpha} d\alpha}
#'
#' The approximation is:
#' \deqn{E[f(\alpha)] \approx \sum_{m=1}^M w_m f(\alpha_m)}
#'
#' where \eqn{\alpha_m} are quadrature nodes and \eqn{w_m} are normalized weights.
#'
#' For polynomial integrands of degree up to \eqn{2M - 1}, the quadrature is exact.
#' For other smooth functions, accuracy improves rapidly with \eqn{M}.
#'
#' @examples
#' # E[alpha] for Gamma(2.5, 1.5) should be 2.5/1.5
#' integrate_gamma(identity, 2.5, 1.5)
#'
#' # E[alpha^2] for Gamma(2.5, 1.5) should be 2.5*3.5/1.5^2
#' integrate_gamma(function(x) x^2, 2.5, 1.5)
#'
#' # More complex function
#' integrate_gamma(function(x) log(x + 1), 2.5, 1.5)
#'
#' @seealso \code{\link{build_gamma_quadrature}} for the underlying quadrature
#'
#' @export
integrate_gamma <- function(f, a, b, M = .QUAD_NODES_DEFAULT) {
  # Build quadrature
  quad <- build_gamma_quadrature(a, b, M)

  # Evaluate function at quadrature nodes
  f_vals <- .quadrature_evaluate_integrand(f, quad$alpha_nodes)

  # Stable weighted sum
  .quadrature_weighted_sum(quad$weights_normalized, f_vals)
}


#' Verify Quadrature Accuracy Against Known Gamma Moments
#'
#' Validates quadrature implementation by comparing computed expectations
#' against known closed-form Gamma distribution moments.
#'
#' @param a Numeric; shape parameter of Gamma distribution.
#' @param b Numeric; rate parameter of Gamma distribution.
#' @param M Integer; number of quadrature nodes.
#' @param tol Numeric; tolerance for verification (default: 1e-10).
#' @param verbose Logical; if \code{TRUE}, print detailed results.
#'
#' @return Logical; \code{TRUE} if all moments match within tolerance.
#'
#' @details
#' For \eqn{\alpha \sim \text{Gamma}(a, b)}:
#' \itemize{
#'   \item \eqn{E[\alpha] = a/b}
#'   \item \eqn{E[\alpha^2] = a(a+1)/b^2}
#'   \item \eqn{Var(\alpha) = a/b^2}
#' }
#'
#' @examples
#' \dontrun{
#' # Should return TRUE
#' verify_quadrature(2.5, 1.5, M = 80)
#'
#' # More challenging case
#' verify_quadrature(0.5, 2.0, M = 100, verbose = TRUE)
#'
#' }
#' @keywords internal
verify_quadrature <- function(a, b, M, tol = 1e-10, verbose = TRUE) {
  assert_positive(a, "a")
  assert_positive(b, "b")

  # Theoretical moments
  E_alpha_true <- a / b
  E_alpha_sq_true <- a * (a + 1) / b^2
  Var_alpha_true <- a / b^2

  # Computed via quadrature
  E_alpha_quad <- integrate_gamma(identity, a, b, M)
  E_alpha_sq_quad <- integrate_gamma(function(x) x^2, a, b, M)
  Var_alpha_quad <- E_alpha_sq_quad - E_alpha_quad^2

  # Errors
  err_mean <- abs(E_alpha_quad - E_alpha_true)
  err_sq <- abs(E_alpha_sq_quad - E_alpha_sq_true)
  err_var <- abs(Var_alpha_quad - Var_alpha_true)

  pass_mean <- err_mean < tol
  pass_sq <- err_sq < tol
  pass_var <- err_var < tol
  all_pass <- pass_mean && pass_sq && pass_var

  if (verbose) {
    cat(sprintf("Quadrature verification (a=%.2f, b=%.2f, M=%d):\n", a, b, M))
    cat(sprintf("  E[alpha]:   true=%.6f, quad=%.10f, error=%.2e [%s]\n",
                E_alpha_true, E_alpha_quad, err_mean,
                if (pass_mean) "PASS" else "FAIL"))
    cat(sprintf("  E[alpha^2]:  true=%.6f, quad=%.10f, error=%.2e [%s]\n",
                E_alpha_sq_true, E_alpha_sq_quad, err_sq,
                if (pass_sq) "PASS" else "FAIL"))
    cat(sprintf("  Var(alpha): true=%.6f, quad=%.10f, error=%.2e [%s]\n",
                Var_alpha_true, Var_alpha_quad, err_var,
                if (pass_var) "PASS" else "FAIL"))
    cat(sprintf("  Overall: %s\n", if (all_pass) "PASS" else "FAIL"))
  }

  invisible(all_pass)
}


#' Get Quadrature Information Summary
#'
#' Returns summary information about the quadrature nodes and weights
#' for a given Gamma(a, b) distribution.
#'
#' @param a Numeric; shape parameter of Gamma distribution.
#' @param b Numeric; rate parameter of Gamma distribution.
#' @param M Integer; number of quadrature nodes.
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{n_nodes}}{Number of quadrature nodes.}
#'     \item{\code{alpha_range}}{Range of alpha nodes (min, max).}
#'     \item{\code{weight_range}}{Range of normalized weights (min, max).}
#'     \item{\code{gamma_mean}}{Theoretical mean of Gamma(a, b).}
#'     \item{\code{gamma_sd}}{Theoretical SD of Gamma(a, b).}
#'     \item{\code{coverage}}{Approximate coverage in terms of SD from mean.}
#'   }
#'
#' @examples
#' \dontrun{
#' summary_quadrature(2.5, 1.5, M = 80)
#'
#' }
#' @keywords internal
summary_quadrature <- function(a, b, M = .QUAD_NODES_DEFAULT) {
  quad <- build_gamma_quadrature(a, b, M)

  gamma_mean <- a / b
  gamma_sd <- sqrt(a) / b

  list(
    n_nodes = M,
    alpha_range = c(min = min(quad$alpha_nodes), max = max(quad$alpha_nodes)),
    weight_range = c(min = min(quad$weights_normalized),
                     max = max(quad$weights_normalized)),
    gamma_mean = gamma_mean,
    gamma_sd = gamma_sd,
    coverage = c(
      lower_sd = (gamma_mean - min(quad$alpha_nodes)) / gamma_sd,
      upper_sd = (max(quad$alpha_nodes) - gamma_mean) / gamma_sd
    )
  )
}


#' Convergence Diagnostic for Quadrature
#'
#' Examines how quadrature accuracy improves with increasing number of nodes.
#'
#' @param f Function to integrate.
#' @param a Numeric; shape parameter of Gamma distribution.
#' @param b Numeric; rate parameter of Gamma distribution.
#' @param M_values Integer vector; number of nodes to test.
#' @param true_value Numeric; known true value (optional).
#'
#' @return A data frame with columns: M, estimate, change, relative_change.
#'
#' @examples
#' \dontrun{
#' # Check convergence for E[alpha]
#' convergence_quadrature(identity, 2.5, 1.5,
#'                        M_values = c(10, 20, 50, 80, 100),
#'                        true_value = 2.5/1.5)
#'
#' }
#' @keywords internal
convergence_quadrature <- function(f, a, b, M_values = c(10, 20, 50, 80, 100),
                                   true_value = NULL) {
  estimates <- sapply(M_values, function(M) integrate_gamma(f, a, b, M))

  result <- data.frame(
    M = M_values,
    estimate = estimates,
    change = c(NA, diff(estimates)),
    relative_change = c(NA, abs(diff(estimates)) / abs(estimates[-length(estimates)]))
  )

  if (!is.null(true_value)) {
    result$error <- abs(estimates - true_value)
    result$relative_error <- result$error / abs(true_value)
  }

  result
}


#' Independently Verify a Gamma Quadrature Integral at Higher Order
#'
#' Internal helper that evaluates the same contracted integrand at selected
#' and higher node counts, records both computation metadata objects, and
#' applies a combined absolute/relative tolerance to their difference.
#'
#' @param f Function to integrate.
#' @param a,b Gamma shape and rate.
#' @param M Selected quadrature order.
#' @param M_verify Independent verification order. The minimum admissible order
#'   is \code{max(2*M, M+40)}. When \code{NULL}, that minimum is used. Because
#'   the supported ceiling is 512 nodes, verification is unavailable when
#'   \code{M > 256}; the function then raises a typed verification error even
#'   if an explicit \code{M_verify} is supplied.
#' @param abs_tol,rel_tol Non-negative finite tolerances.
#'
#' @return A list containing selected and verification estimates, absolute
#'   difference, tolerance, pass flag, and versioned metadata.
#'
#' @keywords internal
.verify_gamma_quadrature_order <- function(
    f, a, b, M = .QUAD_NODES_DEFAULT,
    M_verify = NULL,
    abs_tol = 1e-12, rel_tol = 1e-10) {
  if (!is.function(f)) {
    .dpprior_abort_invalid(
      "f must be a function",
      "dpprior_integrand_error", "f", f, "function", "type"
    )
  }
  M <- .dpprior_validate_count(
    M, "M", minimum = 1L, maximum = .QUADRATURE_MAX_NODES,
    .subclass = "dpprior_control_error"
  )
  M_verification_required <-
    .quadrature_verification_required_order(M)
  verification_available <-
    M_verification_required <= .QUADRATURE_MAX_NODES
  if (!verification_available) {
    .dpprior_abort_invalid(
      sprintf(
        paste(
          "quadrature verification is unavailable for M=%d:",
          "required M_verify=%d exceeds the supported ceiling (%d)"
        ),
        M, M_verification_required, .QUADRATURE_MAX_NODES
      ),
      c("dpprior_quadrature_verification_error", "dpprior_bounds_error"),
      "M_verify", M_verification_required,
      sprintf("required order <= %d", .QUADRATURE_MAX_NODES),
      "verification_unavailable"
    )
  }
  if (is.null(M_verify)) {
    M_verify <- M_verification_required
  }
  M_verify <- .dpprior_validate_count(
    M_verify, "M_verify", minimum = 1L,
    maximum = .QUADRATURE_MAX_NODES,
    .subclass = "dpprior_quadrature_verification_error"
  )
  if (M_verify < M_verification_required) {
    .dpprior_abort_invalid(
      sprintf(
        "M_verify must be at least %d for selected order M=%d",
        M_verification_required, M
      ),
      c("dpprior_quadrature_verification_error", "dpprior_bounds_error"),
      "M_verify", M_verify,
      sprintf(
        "integer in [%d, %d]",
        M_verification_required, .QUADRATURE_MAX_NODES
      ),
      "insufficient_verification_order"
    )
  }
  abs_tol <- .dpprior_validate_scalar(
    abs_tol, "abs_tol", lower = 0,
    .subclass = "dpprior_control_error"
  )
  rel_tol <- .dpprior_validate_scalar(
    rel_tol, "rel_tol", lower = 0,
    .subclass = "dpprior_control_error"
  )

  selected_quad <- build_gamma_quadrature(a, b, M)
  verification_quad <- build_gamma_quadrature(a, b, M_verify)
  selected_values <- .quadrature_evaluate_integrand(
    f, selected_quad$alpha_nodes
  )
  verification_values <- .quadrature_evaluate_integrand(
    f, verification_quad$alpha_nodes
  )
  selected <- .quadrature_weighted_sum(
    selected_quad$weights_normalized, selected_values,
    fail_on_nonfinite = FALSE
  )
  verification <- .quadrature_weighted_sum(
    verification_quad$weights_normalized, verification_values,
    fail_on_nonfinite = FALSE
  )
  abs_difference <- abs(selected - verification)
  tolerance <- abs_tol + rel_tol * max(abs(selected), abs(verification))
  estimates_finite <- is.finite(selected) && is.finite(verification)
  passed <- estimates_finite &&
    is.finite(abs_difference) && abs_difference <= tolerance

  list(
    status = if (!estimates_finite) {
      "failed"
    } else if (passed) {
      "converged"
    } else {
      "approximate"
    },
    passed = passed,
    selected = selected,
    verification = verification,
    abs_difference = abs_difference,
    tolerance = tolerance,
    metadata = list(
      schema_version = 1L,
      M_selected = as.integer(M),
      M_verification = as.integer(M_verify),
      M_verification_required = M_verification_required,
      verification_available = verification_available,
      selected = .quadrature_metadata(selected_quad),
      verification = .quadrature_metadata(verification_quad)
    )
  )
}

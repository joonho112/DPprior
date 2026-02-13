# =============================================================================
# Module 01: Log-Stirling Numbers of the First Kind
# =============================================================================
#
# This module provides numerically stable computation of unsigned Stirling
# numbers of the first kind |s(J,k)| in log-space.
#
# Theory Background:
# -----------------
# The unsigned Stirling numbers of the first kind |s(J,k)| count the number
# of permutations of J elements with exactly k disjoint cycles.
#
# They satisfy the recurrence:
#   |s(J,k)| = |s(J-1,k-1)| + (J-1) * |s(J-1,k)|
#
# with boundary conditions:
#   |s(0,0)| = 1
#   |s(J,0)| = 0  for J >= 1
#   |s(J,J)| = 1  for all J
#
# For the Antoniak distribution (exact PMF of K_J | alpha in a DP):
#   P(K_J = k | alpha) = |s(J,k)| * alpha^k / (alpha)_J
#
# where (alpha)_J = alpha * (alpha+1) * ... * (alpha+J-1) is the rising factorial.
#
# Author: JoonHo Lee
# Date: December 2025
# Part of: DPprior R Package
# Reference: Lee (2026), Sections 2--3
# =============================================================================


#' Compute Log Stirling Numbers (First Kind, Unsigned)
#'
#' Computes the logarithm of unsigned Stirling numbers of the first kind
#' \eqn{|s(J,k)|} for all \eqn{J} from 0 to \code{J_max} and \eqn{k} from 0 to
#' \eqn{J}. Uses log-space recursion to prevent numerical overflow.
#'
#' @param J_max Maximum J value (non-negative integer, at most 500).
#'
#' @return A lower triangular matrix of dimension \code{(J_max+1) x (J_max+1)}.
#'   Entry \code{[J+1, k+1]} contains \eqn{\log|s(J,k)|} (using R's 1-based indexing).
#'   Invalid entries (k > J or k < 1 when J >= 1) contain \code{-Inf}.
#'
#' @details
#' The unsigned Stirling numbers of the first kind \eqn{|s(J,k)|} count the
#' number of permutations of \eqn{J} elements into exactly \eqn{k} disjoint
#' cycles.
#'
#' Uses the log-space recurrence relation:
#' \deqn{L_{J,k} = \text{logsumexp}(L_{J-1,k-1}, \log(J-1) + L_{J-1,k})}
#'
#' where \eqn{L_{J,k} = \log|s(J,k)|}.
#'
#' Boundary conditions:
#' \itemize{
#'   \item \eqn{|s(0,0)| = 1} \eqn{\rightarrow} \code{L[1,1] = 0}
#'   \item \eqn{|s(J,0)| = 0} for \eqn{J \geq 1} \eqn{\rightarrow} \code{L[J+1,1] = -Inf}
#'   \item \eqn{|s(J,J)| = 1} for all \eqn{J} \eqn{\rightarrow} \code{L[J+1,J+1] = 0}
#' }
#'
#' Complexity: O(\code{J_max}^2) time and space. Results should be cached for
#' repeated use within a session.
#'
#' @examples
#' # Compute Stirling numbers up to J=10
#' logS <- compute_log_stirling(10)
#'
#' # Access |s(4,2)| = 11
#' exp(logS[5, 3])
#'
#' # Access |s(5,3)| = 35
#' exp(logS[6, 4])
#'
#' @references
#' Antoniak, C. E. (1974). Mixtures of Dirichlet Processes with Applications
#' to Bayesian Nonparametric Problems. \emph{The Annals of Statistics},
#' 2(6), 1152-1174.
#'
#' @seealso \code{\link{get_log_stirling}} for safe accessor with bounds checking
#'
#' @export
compute_log_stirling <- function(J_max) {
  # Input validation
  if (!is.numeric(J_max) || length(J_max) != 1L || !is.finite(J_max) ||
      J_max != floor(J_max) || J_max < 0L) {
    stop("J_max must be a non-negative integer", call. = FALSE)
  }
  if (J_max > .MAX_J_DEFAULT) {
    stop(sprintf("J_max exceeds maximum supported value (%d)", .MAX_J_DEFAULT),
         call. = FALSE)
  }
  
  J_max <- as.integer(J_max)
  
  # Initialize matrix with -Inf (log(0))
  logS <- matrix(-Inf, nrow = J_max + 1L, ncol = J_max + 1L)
  
  # Boundary: s(0,0) = 1 -> log = 0
  logS[1L, 1L] <- 0
  
  # Handle J_max = 0 case
  if (J_max == 0L) {
    return(logS)
  }
  
  # Fill recursively using log-space operations
  for (J in seq_len(J_max)) {
    # s(J, J) = 1 -> log = 0
    logS[J + 1L, J + 1L] <- 0
    
    # s(J, k) for 1 <= k < J using recurrence:
    # |s(J,k)| = |s(J-1,k-1)| + (J-1)|s(J-1,k)|
    # In log-space: logsumexp(L[J-1,k-1], log(J-1) + L[J-1,k])
    if (J >= 2L) {
      log_Jm1 <- log(J - 1L)
      for (k in seq_len(J - 1L)) {
        # term1: log|s(J-1, k-1)| = logS[J, k]
        # term2: log((J-1) * |s(J-1, k)|) = log(J-1) + logS[J, k+1]
        term1 <- logS[J, k]
        term2 <- log_Jm1 + logS[J, k + 1L]
        logS[J + 1L, k + 1L] <- logsumexp(term1, term2)
      }
    }
    
    # s(J, 0) = 0 (already -Inf)
  }
  
  logS
}


#' Get Single Log-Stirling Value with Bounds Checking
#'
#' Safe accessor for the log-Stirling matrix with automatic bounds checking.
#' Returns \code{-Inf} for invalid indices (k > J, k < 1).
#'
#' @param J Sample size (non-negative integer).
#' @param k Number of clusters (integer).
#' @param logS Pre-computed log-Stirling matrix from \code{\link{compute_log_stirling}}.
#'
#' @return The value \eqn{\log|s(J,k)|}, or \code{-Inf} if k > J or k < 1.
#'
#' @examples
#' \dontrun{
#' logS <- compute_log_stirling(10)
#'
#' # Valid access
#' get_log_stirling(4, 2, logS)
#'
#' # Invalid access returns -Inf
#' get_log_stirling(4, 5, logS)
#' get_log_stirling(4, 0, logS)
#'
#' }
#' @seealso \code{\link{compute_log_stirling}} for matrix computation
#'
#' @keywords internal
get_log_stirling <- function(J, k, logS) {
  # Validate logS
  if (!is.matrix(logS) || nrow(logS) != ncol(logS)) {
    stop("logS must be a square matrix returned by compute_log_stirling()",
         call. = FALSE)
  }
  
  # Validate J
  if (!is.numeric(J) || length(J) != 1L || !is.finite(J) ||
      J != floor(J) || J < 0L) {
    stop("J must be a non-negative integer", call. = FALSE)
  }
  
  # Validate k
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k != floor(k)) {
    stop("k must be an integer", call. = FALSE)
  }
  
  # Per package convention, k must be in {1,...,J} for valid Stirling numbers
  if (k < 1L || k > J) {
    return(-Inf)
  }
  
  # Check matrix bounds
  if (J + 1L > nrow(logS)) {
    stop(sprintf("J=%d exceeds precomputed matrix size (J_max=%d)", 
                 J, nrow(logS) - 1L), call. = FALSE)
  }
  
  # Return the value (R is 1-indexed)
  logS[J + 1L, k + 1L]
}


#' Validate Stirling Number Computation
#'
#' Validates a log-Stirling matrix against known reference values.
#' Useful for testing and verification.
#'
#' @param logS Pre-computed log-Stirling matrix from \code{\link{compute_log_stirling}}.
#' @param verbose Logical; if \code{TRUE}, print validation results.
#'
#' @return Logical indicating whether all validations passed.
#'
#' @details
#' Checks against known values:
#' \itemize{
#'   \item \eqn{|s(4,2)| = 11}
#'   \item \eqn{|s(5,3)| = 35}
#'   \item \eqn{|s(10,5)| = 269325}
#' }
#'
#' @examples
#' logS <- compute_log_stirling(10)
#' validate_stirling(logS)
#'
#' @export
validate_stirling <- function(logS, verbose = TRUE) {
  if (!is.matrix(logS) || nrow(logS) != ncol(logS)) {
    if (verbose) message("FAIL: logS is not a valid square matrix")
    return(FALSE)
  }
  
  # Known reference values
  known_values <- data.frame(
    J = c(4L, 5L, 6L, 10L),
    k = c(2L, 3L, 3L, 5L),
    expected = c(11, 35, 225, 269325)
  )
  
  J_max <- nrow(logS) - 1L
  all_passed <- TRUE
  
  for (i in seq_len(nrow(known_values))) {
    J <- known_values$J[i]
    k <- known_values$k[i]
    expected <- known_values$expected[i]
    
    # Skip if matrix is too small
    if (J > J_max) {
      if (verbose) {
        message(sprintf("SKIP: |s(%d,%d)| - matrix too small", J, k))
      }
      next
    }
    
    computed <- round(exp(logS[J + 1L, k + 1L]))
    passed <- (abs(computed - expected) < 0.5)
    all_passed <- all_passed && passed
    
    if (verbose) {
      status <- if (passed) "PASS" else "FAIL"
      message(sprintf("%s: |s(%d,%d)| = %d (expected %d)", 
                      status, J, k, computed, expected))
    }
  }
  
  all_passed
}


#' Verify Row Sum Identity for Stirling Numbers
#'
#' Verifies that the row sum of Stirling numbers equals J! for each row.
#' This is the fundamental identity: \eqn{\sum_{k=1}^{J} |s(J,k)| = J!}
#'
#' @param logS Pre-computed log-Stirling matrix from \code{\link{compute_log_stirling}}.
#' @param J_values Vector of J values to verify (default: 2:10).
#' @param tolerance Numerical tolerance for comparison (default: 1e-10).
#' @param verbose Logical; if \code{TRUE}, print verification results.
#'
#' @return Logical indicating whether all verifications passed.
#'
#' @details
#' The unsigned Stirling numbers of the first kind satisfy:
#' \deqn{\sum_{k=1}^{J} |s(J,k)| = J!}
#'
#' This identity follows from the fact that \eqn{|s(J,k)|} counts permutations
#' of J elements with exactly k cycles, and the total number of permutations
#' is J!.
#'
#' @examples
#' \dontrun{
#' logS <- compute_log_stirling(15)
#' verify_stirling_row_sum(logS, J_values = 2:10)
#'
#' }
#' @keywords internal
verify_stirling_row_sum <- function(logS, J_values = 2:10,
                                     tolerance = 1e-10, verbose = TRUE) {
  if (!is.matrix(logS) || nrow(logS) != ncol(logS)) {
    stop("logS must be a square matrix returned by compute_log_stirling()",
         call. = FALSE)
  }
  
  all_passed <- TRUE
  J_max <- nrow(logS) - 1L
  
  for (J in J_values) {
    if (J > J_max) {
      if (verbose) {
        message(sprintf("SKIP: J=%d exceeds matrix dimensions", J))
      }
      next
    }
    
    # log(sum_{k=1}^{J} |s(J,k)|)
    log_row_sum <- logsumexp_vec(logS[J + 1L, 2L:(J + 1L)])
    
    # log(J!)
    log_factorial <- sum(log(seq_len(J)))
    
    error <- abs(log_row_sum - log_factorial)
    passed <- error < tolerance
    all_passed <- all_passed && passed
    
    if (verbose) {
      status <- if (passed) "PASS" else "FAIL"
      message(sprintf("%s: sum_k |s(%d,k)| = %d! [error = %.2e]",
                      status, J, J, error))
    }
  }
  
  all_passed
}


#' Get Stirling Numbers for a Fixed J
#'
#' Returns a vector of \eqn{\log|s(J,k)|} for \eqn{k = 1, \ldots, J}.
#'
#' @param J Sample size (positive integer).
#' @param logS Pre-computed log-Stirling matrix from \code{\link{compute_log_stirling}}.
#'
#' @return Numeric vector of length J containing \eqn{\log|s(J,k)|} for \eqn{k = 1, \ldots, J}.
#'
#' @examples
#' \dontrun{
#' logS <- compute_log_stirling(10)
#' log_s_row <- get_stirling_row(5, logS)
#' exp(log_s_row)
#'
#' }
#' @keywords internal
get_stirling_row <- function(J, logS) {
  if (!is.matrix(logS) || nrow(logS) != ncol(logS)) {
    stop("logS must be a square matrix returned by compute_log_stirling()",
         call. = FALSE)
  }
  if (!is.numeric(J) || length(J) != 1L || !is.finite(J) ||
      J != floor(J) || J < 1L) {
    stop("J must be a positive integer", call. = FALSE)
  }
  if (J + 1L > nrow(logS)) {
    stop(sprintf("J=%d exceeds precomputed matrix size (J_max=%d)",
                 J, nrow(logS) - 1L), call. = FALSE)
  }
  
  logS[J + 1L, 2L:(J + 1L)]
}

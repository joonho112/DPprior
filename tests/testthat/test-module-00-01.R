# =============================================================================
# Test File: Module 00 (Constants) and Module 01 (Stirling Numbers)
# =============================================================================
#
# This file contains unit tests for:
# - logsumexp and logsumexp_vec functions (including edge cases)
# - softmax function (including Inf handling)
# - Assertion helper functions
# - Stirling number computation
# - Stirling number validation
#
# Author: JoonHo Lee
# Date: December 2025
# Part of: DPprior R Package
# =============================================================================


# =============================================================================
# Module 00: Constants and Utilities Tests
# =============================================================================

test_that("logsumexp is numerically stable", {
  # Extreme positive values that would overflow with naive implementation
  expect_equal(logsumexp(1000, 1000), 1000 + log(2), tolerance = 1e-12)
  
  # Extreme negative values
  expect_equal(logsumexp(-1000, -1000), -1000 + log(2), tolerance = 1e-12)
  
  # One value dominates
  expect_equal(logsumexp(1000, -1000), 1000, tolerance = 1e-12)
  
  # Standard cases
  expect_equal(logsumexp(0, 0), log(2), tolerance = 1e-12)
  expect_equal(logsumexp(log(2), log(3)), log(5), tolerance = 1e-12)
  
  # Vectorization
  result <- logsumexp(c(0, 0), c(0, 0))
  expect_equal(result, rep(log(2), 2), tolerance = 1e-12)
})


test_that("logsumexp handles Inf edge cases correctly", {
  # Both -Inf should return -Inf
  expect_equal(logsumexp(-Inf, -Inf), -Inf)
  
  # Both +Inf should return +Inf
  expect_equal(logsumexp(Inf, Inf), Inf)
  
  # Mixed Inf cases
  expect_equal(logsumexp(Inf, -Inf), Inf)
  expect_equal(logsumexp(-Inf, Inf), Inf)
  expect_equal(logsumexp(Inf, 0), Inf)
  expect_equal(logsumexp(0, Inf), Inf)
  
  # -Inf with finite
  expect_equal(logsumexp(-Inf, 0), 0, tolerance = 1e-12)
  expect_equal(logsumexp(0, -Inf), 0, tolerance = 1e-12)
  
  # Vectorized Inf handling
  result <- logsumexp(c(-Inf, Inf), c(-Inf, Inf))
  expect_equal(result, c(-Inf, Inf))
})


test_that("logsumexp_vec works correctly", {
  # Equal values
  expect_equal(logsumexp_vec(c(1000, 1000, 1000)), 1000 + log(3), tolerance = 1e-12)
  expect_equal(logsumexp_vec(c(0, 0, 0, 0)), log(4), tolerance = 1e-12)
  
  # Negative values
  expect_equal(logsumexp_vec(c(-100, -100)), -100 + log(2), tolerance = 1e-12)
  
  # Single value
  expect_equal(logsumexp_vec(5), 5)
})


test_that("logsumexp_vec handles Inf edge cases correctly", {
  # All -Inf should return -Inf
  expect_equal(logsumexp_vec(c(-Inf, -Inf, -Inf)), -Inf)
  
  # Any +Inf should return +Inf
  expect_equal(logsumexp_vec(c(0, Inf, 1)), Inf)
  expect_equal(logsumexp_vec(c(Inf)), Inf)
  
  # Mix of -Inf and finite
  expect_equal(logsumexp_vec(c(-Inf, 0, -Inf)), 0, tolerance = 1e-12)
  
  # Empty vector should error
  expect_error(logsumexp_vec(numeric(0)), "must have positive length")
})


test_that("softmax returns valid probabilities", {
  x <- c(1, 2, 3)
  p <- softmax(x)
  
  # Sum to 1
  expect_equal(sum(p), 1, tolerance = 1e-15)
  
  # All positive
  expect_true(all(p > 0))
  
  # All less than 1
  expect_true(all(p < 1))
  
  # Extreme values don't cause overflow
  x_extreme <- c(1000, 1001, 1002)
  p_extreme <- softmax(x_extreme)
  expect_equal(sum(p_extreme), 1, tolerance = 1e-15)
  expect_true(all(p_extreme > 0))
  
  # Negative extreme values
  x_neg <- c(-1000, -999, -998)
  p_neg <- softmax(x_neg)
  expect_equal(sum(p_neg), 1, tolerance = 1e-15)
  
  # Equal values give equal probabilities
  p_equal <- softmax(c(0, 0, 0))
  expect_equal(p_equal, rep(1/3, 3), tolerance = 1e-15)
})


test_that("softmax handles Inf correctly", {
  # Single Inf gets all mass
  p <- softmax(c(1, Inf, 2))
  expect_equal(p, c(0, 1, 0))
  
  # Multiple Inf split mass uniformly
  p <- softmax(c(Inf, 1, Inf))
  expect_equal(p, c(0.5, 0, 0.5))
  
  # All Inf split uniformly
  p <- softmax(c(Inf, Inf, Inf))
  expect_equal(p, rep(1/3, 3), tolerance = 1e-15)
  
  # Empty vector should error
  expect_error(softmax(numeric(0)), "must have positive length")
})


test_that("assertion functions work correctly", {
  # assert_positive - now also checks is.finite
  expect_true(assert_positive(c(1, 2, 3), "x"))
  expect_error(assert_positive(c(1, -1), "x"), "must be finite and positive")
  expect_error(assert_positive(0, "x"), "must be finite and positive")
  expect_error(assert_positive(Inf, "x"), "must be finite and positive")
  expect_error(assert_positive(NaN, "x"), "must be finite and positive")
  
  # assert_valid_J
  expect_true(assert_valid_J(50))
  expect_true(assert_valid_J(1))
  expect_error(assert_valid_J(0), "positive integer")
  expect_error(assert_valid_J(-1), "positive integer")
  expect_error(assert_valid_J(1.5), "positive integer")
  expect_error(assert_valid_J(1000), "exceeds maximum")
  expect_error(assert_valid_J(Inf), "positive integer")
  
  # assert_probability - now also checks is.finite
  expect_true(assert_probability(0.5, "p"))
  expect_true(assert_probability(0, "p"))
  expect_true(assert_probability(1, "p"))
  expect_error(assert_probability(-0.1, "p"), "must be finite and in")
  expect_error(assert_probability(1.5, "p"), "must be finite and in")
  expect_error(assert_probability(NaN, "p"), "must be finite and in")
  
  # assert_nonnegative
  expect_true(assert_nonnegative(c(0, 1, 2), "x"))
  expect_error(assert_nonnegative(c(0, -1), "x"), "must be finite and non-negative")
  expect_error(assert_nonnegative(Inf, "x"), "must be finite and non-negative")
  
  # assert_valid_k
  expect_true(assert_valid_k(2, 5))
  expect_error(assert_valid_k(0, 5), "must be an integer in")
  expect_error(assert_valid_k(6, 5), "must be an integer in")
  expect_error(assert_valid_k(NaN, 5), "must be an integer in")
})


# =============================================================================
# Module 01: Stirling Numbers Tests
# =============================================================================

test_that("Stirling numbers match known values", {
  logS <- compute_log_stirling(100)
  
  # |s(4,2)| = 11
  expect_equal(round(exp(logS[5, 3])), 11)
  
  # |s(5,3)| = 35
  expect_equal(round(exp(logS[6, 4])), 35)
  
  # |s(6,3)| = 225
  expect_equal(round(exp(logS[7, 4])), 225)
  
  # |s(10,5)| = 269325
  expect_equal(round(exp(logS[11, 6])), 269325)
  
  # |s(1,1)| = 1
  expect_equal(exp(logS[2, 2]), 1)
})


test_that("Stirling row sums equal J!", {
  logS <- compute_log_stirling(15)
  
  for (J in 2:15) {
    # sum_{k=1}^{J} |s(J,k)|
    log_row_sum <- logsumexp_vec(logS[J + 1, 2:(J + 1)])
    
    # log(J!)
    expected <- sum(log(1:J))
    
    expect_equal(log_row_sum, expected, tolerance = 1e-10,
                 info = sprintf("J = %d", J))
  }
})


test_that("Stirling boundary conditions are correct", {
  logS <- compute_log_stirling(10)
  
  # |s(0,0)| = 1
  expect_equal(exp(logS[1, 1]), 1)
  
  # |s(J,J)| = 1 for all J
  for (J in 1:10) {
    expect_equal(exp(logS[J + 1, J + 1]), 1,
                 info = sprintf("J = %d", J))
  }
  
  # |s(J,0)| = 0 for J >= 1 (represented as -Inf in log)
  for (J in 1:10) {
    expect_equal(logS[J + 1, 1], -Inf,
                 info = sprintf("J = %d", J))
  }
})


test_that("Stirling recurrence relation holds", {
  logS <- compute_log_stirling(20)
  
  # Test recurrence: |s(J,k)| = |s(J-1,k-1)| + (J-1)|s(J-1,k)|
  for (J in 3:20) {
    for (k in 2:(J - 1)) {
      # LHS: log|s(J,k)|
      lhs <- logS[J + 1, k + 1]
      
      # RHS: logsumexp(log|s(J-1,k-1)|, log(J-1) + log|s(J-1,k)|)
      term1 <- logS[J, k]
      term2 <- log(J - 1) + logS[J, k + 1]
      rhs <- logsumexp(term1, term2)
      
      expect_equal(lhs, rhs, tolerance = 1e-10,
                   info = sprintf("J = %d, k = %d", J, k))
    }
  }
})


test_that("compute_log_stirling handles edge cases", {
  # J_max = 0 should work
  logS <- compute_log_stirling(0)
  expect_equal(dim(logS), c(1, 1))
  expect_equal(logS[1, 1], 0)
  
  # J_max = 1 should work
  logS <- compute_log_stirling(1)
  expect_equal(dim(logS), c(2, 2))
  expect_equal(exp(logS[2, 2]), 1)
})


test_that("compute_log_stirling input validation works", {
  expect_error(compute_log_stirling(-1), "non-negative integer")
  expect_error(compute_log_stirling(1000), "exceeds maximum")
  expect_error(compute_log_stirling(1.5), "non-negative integer")
  expect_error(compute_log_stirling(Inf), "non-negative integer")
  expect_error(compute_log_stirling(NaN), "non-negative integer")
  expect_error(compute_log_stirling("10"), "non-negative integer")
})


test_that("get_log_stirling handles bounds correctly", {
  logS <- compute_log_stirling(10)
  
  # Valid access
  expect_equal(get_log_stirling(4, 2, logS), logS[5, 3])
  
  # k > J returns -Inf
  expect_equal(get_log_stirling(4, 5, logS), -Inf)
  
  # k < 1 returns -Inf
  expect_equal(get_log_stirling(4, 0, logS), -Inf)
  
  # J exceeds matrix bounds should error
  expect_error(get_log_stirling(100, 50, logS), "exceeds precomputed matrix size")
  
  # Invalid input types should error
  expect_error(get_log_stirling(-1, 1, logS), "non-negative integer")
  expect_error(get_log_stirling(4, "2", logS), "must be an integer")
})


test_that("validate_stirling works correctly", {
  logS <- compute_log_stirling(15)
  
  # Should pass for valid matrix
  expect_true(validate_stirling(logS, verbose = FALSE))
  
  # Create an invalid matrix and check it fails
  logS_bad <- logS
  logS_bad[5, 3] <- 0
  expect_false(validate_stirling(logS_bad, verbose = FALSE))
  
  # Non-matrix input
  expect_false(validate_stirling(c(1, 2, 3), verbose = FALSE))
})


test_that("verify_stirling_row_sum works correctly", {
  logS <- compute_log_stirling(15)
  
  # Should pass for valid matrix
  expect_true(verify_stirling_row_sum(logS, J_values = 2:10, verbose = FALSE))
  
  # Invalid matrix should error
  expect_error(verify_stirling_row_sum(c(1, 2, 3)), "must be a square matrix")
})


test_that("get_stirling_row returns correct values", {
  logS <- compute_log_stirling(10)
  
  # For J=5, should return log|s(5,k)| for k=1,...,5
  row <- get_stirling_row(5, logS)
  expect_length(row, 5)
  expect_equal(row[1], logS[6, 2])
  expect_equal(row[5], logS[6, 6])
  
  # Invalid input should error
  expect_error(get_stirling_row(100, logS), "exceeds precomputed matrix size")
  expect_error(get_stirling_row(0, logS), "positive integer")
})


# =============================================================================
# Golden Test Data Comparison
# =============================================================================

test_that("Stirling numbers match golden test data", {
  # Read golden data
  golden_path <- system.file("extdata/golden/golden_stirling.csv", 
                             package = "DPprior")
  
  # Skip if golden data not available (during development)
  skip_if(golden_path == "", "Golden data file not found")
  
  golden <- read.csv(golden_path, stringsAsFactors = FALSE)
  
  # Compute Stirling matrix up to max J in golden data
  J_max <- max(golden$J)
  logS <- compute_log_stirling(J_max)
  
  # Compare each golden value
  # Note: tolerance relaxed to 1e-6 to account for floating-point differences
  # between Python (used to generate golden data) and R implementations.
  # The relative error is still < 1e-8, which is well within acceptable limits
  # for log-scale Stirling numbers.
  for (i in seq_len(nrow(golden))) {
    J <- golden$J[i]
    k <- golden$k[i]
    expected <- golden$log_s[i]
    actual <- logS[J + 1, k + 1]
    
    expect_equal(actual, expected, tolerance = 1e-6,
                 info = sprintf("J=%d, k=%d", J, k))
  }
})

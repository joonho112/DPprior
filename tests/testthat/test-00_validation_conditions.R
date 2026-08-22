# =============================================================================
# Shared validation and typed-condition contract
# =============================================================================

test_that("invalid-input conditions expose stable classes and fields", {
  condition <- tryCatch(
    .dpprior_validate_scalar(NA_real_, "alpha"),
    dpprior_invalid_input = identity
  )

  expect_s3_class(condition, "dpprior_invalid_input")
  expect_s3_class(condition, "dpprior_scalar_error")
  expect_s3_class(condition, "dpprior_missing_error")
  expect_s3_class(condition, "error")
  expect_identical(condition$argument, "alpha")
  expect_identical(condition$expected, "non-missing")
  expect_identical(condition$code, "missing")
  expect_null(condition$call)
})


test_that("scalar validation distinguishes type, length, finiteness, and bounds", {
  expect_identical(.dpprior_validate_scalar(2, "x"), 2)
  expect_identical(.dpprior_validate_scalar(0, "x", 0, 1), 0)
  expect_error(
    .dpprior_validate_scalar("2", "x"),
    "x must be numeric", class = "dpprior_type_error"
  )
  expect_error(
    .dpprior_validate_scalar(c(1, 2), "x"),
    "x must have length 1", class = "dpprior_length_error"
  )
  expect_error(
    .dpprior_validate_scalar(NA_real_, "x"),
    "missing", class = "dpprior_missing_error"
  )
  expect_error(
    .dpprior_validate_scalar(Inf, "x"),
    "finite", class = "dpprior_nonfinite_error"
  )
  expect_error(
    .dpprior_validate_scalar(0, "x", lower = 0, lower_open = TRUE),
    "x must be > 0", class = "dpprior_bounds_error"
  )
  expect_error(
    .dpprior_validate_scalar(2, "x", lower = 0, upper = 1),
    "x must be in \\[0, 1\\]", class = "dpprior_bounds_error"
  )
})


test_that("probability validation supports scalar and vector contracts", {
  expect_identical(.dpprior_validate_probability(0), 0)
  expect_identical(.dpprior_validate_probability(1), 1)
  expect_equal(
    .dpprior_validate_probability(c(0.2, 0.8), scalar = FALSE),
    c(0.2, 0.8)
  )
  expect_error(
    .dpprior_validate_probability(c(0.2, 0.8)),
    class = "dpprior_length_error"
  )
  expect_error(
    .dpprior_validate_probability(NA_real_),
    class = "dpprior_missing_error"
  )
  expect_error(
    .dpprior_validate_probability(-Inf),
    class = "dpprior_nonfinite_error"
  )
  expect_error(
    .dpprior_validate_probability(1.1),
    "must be in", class = "dpprior_bounds_error"
  )
  expect_error(
    .dpprior_validate_probability(0, open = TRUE),
    "must be in", class = "dpprior_probability_error"
  )
})


test_that("count validation rejects non-integers and unsafe bounds", {
  expect_identical(.dpprior_validate_count(4, minimum = 1L), 4L)
  expect_error(.dpprior_validate_count("4"), class = "dpprior_type_error")
  expect_error(.dpprior_validate_count(c(1, 2)), class = "dpprior_length_error")
  expect_error(.dpprior_validate_count(NA_real_), class = "dpprior_missing_error")
  expect_error(.dpprior_validate_count(Inf), class = "dpprior_nonfinite_error")
  expect_error(.dpprior_validate_count(1.5), class = "dpprior_integer_error")
  expect_error(
    .dpprior_validate_count(0, minimum = 1L),
    class = "dpprior_bounds_error"
  )
  expect_error(
    .dpprior_validate_count(.Machine$integer.max + 1),
    class = "dpprior_bounds_error"
  )
})


test_that("PMF validation is strict and never silently normalizes", {
  pmf <- c(0.2, 0.3, 0.5)
  expect_identical(.dpprior_validate_pmf(pmf, expected_length = 3L), pmf)
  expect_identical(
    .dpprior_validate_pmf(c(2, 3), require_sum_one = FALSE),
    c(2, 3)
  )

  expect_error(.dpprior_validate_pmf("pmf"), class = "dpprior_type_error")
  expect_error(.dpprior_validate_pmf(numeric()), class = "dpprior_length_error")
  expect_error(
    .dpprior_validate_pmf(c(0.5, NA_real_)),
    class = "dpprior_missing_error"
  )
  expect_error(
    .dpprior_validate_pmf(c(0.5, Inf)),
    class = "dpprior_nonfinite_error"
  )
  expect_error(
    .dpprior_validate_pmf(c(0.5, 0.5), expected_length = 3L),
    class = "dpprior_length_error"
  )
  expect_error(
    .dpprior_validate_pmf(c(1.1, -0.1)),
    class = "dpprior_bounds_error"
  )
  expect_error(
    .dpprior_validate_pmf(c(0, 0)),
    class = "dpprior_pmf_mass_error"
  )
  expect_error(
    .dpprior_validate_pmf(c(2, 3)),
    "must sum to 1", class = "dpprior_normalization_error"
  )
})


test_that("control validation has one scalar contract across control types", {
  expect_identical(
    .dpprior_validate_control(1e-8, "tol", lower = 0, lower_open = TRUE),
    1e-8
  )
  expect_identical(
    .dpprior_validate_control(20, "M", type = "count", lower = 10),
    20L
  )
  expect_identical(
    .dpprior_validate_control(0.5, "p", type = "probability"),
    0.5
  )
  expect_identical(
    .dpprior_validate_control(
      0, "p", type = "probability",
      lower_open = FALSE, upper_open = TRUE
    ),
    0
  )
  expect_true(.dpprior_validate_control(TRUE, "verbose", type = "logical"))

  expect_error(
    .dpprior_validate_control("1e-8", "tol"),
    class = "dpprior_control_error"
  )
  expect_error(
    .dpprior_validate_control(c(1, 2), "tol"),
    class = "dpprior_length_error"
  )
  expect_error(
    .dpprior_validate_control(Inf, "tol"),
    class = "dpprior_nonfinite_error"
  )
  expect_error(
    .dpprior_validate_control(0, "tol", lower = 0, lower_open = TRUE),
    class = "dpprior_bounds_error"
  )
  expect_error(
    .dpprior_validate_control(NA, "verbose", type = "logical"),
    class = "dpprior_missing_error"
  )
  expect_error(
    .dpprior_validate_control(
      1, "p", type = "probability",
      lower_open = FALSE, upper_open = TRUE
    ),
    class = "dpprior_bounds_error"
  )
  expect_error(
    .dpprior_validate_control(
      2, "n", type = "count", lower = 2, lower_open = TRUE
    ),
    class = "dpprior_bounds_error"
  )
})


test_that("legacy assertions preserve messages and gain typed roots", {
  expect_error(
    assert_positive(0, "alpha"),
    "alpha must be finite and positive",
    class = "dpprior_invalid_input"
  )
  expect_error(
    assert_probability(NA_real_, "p"),
    "p must be finite and in \\[0, 1\\]",
    class = "dpprior_invalid_input"
  )
  expect_error(
    assert_valid_J(c(1, 2)),
    "J must be a positive integer",
    class = "dpprior_invalid_input"
  )
  expect_error(
    assert_positive(numeric(), "alpha"),
    "alpha must be finite and positive",
    class = "dpprior_length_error"
  )
  expect_error(
    assert_probability(numeric(), "p"),
    "p must be finite and in \\[0, 1\\]",
    class = "dpprior_length_error"
  )
  expect_error(
    assert_nonnegative(numeric(), "x"),
    "x must be finite and non-negative",
    class = "dpprior_length_error"
  )
  expect_error(assert_positive("bad"), class = "dpprior_type_error")
  expect_error(assert_positive(NA_real_), class = "dpprior_missing_error")
  expect_error(assert_positive(Inf), class = "dpprior_nonfinite_error")
  expect_error(assert_positive(-1), class = "dpprior_bounds_error")
  expect_identical(.as_integer_scalar(10, "M", min = 1L), 10L)
})


test_that("typed warning primitive uses the package warning hierarchy", {
  condition <- NULL
  withCallingHandlers(
    .dpprior_warn(
      "approximation used",
      subclass = "dpprior_approximation_warning",
      code = "approximation"
    ),
    warning = function(w) {
      condition <<- w
      invokeRestart("muffleWarning")
    }
  )
  expect_s3_class(condition, "dpprior_approximation_warning")
  expect_s3_class(condition, "dpprior_warning")
  expect_s3_class(condition, "warning")
  expect_identical(condition$code, "approximation")
})

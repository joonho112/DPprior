# =============================================================================
# Tests for Module 10: A1 Closed-Form Prior Elicitation
# =============================================================================
#
# Author: JoonHo Lee (jlee296@ua.edu)
# Date: December 2025
# Part of: DPprior R Package
# =============================================================================

# Note: context() is deprecated in testthat 3rd edition


# =============================================================================
# Test: compute_scaling_constant
# =============================================================================

test_that("log scaling returns log(J)", {
  expect_equal(compute_scaling_constant(50, "log"), log(50))
  expect_equal(compute_scaling_constant(100, "log"), log(100))
  expect_equal(compute_scaling_constant(2, "log"), log(2))
})

test_that("harmonic scaling returns H_{J-1} = digamma(J) + gamma", {
  euler_gamma <- 0.5772156649015329
  expect_equal(compute_scaling_constant(50, "harmonic"),
               digamma(50) + euler_gamma)
  expect_equal(compute_scaling_constant(10, "harmonic"),
               digamma(10) + euler_gamma)
})

test_that("digamma scaling requires mu_K argument", {
  expect_error(compute_scaling_constant(50, "digamma"),
               "mu_K required for digamma scaling")
  # Should succeed when mu_K is provided
  result <- compute_scaling_constant(50, "digamma", mu_K = 5)
  expect_true(is.finite(result))
  expect_true(result > 0)
})

test_that("compute_scaling_constant validates its documented public domain", {
  for (J in list(1, 2.5, NA_real_, "50")) {
    condition <- tryCatch(
      compute_scaling_constant(J, "log"),
      error = identity
    )
    expect_s3_class(condition, "dpprior_invalid_input")
    expect_s3_class(condition, "dpprior_a1_sample_size_error")
  }

  for (mu_K in list(NA_real_, Inf, c(2, 3), 1, 50)) {
    condition <- tryCatch(
      compute_scaling_constant(50, "digamma", mu_K),
      error = identity
    )
    expect_s3_class(condition, "dpprior_invalid_input")
    expect_s3_class(condition, "dpprior_a1_scaling_error")
  }
})

test_that("all scalings converge toward log(J) for large J", {
  J_large <- 500
  cJ_log <- compute_scaling_constant(J_large, "log")
  cJ_harm <- compute_scaling_constant(J_large, "harmonic")
  cJ_dig <- compute_scaling_constant(J_large, "digamma", mu_K = 5)

  # Harmonic and digamma should be within ~15% of log(J) for large J
  # The harmonic sum H_{J-1} converges to log(J) + gamma, so the relative

  # difference shrinks with J but does not vanish entirely
  expect_equal(cJ_harm / cJ_log, 1, tolerance = 0.15)
  expect_equal(cJ_dig / cJ_log, 1, tolerance = 0.25)
})

test_that("digamma scaling uses correct formula", {
  J <- 50
  mu_K <- 5
  alpha_tilde <- (mu_K - 1) / log(J)
  expected <- digamma(alpha_tilde + J) - digamma(alpha_tilde)
  expect_equal(compute_scaling_constant(J, "digamma", mu_K = mu_K), expected)
})


# =============================================================================
# Test: DPprior_a1 basic functionality
# =============================================================================

test_that("DPprior_a1 returns valid positive parameters", {
  fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)
  expect_true(fit$a > 0)
  expect_true(fit$b > 0)
})

test_that("DPprior_a1 returns correct structure", {
  fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)

  # Check all required components are present
  expect_true("a" %in% names(fit))
  expect_true("b" %in% names(fit))
  expect_true("J" %in% names(fit))
  expect_true("target" %in% names(fit))
  expect_true("method" %in% names(fit))
  expect_true("status" %in% names(fit))
  expect_true("scaling" %in% names(fit))
  expect_true("cJ" %in% names(fit))
  expect_true("var_K_used" %in% names(fit))
  expect_true("converged" %in% names(fit))
  expect_true("iterations" %in% names(fit))
  expect_true(all(c(
    "usable", "verified", "parameters", "achieved", "residuals",
    "tolerances", "attempts", "verification", "provenance"
  ) %in% names(fit)))

  # Check values
  expect_identical(class(fit), c("DPprior_fit", "dpprior_result", "list"))
  expect_identical(fit$schema, .dpprior_schema("result"))
  expect_identical(fit$object_type, "fit")
  expect_identical(fit$mode, "a1_proxy")
  expect_equal(fit$method, "A1")
  expect_equal(fit$J, 50)
  expect_identical(names(fit$target), "K")
  expect_identical(fit$target$K$kind, "moments")
  expect_equal(fit$target$K$request$mu_K, 5)
  expect_equal(fit$target$K$request$var_K, 8)
  expect_identical(fit$target$K$request, list(
    J = 50L, mu_K = 5, var_K = 8
  ))
  expect_identical(fit$target$K$normalized, list(
    J = 50L, mu_K = 5, var_K = 8, interval = NULL, pmf = NULL
  ))
  expect_identical(fit$target$K$used, fit$target$K$normalized)
  expect_identical(
    fit$target$K$derivation$request_to_normalized$rule,
    "canonicalize_direct_moments"
  )
  expect_null(fit$target$K$derivation$normalized_to_used)
  expect_identical(fit$status, "approximate")
  expect_true(fit$usable)
  expect_false(fit$verified)
  expect_false(fit$converged)
  expect_true(fit$mapping_verified)
  expect_identical(fit$iterations, 0L)
  expect_false(fit$verification$performed)
  expect_false(fit$verification$passed)
  expect_null(fit$verification$verifier_snapshot)
  expect_identical(
    fit$verification$reason,
    "not_performed_for_A1_proxy"
  )
  expect_identical(fit$computation$attempts, list())
  expect_identical(fit$computation$selected_attempt_id, NULL)
  expect_identical(fit$computation$termination$code, "closed_form")
  expect_identical(fit$computation$termination$source, "closed_form")
  expect_identical(fit$computation$termination$iterations, 0L)
  expect_identical(fit$computation$scaling$values$cJ, log(50))
  expect_true(fit$mapping_verification$performed)
  expect_true(fit$mapping_verification$passed)
  expect_true(all(fit$mapping_verification$component_pass))
  expect_true(all(
    abs(fit$mapping_verification$residuals) <=
      fit$mapping_verification$tolerances
  ))
  expect_identical(
    fit$mapping_verification$estimand,
    "shifted_negative_binomial_proxy_moments"
  )
  expect_identical(fit$parameters$a, fit$a)
  expect_identical(fit$parameters$b, fit$b)
  expect_identical(fit$proxy$mapping_verification,
                   fit$mapping_verification)
  expect_identical(
    fit$compatibility$deprecations$a1_v0$authority,
    "non_authoritative"
  )
  expect_false(fit$provenance$migration$lossless)
  expect_contains(
    fit$provenance$migration$missing_evidence,
    "independent_finite_J_verifier_not_performed"
  )
  expect_invisible(.dpprior_validate_result_v1(fit))
})

test_that("DPprior_a1 output has class DPprior_fit", {
  fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)
  expect_s3_class(fit, "DPprior_fit")
})

test_that("DPprior_a1 round-trip: NegBin moments recover targets", {
  fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)
  a <- fit$a
  b <- fit$b
  cJ <- fit$cJ

  # Forward NegBin model
  p <- b / (b + cJ)
  mu_S <- a * (1 - p) / p
  var_S <- a * (1 - p) / p^2

  mu_K_recovered <- mu_S + 1
  var_K_recovered <- var_S

  expect_equal(mu_K_recovered, 5, tolerance = 1e-10)
  expect_equal(var_K_recovered, 8, tolerance = 1e-10)
})

test_that("DPprior_a1 works with all scaling methods", {
  for (scaling in c("log", "harmonic", "digamma")) {
    fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 8, scaling = scaling)
    expect_true(fit$a > 0)
    expect_true(fit$b > 0)
    expect_equal(fit$scaling, scaling)
    expect_equal(fit$status, "approximate")
  }
})

test_that("DPprior_a1 different scalings give different parameters", {
  fit_log <- DPprior_a1(J = 50, mu_K = 5, var_K = 8, scaling = "log")
  fit_harm <- DPprior_a1(J = 50, mu_K = 5, var_K = 8, scaling = "harmonic")

  # Different cJ values should lead to different b (but same a)
  expect_false(identical(fit_log$cJ, fit_harm$cJ))
  expect_false(identical(fit_log$b, fit_harm$b))
  # Shape a should be identical since it depends only on mu_K and var_K
  expect_equal(fit_log$a, fit_harm$a, tolerance = 1e-10)
})


# =============================================================================
# Test: DPprior_a1 input validation
# =============================================================================

test_that("DPprior_a1 errors on invalid J", {
  expect_error(DPprior_a1(J = 1, mu_K = 5, var_K = 8),
               "J must be an integer >= 2")
  expect_error(DPprior_a1(J = 0, mu_K = 5, var_K = 8),
               "J must be an integer >= 2")
  expect_error(DPprior_a1(J = -5, mu_K = 5, var_K = 8),
               "J must be an integer >= 2")
  expect_error(DPprior_a1(J = 3.5, mu_K = 2, var_K = 3),
               "J must be an integer >= 2")
})

test_that("DPprior_a1 errors on invalid mu_K", {
  expect_error(DPprior_a1(J = 50, mu_K = 0.5, var_K = 8),
               "mu_K must be > 1")
  expect_error(DPprior_a1(J = 50, mu_K = 1, var_K = 8),
               "mu_K must be > 1")
  expect_error(DPprior_a1(J = 50, mu_K = 51, var_K = 8),
               "mu_K must be < J")
  expect_error(DPprior_a1(J = 50, mu_K = 50, var_K = 1),
               "mu_K must be < J")
})

test_that("DPprior_a1 errors on invalid var_K", {
  expect_error(DPprior_a1(J = 50, mu_K = 5, var_K = 0),
               "var_K must be a positive finite numeric scalar")
  expect_error(DPprior_a1(J = 50, mu_K = 5, var_K = -1),
               "var_K must be a positive finite numeric scalar")
  expect_error(DPprior_a1(J = 10, mu_K = 9, var_K = 9),
               "var_K = 9.*maximum possible variance 8.*K in \\{1,...,10\\}.*mu_K = 9")
})

test_that("DPprior_a1 public input failures have stable typed classes", {
  cases <- list(
    list(
      call = quote(DPprior_a1(J = 1, mu_K = 5, var_K = 8)),
      subclass = "dpprior_a1_sample_size_error"
    ),
    list(
      call = quote(DPprior_a1(J = 50, mu_K = 1, var_K = 8)),
      subclass = "dpprior_a1_mean_error"
    ),
    list(
      call = quote(DPprior_a1(J = 50, mu_K = 50, var_K = 1)),
      subclass = "dpprior_a1_mean_error"
    ),
    list(
      call = quote(DPprior_a1(J = 50, mu_K = 5, var_K = 0)),
      subclass = "dpprior_a1_variance_error"
    )
  )

  for (case in cases) {
    condition <- tryCatch(eval(case$call), error = identity)
    expect_s3_class(condition, "dpprior_invalid_input")
    expect_s3_class(condition, case$subclass)
    expect_s3_class(condition, "dpprior_bounds_error")
  }

  expect_error(
    DPprior_a1(J = 50, mu_K = 5, var_K = 8, scaling = "unknown"),
    class = "dpprior_a1_scaling_error"
  )
  expect_error(
    DPprior_a1(J = 50, mu_K = 5, var_K = 8, projection = "silent"),
    class = "dpprior_a1_projection_policy_error"
  )
})


# =============================================================================
# Test: DPprior_a1 variance feasibility and explicit projection
# =============================================================================

test_that("A1-infeasible target is not projected by default", {
  condition <- tryCatch(
    DPprior_a1(J = 50, mu_K = 5, var_K = 3),
    error = identity
  )

  expect_s3_class(condition, "dpprior_a1_projection_required")
  expect_s3_class(condition, "dpprior_a1_infeasible")
  expect_identical(condition$code, "a1_projection_required")
  expect_equal(condition$original_target, list(mu_K = 5, var_K = 3))
  expected_used <- 4 + .TOL_PROJECTION_BUFFER * (1 + 4^2)
  expect_equal(condition$projected_target$var_K, expected_used)
  expect_equal(condition$projection_distance, expected_used - 3)
  expect_identical(condition$projection_policy, "error")
})

test_that("explicit A1 projection emits exactly one typed condition", {
  warnings <- list()
  fit <- withCallingHandlers(
    DPprior_a1(
      J = 50, mu_K = 5, var_K = 3,
      projection = "nearest"
    ),
    warning = function(condition) {
      warnings[[length(warnings) + 1L]] <<- condition
      invokeRestart("muffleWarning")
    }
  )

  expect_length(warnings, 1L)
  expect_s3_class(warnings[[1L]], "dpprior_a1_projection_warning")
  expect_s3_class(warnings[[1L]], "dpprior_target_projection_warning")
  expect_identical(warnings[[1L]]$code, "a1_target_projected")
  expect_identical(fit$status, "approximate")
  expect_true(fit$usable)
  expect_false(fit$verified)
  expect_true(fit$mapping_verified)
  expect_true(fit$at_target_boundary)
  expect_false(fit$converged)
  expect_true(fit$projection$opt_in)
  expect_true(fit$projection$applied)
  expect_identical(fit$projection$policy, "nearest")
  expect_equal(fit$projection$original_target$var_K, 3)
  expect_equal(
    fit$projection$projected_target$var_K,
    fit$var_K_used
  )
  expect_equal(fit$projection$distance, fit$var_K_used - 3)
  expect_equal(fit$target$K$request$var_K, 3)
  expect_equal(fit$target$K$normalized$var_K, 3)
  expect_equal(fit$target$K$used$var_K, fit$var_K_used)
  expect_identical(
    fit$target$K$derivation$normalized_to_used$rule,
    "project_a1_variance_to_nearest_interior"
  )
  expect_identical(
    fit$target$K$derivation$normalized_to_used$before,
    fit$target$K$normalized
  )
  expect_identical(
    fit$target$K$derivation$normalized_to_used$after,
    fit$target$K$used
  )
  expect_identical(
    fit$provenance$projection,
    fit$target$K$provenance$projection
  )
  expect_invisible(.dpprior_validate_result_v1(fit))
  expect_true(fit$a > 0)
  expect_true(fit$b > 0)
})

test_that("feasible variance does not trigger projection", {
  fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)
  expect_equal(fit$status, "approximate")
  expect_false(fit$verified)
  expect_true(fit$mapping_verified)
  expect_false(fit$at_target_boundary)
  expect_equal(fit$var_K_used, 8)
  expect_false(fit$projection$applied)
  expect_identical(fit$projection$reason, "not_required")
})

test_that("var_K exactly at A1 boundary requires explicit projection", {
  condition <- tryCatch(
    DPprior_a1(J = 50, mu_K = 5, var_K = 4),
    error = identity
  )
  expect_s3_class(condition, "dpprior_a1_projection_required")

  warnings <- list()
  fit <- withCallingHandlers(
    DPprior_a1(
      J = 50, mu_K = 5, var_K = 4,
      projection = "nearest"
    ),
    warning = function(condition) {
      warnings[[length(warnings) + 1L]] <<- condition
      invokeRestart("muffleWarning")
    }
  )
  expect_length(warnings, 1L)
  expect_identical(fit$status, "approximate")
  expect_true(fit$at_target_boundary)
  expect_gt(fit$var_K_used, 4)
})

test_that("near-lower and support-upper A1 targets are visible boundaries", {
  epsilon <- .TOL_PROJECTION_BUFFER
  lower_near <- 4 + 0.5 * epsilon * (1 + 4^2)
  near <- DPprior_a1(50, 5, lower_near)
  upper <- DPprior_a1(10, 8, .max_var_K_fixed_mean(10, 8))

  expect_identical(near$status, "approximate")
  expect_true(near$at_target_boundary)
  expect_false(near$projection$applied)
  expect_identical(near$projection$reason, "near_a1_lower_boundary")
  expect_equal(near$var_K_used, lower_near)
  expect_identical(upper$status, "approximate")
  expect_true(upper$at_target_boundary)
  expect_identical(upper$projection$reason, "support_upper_boundary")
})

test_that("A1 shape caveats never mutate the closed status", {
  fit <- NULL
  expect_no_warning(
    fit <- DPprior_a1(J = 50, mu_K = 1.1, var_K = 4)
  )

  expect_identical(fit$status, "approximate")
  expect_false(grepl("warning|\\[|\\]", fit$status))
  expect_contains(fit$caveats, "quasi_improper_shape")
})

test_that("projected target is rechecked against fixed-support upper bound", {
  warnings <- list()
  condition <- tryCatch(
    withCallingHandlers(
      DPprior_a1(
        J = 10, mu_K = 9.5, var_K = 4,
        projection = "nearest"
      ),
      warning = function(w) {
        warnings[[length(warnings) + 1L]] <<- w
        invokeRestart("muffleWarning")
      }
    ),
    error = identity
  )

  expect_length(warnings, 0L)
  expect_s3_class(condition, "dpprior_a1_projection_impossible")
  expect_s3_class(condition, "dpprior_a1_infeasible")
  expect_identical(condition$code, "a1_projection_exceeds_support")
  expect_lte(condition$support_upper_bound, condition$a1_lower_bound)
  expect_lte(condition$available_gap, 0)
})

test_that("A1 adaptive buffer preserves a narrow feasible projection interval", {
  J <- 500L
  mu_K <- 498.9999
  lower <- mu_K - 1
  upper <- (mu_K - 1) * (J - mu_K)
  request <- lower - 0.01

  denied <- tryCatch(
    DPprior_a1(J, mu_K, request),
    error = identity
  )
  expect_s3_class(denied, "dpprior_a1_projection_required")
  expect_gt(denied$projected_target$var_K, lower)
  expect_lte(denied$projected_target$var_K, upper)

  warnings <- list()
  fit <- withCallingHandlers(
    DPprior_a1(J, mu_K, request, projection = "nearest"),
    warning = function(condition) {
      warnings[[length(warnings) + 1L]] <<- condition
      invokeRestart("muffleWarning")
    }
  )
  expect_length(warnings, 1L)
  expect_s3_class(warnings[[1L]], "dpprior_a1_projection_warning")
  expect_gt(fit$var_K_used, lower)
  expect_lte(fit$var_K_used, upper)
  expect_true(fit$projection$buffer_was_capped)
  expect_gt(fit$projection$requested_buffer,
            fit$projection$effective_buffer)
  expect_equal(fit$projection$effective_buffer,
               fit$projection$available_gap / 2)
  expect_true(fit$projection$representable_interior)
})

test_that("A1 tiny epsilon is raised to a representable interior step", {
  warnings <- list()
  fit <- withCallingHandlers(
    DPprior_a1(
      J = 50, mu_K = 5, var_K = 3,
      projection = "nearest", epsilon = 1e-20
    ),
    warning = function(condition) {
      warnings[[length(warnings) + 1L]] <<- condition
      invokeRestart("muffleWarning")
    }
  )

  expect_length(warnings, 1L)
  expect_s3_class(warnings[[1L]], "dpprior_a1_projection_warning")
  expect_gt(fit$var_K_used, 4)
  expect_lte(fit$var_K_used, fit$projection$support_upper_bound)
  expect_true(fit$projection$buffer_was_floored)
  expect_gt(fit$projection$effective_buffer,
            fit$projection$requested_buffer)
  expect_gte(fit$projection$effective_buffer,
             fit$projection$representability_floor)
  expect_true(fit$projection$representable_interior)
})


# =============================================================================
# Test: DPprior_a1 closed-form algebra
# =============================================================================

test_that("A1 closed-form formulas are algebraically correct", {
  J <- 50
  mu_K <- 5
  var_K <- 8

  cJ <- log(J)
  m <- mu_K - 1       # shifted mean
  D <- var_K - m       # denominator

  a_expected <- m^2 / D
  b_expected <- m * cJ / D

  fit <- DPprior_a1(J = J, mu_K = mu_K, var_K = var_K, scaling = "log")
  expect_equal(fit$a, a_expected, tolerance = 1e-12)
  expect_equal(fit$b, b_expected, tolerance = 1e-12)
})

test_that("A1 canonical producer preserves the frozen numerical baseline", {
  cases <- list(
    log = log(50),
    harmonic = digamma(50) + .EULER_GAMMA,
    digamma = {
      alpha_tilde <- 4 / log(50)
      digamma(alpha_tilde + 50) - digamma(alpha_tilde)
    }
  )

  for (scaling in names(cases)) {
    fit <- DPprior_a1(50, 5, 8, scaling = scaling)
    expect_identical(fit$parameters$a, 4)
    expect_equal(fit$parameters$b, cases[[scaling]], tolerance = 0)
    expect_equal(
      fit$computation$scaling$values$cJ,
      cases[[scaling]], tolerance = 0
    )
    expect_identical(
      fit$target$K$implied,
      list(mean = 5, variance = 8)
    )
    expect_equal(fit$achieved$K$mean, 5, tolerance = 1e-12)
    expect_equal(fit$achieved$K$variance, 8, tolerance = 1e-12)
    expect_true(all(
      abs(unlist(fit$residuals$K, use.names = TRUE)) <=
        fit$tolerances$mapping
    ))
  }
})

test_that("A1 canonical schema rejects authority mutations", {
  fit <- DPprior_a1(50, 5, 8)

  wrong_mode <- fit
  wrong_mode$mode <- "a2_moment"
  expect_error(
    .dpprior_validate_result_v1(wrong_mode),
    class = "dpprior_schema_error"
  )

  forged_target <- fit
  forged_target$target$K$used$var_K <- 9
  expect_error(
    .dpprior_validate_result_v1(forged_target),
    class = "dpprior_schema_error"
  )

  forged_verification <- fit
  forged_verification$verified <- TRUE
  expect_error(
    .dpprior_validate_result_v1(forged_verification),
    class = "dpprior_schema_error"
  )

  forged_alias <- fit
  forged_alias$a <- forged_alias$a + 1
  expect_error(
    .dpprior_validate_result_v1(forged_alias),
    class = "dpprior_schema_error"
  )

  projected <- suppressWarnings(DPprior_a1(
    50, 5, 3, projection = "nearest"
  ))
  projected$provenance$projection$record$after$var_K <-
    projected$provenance$projection$record$after$var_K + 1e-6
  expect_error(
    .dpprior_validate_result_v1(projected),
    class = "dpprior_schema_error"
  )
})

test_that("A1 defensive mapping failure is retained but not usable", {
  resolution <- .dpprior_a1_resolve_target(
    50, 5, 8, projection = "error", signal_projection = FALSE
  )
  target <- c(mean = 5, variance = 8)
  achieved <- c(mean = 5, variance = 9)
  residual <- achieved - target
  tolerance <- stats::setNames(
    1e-12 + 1e-10 * pmax(1, abs(target)), names(target)
  )
  component_pass <- abs(residual) <= tolerance

  fit <- .dpprior_a1_result_v1(
    J = 50,
    mu_K = 5,
    var_K_requested = 8,
    var_K_used = 8,
    scaling = "log",
    cJ = log(50),
    epsilon = .TOL_PROJECTION_BUFFER,
    projection_policy = "error",
    target_resolution = resolution,
    a = 4,
    b = log(50),
    mapping_target = target,
    mapping_achieved = achieved,
    mapping_residual = residual,
    mapping_tolerance = tolerance,
    mapping_component_pass = component_pass,
    mapping_passed = FALSE,
    caveats = character(),
    message = "A1 proxy round-trip failed."
  )

  expect_identical(fit$status, "approximate")
  expect_false(fit$usable)
  expect_false(fit$verified)
  expect_identical(fit$parameters$a, 4)
  expect_identical(fit$parameters$b, log(50))
  expect_identical(fit$achieved$K$variance, 9)
  expect_false(is.null(fit$verification$selected_snapshot))
  expect_identical(fit$computation$termination$code, "closed_form")
  expect_identical(fit$computation$termination$source, "closed_form")
  expect_true("a" %in% names(fit))
  expect_identical(fit$compatibility$views$a1_v0$a, 4)
  expect_false(fit$proxy$mapping_verification$passed)
  expect_invisible(.dpprior_validate_result_v1(fit))
})

test_that("A1 canonical result round-trips through base serialization", {
  fit <- suppressWarnings(DPprior_a1(
    50, 5, 3, scaling = "harmonic", projection = "nearest"
  ))
  raw_roundtrip <- unserialize(serialize(fit, NULL, version = 3L))
  expect_identical(raw_roundtrip, fit)
  expect_invisible(.dpprior_validate_result_v1(raw_roundtrip))

  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(fit, path, version = 3L)
  disk_roundtrip <- readRDS(path)
  expect_identical(disk_roundtrip, fit)
  expect_invisible(.dpprior_validate_result_v1(disk_roundtrip))
})


# =============================================================================
# Test: vif_to_variance
# =============================================================================

test_that("VIF = 1 gives Poisson baseline variance (mu_K - 1)", {
  expect_equal(vif_to_variance(mu_K = 5, vif = 1), 4)
  expect_equal(vif_to_variance(mu_K = 10, vif = 1), 9)
})

test_that("VIF = 2 gives twice the Poisson baseline variance", {
  expect_equal(vif_to_variance(mu_K = 5, vif = 2), 8)
  expect_equal(vif_to_variance(mu_K = 10, vif = 2), 18)
})

test_that("VIF < 1 errors", {
  expect_error(vif_to_variance(mu_K = 5, vif = 0.5),
               "vif must be >= 1")
  expect_error(vif_to_variance(mu_K = 5, vif = 0),
               "vif must be >= 1")
})

test_that("vif_to_variance errors on invalid mu_K", {
  expect_error(vif_to_variance(mu_K = 0.5, vif = 2),
               "mu_K must be > 1")
})

test_that("vif_to_variance finite-input failures are typed", {
  missing_condition <- tryCatch(
    vif_to_variance(mu_K = 5, vif = NA_real_),
    error = identity
  )
  expect_s3_class(missing_condition, "dpprior_invalid_input")
  expect_s3_class(missing_condition, "dpprior_a1_vif_error")
  expect_s3_class(missing_condition, "dpprior_missing_error")

  nonfinite_condition <- tryCatch(
    vif_to_variance(mu_K = 5, vif = Inf),
    error = identity
  )
  expect_s3_class(nonfinite_condition, "dpprior_invalid_input")
  expect_s3_class(nonfinite_condition, "dpprior_a1_vif_error")
  expect_s3_class(nonfinite_condition, "dpprior_nonfinite_error")
})


# =============================================================================
# Test: confidence_to_vif
# =============================================================================

test_that("confidence levels map to decreasing VIF (low > medium > high)", {
  vif_low <- confidence_to_vif("low")
  vif_med <- confidence_to_vif("medium")
  vif_high <- confidence_to_vif("high")

  expect_true(vif_low > vif_med)
  expect_true(vif_med > vif_high)
})

test_that("confidence_to_vif returns known values", {
  expect_equal(confidence_to_vif("low"), 5.0)
  expect_equal(confidence_to_vif("medium"), 2.5)
  expect_equal(confidence_to_vif("high"), 1.5)
})

test_that("confidence_to_vif returns scalar numeric", {
  result <- confidence_to_vif("medium")
  expect_true(is.numeric(result))
  expect_length(result, 1)
})

test_that("confidence_to_vif rejects invalid input", {
  expect_error(
    confidence_to_vif("extreme"),
    class = "dpprior_a1_confidence_error"
  )
  expect_error(confidence_to_vif("extreme"), class = "dpprior_invalid_input")
  expect_error(confidence_to_vif("extreme"), class = "dpprior_choice_error")
})


# =============================================================================
# Test: cv_alpha_to_variance
# =============================================================================

test_that("cv_alpha_to_variance follows correct algebra", {
  mu_K <- 5
  cv_alpha <- 0.5
  m <- mu_K - 1
  expected <- m + (cv_alpha * m)^2   # m(1 + cv^2 * m) = m + cv^2 * m^2
  expect_equal(cv_alpha_to_variance(mu_K, cv_alpha), expected)
})

test_that("larger cv_alpha gives larger variance", {
  var1 <- cv_alpha_to_variance(mu_K = 5, cv_alpha = 0.5)
  var2 <- cv_alpha_to_variance(mu_K = 5, cv_alpha = 1.0)
  var3 <- cv_alpha_to_variance(mu_K = 5, cv_alpha = 2.0)

  expect_true(var1 < var2)
  expect_true(var2 < var3)
})

test_that("cv_alpha round-trip: recovered CV matches target", {
  mu_K <- 5
  cv_target <- 0.5

  var_K <- cv_alpha_to_variance(mu_K, cv_target)
  fit <- DPprior_a1(J = 50, mu_K = mu_K, var_K = var_K)
  cv_recovered <- 1 / sqrt(fit$a)

  expect_equal(cv_recovered, cv_target, tolerance = 1e-8)
})

test_that("cv_alpha_to_variance errors on invalid inputs", {
  expect_error(cv_alpha_to_variance(mu_K = 0.5, cv_alpha = 1),
               "mu_K must be > 1")
  expect_error(cv_alpha_to_variance(mu_K = 5, cv_alpha = -0.5),
               "cv_alpha must be positive")
  expect_error(cv_alpha_to_variance(mu_K = 5, cv_alpha = 0),
               "cv_alpha must be positive")
})

test_that("cv_alpha_to_variance finite-input failures are typed", {
  missing_condition <- tryCatch(
    cv_alpha_to_variance(mu_K = 5, cv_alpha = NA_real_),
    error = identity
  )
  expect_s3_class(missing_condition, "dpprior_invalid_input")
  expect_s3_class(missing_condition, "dpprior_a1_cv_error")
  expect_s3_class(missing_condition, "dpprior_missing_error")

  nonfinite_condition <- tryCatch(
    cv_alpha_to_variance(mu_K = 5, cv_alpha = Inf),
    error = identity
  )
  expect_s3_class(nonfinite_condition, "dpprior_invalid_input")
  expect_s3_class(nonfinite_condition, "dpprior_a1_cv_error")
  expect_s3_class(nonfinite_condition, "dpprior_nonfinite_error")
})


# =============================================================================
# Test: compare_a1_a2
# =============================================================================

test_that("compare_a1_a2 returns list with a1 and a2 components", {
  result <- compare_a1_a2(J = 50, mu_K = 5, var_K = 8, verbose = FALSE)

  expect_true(is.list(result))
  expect_true("a1" %in% names(result))
  expect_true("a2" %in% names(result))
  expect_true("improvement_ratio" %in% names(result))

  # Both components have expected fields
  expect_true("a" %in% names(result$a1))
  expect_true("b" %in% names(result$a1))
  expect_true("residual" %in% names(result$a1))
  expect_true("a" %in% names(result$a2))
  expect_true("b" %in% names(result$a2))
  expect_true("residual" %in% names(result$a2))
})

test_that("A2 achieves smaller or equal residual than A1", {
  result <- compare_a1_a2(J = 50, mu_K = 5, var_K = 8, verbose = FALSE)
  expect_true(result$a2$residual <= result$a1$residual)
  expect_true(result$improvement_ratio >= 1)
})

test_that("compare_a1_a2 works for various parameter combos", {
  # Moderate case

  result <- compare_a1_a2(J = 100, mu_K = 10, var_K = 20, verbose = FALSE)
  expect_true(is.list(result))
  expect_true(result$a1$a > 0)
  expect_true(result$a2$a > 0)
})


# =============================================================================
# Test: S3 methods (print, summary, as.data.frame)
# =============================================================================

.a1_without_top_level_aliases <- function(fit) {
  raw <- unclass(fit)
  aliases <- names(raw$compatibility$top_level_aliases)
  raw[aliases] <- NULL
  raw$compatibility$top_level_aliases <-
    stats::setNames(character(), character())
  class(raw) <- class(fit)
  .dpprior_validate_result_v1(raw)
  raw
}

test_that("print.DPprior_fit runs without error", {
  fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)
  expect_output(print(fit), "DPprior")
  expect_output(print(fit), "Method: A1")
})

test_that("summary.DPprior_fit returns correct structure", {
  fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)
  s <- summary(fit, print_output = FALSE)

  expect_true(is.list(s))
  expect_equal(s$method, "A1")
  expect_equal(s$gamma_prior$a, fit$a)
  expect_equal(s$gamma_prior$b, fit$b)
  # Module 17 summary uses E_alpha and CV_alpha field names
  expect_equal(s$alpha_summary$E_alpha, fit$a / fit$b)
  expect_equal(s$alpha_summary$CV_alpha, 1 / sqrt(fit$a))
  expect_false(s$converged)
  expect_equal(s$iterations, 0L)
})

test_that("as.data.frame.DPprior_fit returns one-row data.frame", {
  fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)
  df <- as.data.frame(fit)

  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 1)
  expect_equal(df$method, "A1")
  expect_equal(df$a, fit$a)
  expect_equal(df$b, fit$b)
  expect_equal(df$J, 50)
  expect_equal(df$mu_K, 5)
  expect_equal(df$var_K, 8)
})

test_that("A1 S3 and round-trip consumers use canonical nested authority", {
  fit <- .a1_without_top_level_aliases(DPprior_a1(50, 5, 8))

  expect_false(any(c("a", "b", "cJ", "var_K_used") %in% names(fit)))
  expect_output(print(fit), "Method: A1")
  summary_fit <- summary(fit, print_output = FALSE)
  expect_identical(summary_fit$gamma_prior$a, fit$parameters$a)
  expect_identical(summary_fit$gamma_prior$b, fit$parameters$b)
  frame <- as.data.frame(fit)
  expect_equal(frame$a, fit$parameters$a)
  expect_equal(frame$b, fit$parameters$b)
  expect_equal(frame$mu_K, fit$target$K$request$mu_K)
  expect_equal(frame$var_K, fit$target$K$used$var_K)
  expect_true(verify_a1_roundtrip(fit, verbose = FALSE))
})


# =============================================================================
# Test: verify_a1_roundtrip (internal verification function)
# =============================================================================

test_that("verify_a1_roundtrip passes for feasible cases", {
  fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)
  expect_true(verify_a1_roundtrip(fit, verbose = FALSE))

  fit2 <- DPprior_a1(J = 100, mu_K = 10, var_K = 20)
  expect_true(verify_a1_roundtrip(fit2, verbose = FALSE))
})

test_that("verify_a1_roundtrip passes for projected cases", {
  fit <- suppressWarnings(DPprior_a1(
    J = 50, mu_K = 5, var_K = 3,
    projection = "nearest"
  ))
  # Round-trip should still pass using var_K_used (projected value)
  expect_true(verify_a1_roundtrip(fit, verbose = FALSE))
})

test_that("A1 choice inputs reject dimensions and custom classes", {
  malformed <- list(
    matrix("error", nrow = 1L),
    array("nearest", dim = c(1L, 1L, 1L)),
    structure("error", class = "dpprior_test_character")
  )
  for (value in malformed) {
    expect_error(
      DPprior_a1(20, 5, 8, projection = value),
      class = "dpprior_type_error"
    )
  }
  expect_error(
    compute_scaling_constant(20, matrix("log", nrow = 1L)),
    class = "dpprior_type_error"
  )
})

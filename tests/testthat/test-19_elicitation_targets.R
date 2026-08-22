# =============================================================================
# Module 19: bounded-discrete elicitation targets
# =============================================================================

test_that("interval validation rejects ambiguous and malformed requests", {
  bare <- tryCatch(
    .dp_validate_K_interval(c(3, 10), J = 20),
    dpprior_ambiguous_interval_error = identity
  )
  expect_s3_class(bare, "dpprior_ambiguous_interval_error")
  expect_s3_class(bare, "dpprior_invalid_input")
  expect_identical(bare$code, "bare_range")
  expect_identical(bare$inferred_variance, NULL)

  expect_error(
    .dp_validate_K_interval(
      list(lower = 3L, upper = 10L, type = "plausible_range",
           family = "maxent"),
      J = 20
    ),
    "ambiguous", class = "dpprior_ambiguous_interval_error"
  )
  plausible_without_family <- tryCatch(
    .dp_validate_K_interval(
      list(lower = 3L, upper = 10L, type = "plausible"), J = 20
    ),
    dpprior_ambiguous_interval_error = identity
  )
  expect_s3_class(
    plausible_without_family, "dpprior_ambiguous_interval_error"
  )
  expect_identical(plausible_without_family$code, "plausible_range")
  expect_error(
    .dp_validate_K_interval(
      list(lower = 3, upper = 10, type = "central_mass",
           family = "maxent"),
      J = 20
    ),
    "coverage", class = "dpprior_interval_missing_coverage_error"
  )
  expect_error(
    .dp_validate_K_interval(
      list(lower = 3, upper = 10, type = "hard_bounds",
           coverage = 1, family = "maxent"),
      J = 20
    ),
    class = "dpprior_interval_semantic_conflict_error"
  )
  expect_error(
    .dp_validate_K_interval(
      list(lower = 3.5, upper = 10, type = "hard_bounds",
           family = "maxent"),
      J = 20
    ),
    class = "dpprior_interval_support_error"
  )
  expect_error(
    .dp_validate_K_interval(
      list(lower = 3L, upper = 10L, type = "equal_tail",
           coverage = 1, family = "maxent"),
      J = 20
    ),
    "hard bounds", class = "dpprior_interval_semantic_conflict_error"
  )
  expect_error(
    .dp_validate_K_interval(
      list(lower = 3L, upper = 10L, type = "equal_tail",
           coverage = 0, family = "maxent"),
      J = 20
    ),
    class = "dpprior_interval_probability_error"
  )
})

test_that("canonical moment adapters preserve assumptions without inventing a PMF", {
  moments <- .dp_target_K(J = 20, mu_K = 6, var_K = 5)
  expect_s3_class(moments, "dpprior_K_target")
  expect_identical(moments$schema, .dpprior_schema("target"))
  expect_identical(moments$kind, "moments")
  expect_true(.dpprior_validate_target_v1(moments, collect = TRUE)$valid)
  expect_identical(moments$status, "converged")
  expect_true(moments$usable)
  expect_null(moments$pmf)
  expect_equal(moments$implied$mean, 6)
  expect_equal(moments$implied$variance, 5)
  expect_identical(
    moments$compatibility$views$target_v0$provenance$pmf_identified,
    FALSE
  )
  expect_identical(moments$attempts, moments$computation$attempts)
  expect_identical(moments$constraint_residuals, moments$residuals)

  confidence <- .dp_target_K(J = 20, mu_K = 6, confidence = "medium")
  expect_null(confidence$pmf)
  expect_equal(confidence$implied$variance, 2.5 * (6 - 1))
  expect_identical(confidence$request$confidence, "medium")

  cv <- .dp_target_K(J = 20, mu_K = 6, cv_K = 0.25)
  expect_null(cv$pmf)
  expect_equal(cv$implied$variance, (0.25 * 6)^2)
  expect_equal(
    cv$compatibility$views$target_v0$provenance$derived_variance,
    2.25
  )

  expect_error(
    .dp_target_K(J = 20, mu_K = 6, var_K = 5, cv_K = 0.25),
    class = "dpprior_conflicting_input"
  )
  expect_error(
    .dp_target_K(J = 20, mu_K = 6, cv_K = 10),
    class = "dpprior_moment_feasibility_error"
  )
  expect_error(
    .dp_target_K(
      J = 20, mu_K = 6, var_K = 5,
      root_control = list(unknown = 1)
    ),
    class = "dpprior_unknown_argument_error"
  )
})

test_that("strict custom PMF is authoritative and never normalized", {
  pmf <- c(0.1, 0.2, 0.3, 0.4)
  target <- .dp_target_K(J = 4, target_pmf = pmf)
  expect_identical(target$pmf, pmf)
  expect_equal(target$implied$mean, sum((1:4) * pmf))
  expect_null(target$family)
  expect_identical(
    target$compatibility$views$target_v0$family$name, "custom_pmf"
  )
  expect_identical(
    target$compatibility$views$target_v0$provenance$normalization,
    "validated_not_modified"
  )

  expect_error(
    .dp_target_K(J = 4, target_pmf = c(1, 2, 3, 4)),
    class = "dpprior_normalization_error"
  )

  legacy <- .dp_target_K(J = 4, target_pmf = c(0, pmf))
  expect_identical(legacy$pmf, pmf)
  expect_true(
    legacy$compatibility$views$target_v0$provenance$k0_entry_removed
  )

  matrix_input <- matrix(rep(0.25, 4), nrow = 2L)
  matrix_error <- tryCatch(
    .dp_target_K(J = 4, target_pmf = matrix_input),
    error = identity
  )
  expect_s3_class(matrix_error, "dpprior_pmf_error")
  expect_s3_class(matrix_error, "dpprior_type_error")
  expect_identical(matrix_error$code, "vector_required")

  legacy_matrix <- matrix(c(0, 0.5, 0.5), nrow = 3L)
  legacy_matrix_error <- tryCatch(
    .dp_target_K(J = 2, target_pmf = legacy_matrix),
    error = identity
  )
  expect_s3_class(legacy_matrix_error, "dpprior_pmf_error")
  expect_identical(legacy_matrix_error$code, "vector_required")

  legacy_array <- array(c(0, 0.5, 0.5), dim = c(1L, 3L, 1L))
  legacy_array_error <- tryCatch(
    .dp_target_K(J = 2, target_pmf = legacy_array),
    error = identity
  )
  expect_s3_class(legacy_array_error, "dpprior_pmf_error")
  expect_identical(legacy_array_error$code, "vector_required")
})

test_that("elicitation tolerance cannot undercut the numerical error budget", {
  condition <- tryCatch(
    .dp_target_K(10, target_pmf = rep(0.1, 10), tolerance = 1e-300),
    error = identity
  )
  expect_s3_class(condition, "dpprior_bounds_error")
  expect_identical(condition$argument, "tolerance")
})

test_that("explicit maxent hard bounds reproduce the reviewer calculation", {
  interval <- list(
    lower = 3L, upper = 10L, type = "hard_bounds", family = "maxent"
  )
  target <- .dp_target_K(J = 20, K_interval = interval)

  expect_identical(target$status, "converged")
  expect_true(target$verified)
  expect_equal(target$pmf[3:10], rep(1 / 8, 8), tolerance = 1e-14)
  expect_equal(target$pmf[-(3:10)], numeric(12), tolerance = 0)
  expect_equal(target$implied$mean, 6.5, tolerance = 1e-14)
  expect_equal(target$implied$variance, 5.25, tolerance = 1e-14)
  expect_equal(target$achieved_interval$inside_mass, 1, tolerance = 1e-14)

  singleton <- .dp_target_K(
    J = 20,
    K_interval = list(
      lower = 5L, upper = 5L, type = "hard_bounds", family = "maxent"
    )
  )
  expect_identical(singleton$status, "boundary")
  expect_equal(singleton$pmf, replace(numeric(20), 5, 1), tolerance = 0)
})

test_that("maxent central mass uses analytic group masses and common tilt", {
  request <- list(
    lower = 3L, upper = 10L, type = "central_mass",
    coverage = 0.8, family = "maxent"
  )
  target <- .dp_target_K(J = 20, mu_K = 6.5, K_interval = request)
  support <- 1:20

  expect_identical(target$status, "converged")
  expect_true(target$usable)
  expect_true(target$verified)
  expect_equal(sum(support * target$pmf), 6.5, tolerance = 1e-9)
  expect_equal(sum(target$pmf[3:10]), 0.8, tolerance = 1e-12)
  expect_equal(sum(target$pmf), 1, tolerance = 1e-14)
  expect_true(all(target$pmf >= 0))
  expect_identical(target$family$name, "maxent")

  expect_error(
    .dp_target_K(J = 20, K_interval = request),
    class = "dpprior_interval_missing_mean_error"
  )
})

test_that("maxent equal tails enforces distinct two-tail semantics", {
  request <- list(
    lower = 3L, upper = 10L, type = "equal_tail",
    coverage = 0.8, family = "maxent"
  )
  target <- .dp_target_K(J = 20, K_interval = request)

  expect_identical(target$status, "converged")
  expect_equal(sum(target$pmf[1:2]), 0.1, tolerance = 1e-14)
  expect_equal(sum(target$pmf[3:10]), 0.8, tolerance = 1e-14)
  expect_equal(sum(target$pmf[11:20]), 0.1, tolerance = 1e-14)

  tilted <- .dp_target_K(J = 20, mu_K = 6.5, K_interval = request)
  expect_equal(tilted$implied$mean, 6.5, tolerance = 1e-9)
  expect_equal(tilted$achieved_interval$left_mass, 0.1, tolerance = 1e-12)
  expect_equal(tilted$achieved_interval$right_mass, 0.1, tolerance = 1e-12)
})

test_that("analytic infeasibility is returned without projection", {
  empty_tail <- .dp_target_K(
    J = 10,
    K_interval = list(
      lower = 3L, upper = 10L, type = "equal_tail",
      coverage = 0.8, family = "maxent"
    )
  )
  expect_identical(empty_tail$status, "infeasible")
  expect_false(empty_tail$usable)
  expect_true(empty_tail$verified)
  expect_true(empty_tail$verification$passed)
  expect_false(
    empty_tail$compatibility$views$target_v0$verification$constraint_feasible
  )
  expect_true(
    empty_tail$compatibility$views$target_v0$verification$
      infeasibility_certified
  )
  expect_null(empty_tail$pmf)
  expect_true(
    "dpprior_interval_empty_tail_error" %in%
      empty_tail$compatibility$views$target_v0$condition$classes
  )
  expect_error(
    .dp_target_K_stop_unusable(empty_tail),
    class = "dpprior_interval_empty_tail_error"
  )

  impossible_mean <- .dp_target_K(
    J = 20, mu_K = 10.3,
    K_interval = list(
      lower = 3L, upper = 10L, type = "equal_tail",
      coverage = 0.8, family = "maxent"
    )
  )
  expect_identical(impossible_mean$status, "infeasible")
  expect_identical(
    impossible_mean$compatibility$views$target_v0$provenance$projection_used,
    FALSE
  )
  expect_true(
    "dpprior_interval_mean_infeasible" %in%
      impossible_mean$compatibility$views$target_v0$condition$classes
  )
  mean_condition <- tryCatch(
    .dp_target_K_stop_unusable(impossible_mean), error = identity
  )
  expect_s3_class(mean_condition, "dpprior_interval_mean_infeasible")
  expect_identical(mean_condition$code, "mean_outside_group_mass_hull")
  expect_identical(mean_condition$result, impossible_mean)

  exact_boundary <- .dp_target_K(
    J = 20, mu_K = 3.6,
    K_interval = list(
      lower = 3L, upper = 10L, type = "equal_tail",
      coverage = 0.8, family = "maxent"
    )
  )
  expect_identical(exact_boundary$status, "boundary")
  expect_equal(exact_boundary$implied$mean, 3.6, tolerance = 1e-14)

  outside_by_small_amount <- .dp_target_K(
    J = 20, mu_K = 3.6 - 1e-10,
    K_interval = list(
      lower = 3L, upper = 10L, type = "equal_tail",
      coverage = 0.8, family = "maxent"
    )
  )
  expect_identical(outside_by_small_amount$status, "infeasible")
  expect_null(outside_by_small_amount$pmf)
  expect_identical(
    outside_by_small_amount$compatibility$views$target_v0$provenance$
      projection_used,
    FALSE
  )
})

test_that("maxent numerical solver failures are not certified infeasibility", {
  request <- list(
    lower = 3L, upper = 10L, type = "central_mass",
    coverage = 0.8, family = "maxent"
  )
  nonconverged <- .dp_target_K(
    J = 20, mu_K = 6.5, K_interval = request,
    root_control = list(max_iterations = 1L)
  )
  expect_identical(nonconverged$status, "approximate")
  expect_false(nonconverged$verified)
  expect_identical(nonconverged$attempts[[1L]]$exit_code, 1L)
  expect_false(
    nonconverged$compatibility$views$target_v0$attempts[[1L]]$converged
  )
  expect_true(length(nonconverged$attempts[[1L]]$warnings) >= 1L)

  testthat::local_mocked_bindings(
    .dp_interval_group_pmf = function(J, groups, tilt = 0,
                                      boundary = NULL) rep(1 / J, J),
    .package = "DPprior"
  )
  bracket_failed <- .dp_target_K(
    J = 20, mu_K = 6.5, K_interval = request
  )
  expect_identical(bracket_failed$status, "failed")
  expect_false(bracket_failed$verified)
  expect_false(bracket_failed$verification$passed)
  expect_true(
    bracket_failed$compatibility$views$target_v0$verification$
      constraint_feasible
  )
  expect_false(
    bracket_failed$compatibility$views$target_v0$verification$
      infeasibility_certified
  )
  expect_true(
    "dpprior_interval_solver_failed" %in%
      bracket_failed$compatibility$views$target_v0$condition$classes
  )
})

test_that("positive interval masses cannot disappear under floating-point resolution", {
  smallest_positive <- .Machine$double.xmin * .Machine$double.eps
  request <- list(
    lower = 3L, upper = 10L, type = "central_mass",
    coverage = smallest_positive, family = "maxent"
  )
  target <- .dp_target_K(
    J = 20, mu_K = 6.5, K_interval = request
  )

  expect_identical(target$status, "approximate")
  expect_false(target$usable)
  expect_false(target$verified)
  target_v0_verification <-
    target$compatibility$views$target_v0$verification
  expect_false(target_v0_verification$positive_requested_groups_represented)
  expect_gt(target_v0_verification$requested_group_masses[["inside_mass"]], 0)
  expect_identical(
    target_v0_verification$verification_group_masses[["inside_mass"]], 0
  )
  expect_match(target$message, "positive-mass group")

  equal_tail <- .dp_target_K(
    J = 20,
    K_interval = list(
      lower = 10L, upper = 19L, type = "equal_tail",
      coverage = 5 * smallest_positive, family = "maxent"
    )
  )
  expect_identical(equal_tail$status, "approximate")
  expect_false(equal_tail$usable)
  expect_false(equal_tail$verified)
  equal_tail_v0_verification <-
    equal_tail$compatibility$views$target_v0$verification
  expect_false(equal_tail_v0_verification$
    positive_requested_groups_representable)
  expect_identical(
    equal_tail_v0_verification$selected_group_masses[["inside_mass"]], 0
  )
  expect_gt(
    equal_tail_v0_verification$verification_group_masses[["inside_mass"]], 0
  )

  coarse_subnormal <- .dp_target_K(
    J = 10,
    K_interval = list(
      lower = 2L, upper = 3L, type = "equal_tail",
      coverage = 3 * smallest_positive, family = "maxent"
    )
  )
  expect_identical(coarse_subnormal$status, "approximate")
  coarse_v0_verification <-
    coarse_subnormal$compatibility$views$target_v0$verification
  expect_false(coarse_v0_verification$
    positive_requested_groups_representable)
  requested_inside <-
    coarse_v0_verification$requested_group_masses[["inside_mass"]]
  selected_inside <-
    coarse_v0_verification$selected_group_masses[["inside_mass"]]
  expect_gt(abs(selected_inside / requested_inside - 1), 0.1)
})

test_that("non-reference interval families are explicitly deferred", {
  for (family in c("discrete_normal", "beta_binomial")) {
    expect_error(
      .dp_target_K(
        J = 20, mu_K = 6.5,
        K_interval = list(
          lower = 3L, upper = 10L, type = "central_mass",
          coverage = 0.8, family = family
        )
      ),
      "deferred", class = "dpprior_interval_family_deferred"
    )
  }
})

test_that("interval backcheck separates selected and verification PMFs", {
  interval <- list(
    lower = 3L, upper = 10L, type = "equal_tail",
    coverage = 0.8, family = "maxent"
  )
  target <- .dp_target_K(J = 20, K_interval = interval)
  check <- .dp_backcheck_K_interval(
    target$pmf, target$interval,
    verification_pmf = target$pmf,
    J = 20, source = "calibrated_prior"
  )
  expect_s3_class(check, "dpprior_K_interval_backcheck")
  expect_identical(check$status, "converged")
  expect_equal(check$selected$achieved$inside_mass, 0.8, tolerance = 1e-14)
  expect_equal(check$verification$stability_max_abs, 0, tolerance = 0)
  expect_identical(
    check$provenance$moment_matching_assumed_interval_match, FALSE
  )

  unverified <- .dp_backcheck_K_interval(
    target$pmf, target$interval, J = 20
  )
  expect_identical(unverified$status, "approximate")
  expect_false(unverified$usable)
  expect_false(unverified$verified)
  expect_false(unverified$verification$performed)
  expect_identical(
    unverified$verification$reason,
    "independent_verification_not_supplied"
  )

  shifted <- target$pmf
  shifted[1] <- shifted[1] + 0.02
  shifted[3] <- shifted[3] - 0.02
  miss <- .dp_backcheck_K_interval(
    shifted, target$interval, verification_pmf = shifted, J = 20
  )
  expect_identical(miss$status, "approximate")
  expect_false(miss$verified)
  expect_equal(miss$selected$residuals$inside_mass, -0.02,
               tolerance = 1e-12)

  # A structural zero at K=1 is not guessed to be a legacy K=0 entry when J
  # is omitted.
  zero_at_one <- .dp_backcheck_K_interval(
    c(0, 0.5, 0.5),
    list(lower = 2L, upper = 3L, type = "hard_bounds", family = "maxent"),
    verification_pmf = c(0, 0.5, 0.5)
  )
  expect_identical(zero_at_one$J, 3L)
  expect_identical(zero_at_one$status, "converged")

  legacy_backcheck_matrix <- matrix(c(0, 0.5, 0.5), nrow = 3L)
  expect_error(
    .dp_backcheck_K_interval(
      legacy_backcheck_matrix,
      list(lower = 1L, upper = 2L, type = "hard_bounds",
           family = "maxent"),
      verification_pmf = c(0.5, 0.5), J = 2
    ),
    class = "dpprior_pmf_error"
  )

  forged_hard <- list(
    lower = 1L, upper = 5L, type = "hard_bounds", coverage = 999,
    family = 666, support = c(-1, 999), endpoints = "exclusive"
  )
  expect_error(
    .dp_backcheck_K_interval(
      rep(0.2, 5), forged_hard,
      verification_pmf = rep(0.2, 5), J = 5
    ),
    class = "dpprior_interval_backcheck_error"
  )
  forged_equal_tail <- list(
    lower = 2L, upper = 4L, type = "equal_tail", coverage = 2,
    family = "maxent", support = c(lower = 1L, upper = 5L),
    endpoints = "inclusive", mu_K = NULL
  )
  expect_error(
    .dp_backcheck_K_interval(
      rep(0.2, 5), forged_equal_tail,
      verification_pmf = rep(0.2, 5), J = 5
    ),
    class = "dpprior_interval_probability_error"
  )
})

test_that("interval type and family require plain character scalars", {
  base <- list(
    lower = 3L, upper = 10L, type = "central_mass",
    coverage = 0.8, family = "maxent"
  )
  malformed <- list(
    matrix("central_mass", nrow = 1L),
    array("central_mass", dim = c(1L, 1L, 1L)),
    structure("central_mass", class = "dpprior_test_character")
  )
  for (value in malformed) {
    interval <- base
    interval$type <- value
    expect_error(
      .dp_target_K(20, mu_K = 6.5, K_interval = interval),
      class = "dpprior_interval_structure_error"
    )
  }
  interval <- base
  interval$family <- matrix("maxent", nrow = 1L)
  expect_error(
    .dp_target_K(20, mu_K = 6.5, K_interval = interval),
    class = "dpprior_interval_structure_error"
  )

  interval <- base
  interval$coverage <- matrix(0.8, nrow = 1L)
  expect_error(
    .dp_target_K(20, mu_K = 6.5, K_interval = interval),
    class = "dpprior_type_error"
  )

  target <- .dp_target_K(20, mu_K = 6.5, K_interval = base)
  expect_error(
    .dp_backcheck_K_interval(
      target$pmf, target$interval, verification_pmf = target$pmf,
      J = 20, source = matrix("calibrated_prior", nrow = 1L)
    ),
    class = "dpprior_type_error"
  )
})

test_that("public target routes publish the frozen producer authority spine", {
  direct <- DPprior_target_K(J = 20, mu_K = 6, var_K = 5)
  expect_identical(names(direct$request), c("J", "mu_K", "var_K"))
  expect_identical(
    names(direct$normalized), c("J", "mu_K", "var_K", "interval", "pmf")
  )
  expect_identical(direct$used, direct$normalized)
  expect_identical(
    direct$derivation$request_to_normalized$rule,
    "canonicalize_direct_moments"
  )
  expect_identical(
    direct$derivation$request_to_normalized$before, direct$request
  )
  expect_identical(
    direct$derivation$request_to_normalized$after, direct$normalized
  )
  expect_null(direct$derivation$normalized_to_used)
  expect_null(direct$interval)
  expect_null(direct$family)
  expect_null(direct$pmf)
  expect_null(direct$computation$orders$M_selected)
  expect_null(direct$verification$selected_snapshot$M)
  expect_null(direct$verification$verifier_snapshot$M)
  expect_identical(
    direct$verification$selected_snapshot$achieved$implied, direct$implied
  )
  expect_identical(
    direct$verification$verifier_snapshot$achieved$implied, direct$implied
  )
  expect_false(direct$provenance$migration$lossless)
  expect_true("source_commit_not_embedded" %in%
                direct$provenance$migration$missing_evidence)
  expect_true(.dpprior_validate_target_v1(direct, collect = TRUE)$valid)

  confidence <- DPprior_target_K(
    J = 20, mu_K = 6, confidence = "medium"
  )
  expect_identical(
    names(confidence$request), c("J", "mean", "confidence")
  )
  expect_identical(
    confidence$derivation$request_to_normalized$rule,
    "canonicalize_confidence_target"
  )
  expect_identical(
    confidence$derivation$normalized_to_used$rule,
    "derive_variance_from_confidence_vif"
  )
  expect_identical(
    confidence$derivation$normalized_to_used$evidence,
    list(
      confidence = "medium", vif = 2.5,
      formula = "variance = vif * (mean - 1)"
    )
  )
  expect_identical(confidence$used$variance, 12.5)
  expect_null(confidence$family)

  cv <- DPprior_target_K(J = 20, mu_K = 6, cv_K = 0.25)
  expect_identical(names(cv$request), c("J", "mean", "cv"))
  expect_identical(
    cv$derivation$request_to_normalized$rule, "canonicalize_cv_target"
  )
  expect_identical(
    cv$derivation$normalized_to_used$rule, "derive_variance_from_cv"
  )
  expect_identical(
    cv$derivation$normalized_to_used$evidence$formula,
    "variance = (cv * mean)^2"
  )
  expect_identical(cv$used$variance, 2.25)
  expect_null(cv$family)

  pmf <- c(0.1, 0.2, 0.3, 0.4)
  strict <- DPprior_target_K(J = 4, target_pmf = pmf)
  structural <- DPprior_target_K(J = 4, target_pmf = c(0, pmf))
  expect_identical(strict$request, list(J = 4L, pmf = pmf))
  expect_identical(strict$normalized$pmf, strict$pmf)
  expect_identical(strict$used$pmf, strict$pmf)
  expect_identical(
    strict$derivation$request_to_normalized$rule, "validate_strict_pmf"
  )
  expect_identical(
    structural$derivation$request_to_normalized$rule,
    "drop_structural_k0_zero"
  )
  expect_identical(structural$request$pmf, c(0, pmf))
  expect_identical(structural$normalized$pmf, pmf)
  expect_null(strict$family)
  expect_null(structural$family)
  expect_identical(
    strict$verification$verifier_snapshot$achieved$pmf, strict$pmf
  )
  expect_true(.dpprior_validate_target_v1(strict, collect = TRUE)$valid)
  expect_true(.dpprior_validate_target_v1(structural, collect = TRUE)$valid)
})

test_that("verified interval routes retain closed construction truth", {
  cases <- list(
    hard_bounds = list(
      mu_K = NULL,
      interval = list(
        lower = 3L, upper = 10L, type = "hard_bounds", family = "maxent"
      ),
      rule = "construct_maxent_hard_bounds_pmf"
    ),
    equal_tail = list(
      mu_K = NULL,
      interval = list(
        lower = 3L, upper = 10L, type = "equal_tail",
        coverage = 0.8, family = "maxent"
      ),
      rule = "construct_maxent_equal_tail_pmf"
    ),
    central_mass = list(
      mu_K = 6.5,
      interval = list(
        lower = 3L, upper = 10L, type = "central_mass",
        coverage = 0.8, family = "maxent"
      ),
      rule = "construct_maxent_central_mass_pmf"
    )
  )
  for (case in cases) {
    target <- DPprior_target_K(
      J = 20, mu_K = case$mu_K, K_interval = case$interval
    )
    expect_true(.dpprior_validate_target_v1(target, collect = TRUE)$valid)
    expect_identical(
      names(target$interval),
      c(
        "lower", "upper", "type", "coverage", "family", "mu_K",
        "support", "endpoints"
      )
    )
    expect_identical(
      names(target$request), c("J", "K_interval", "mu_K")
    )
    expect_identical(target$request$K_interval, target$interval)
    expect_identical(target$normalized$interval, target$interval)
    expect_null(target$normalized$pmf)
    expect_identical(target$used$pmf, target$pmf)
    expect_identical(
      target$derivation$request_to_normalized$before, target$request
    )
    expect_identical(
      target$derivation$request_to_normalized$after, target$normalized
    )
    expect_identical(
      target$derivation$normalized_to_used$rule, case$rule
    )
    expect_identical(
      target$derivation$normalized_to_used$before, target$normalized
    )
    expect_identical(
      target$derivation$normalized_to_used$after, target$used
    )
    expect_identical(
      names(target$tolerances),
      c("constraint", "pmf_l1", "moment_relative", "moment_scale_floor")
    )
    expect_identical(
      names(target$computation$used$controls),
      c(
        "constraint_tolerance", "root_tolerance", "max_iterations",
        "pmf_l1_tolerance", "moment_relative_tolerance",
        "moment_scale_floor"
      )
    )
    expect_identical(
      target$computation$request$controls,
      target$computation$used$controls
    )
    expect_identical(
      target$verification$method, "independent_target_reconstruction"
    )
    expect_identical(target$verification$reason, "verified")
    expect_identical(
      names(target$verification$components),
      c("target_reconstruction", "order_stability")
    )
    expect_identical(
      names(target$verification$invariants),
      c("support_identity", "authority_identity")
    )
    expect_identical(
      target$verification$selected_snapshot$source, "target_constructor"
    )
    expect_identical(
      target$verification$verifier_snapshot$source,
      "independent_target_reconstruction"
    )
    expect_identical(
      target$verification$selected_snapshot$achieved$pmf, target$pmf
    )
    expect_identical(
      target$verification$stability$source,
      "independent_target_reconstruction"
    )
    expect_true(target$verification$stability$passed)
    expect_null(target$verification$selected_snapshot$M)
    expect_null(target$verification$verifier_snapshot$M)
    expect_null(target$computation$orders$M_selected)
  }

  singleton <- DPprior_target_K(
    J = 20,
    K_interval = list(
      lower = 5L, upper = 5L, type = "hard_bounds", family = "maxent"
    )
  )
  expect_identical(singleton$request$mu_K, 5)
  expect_true(
    singleton$derivation$request_to_normalized$evidence$
      singleton_support_implied_mean
  )
  expect_null(
    singleton$compatibility$views$target_v0$request$K_interval$mu_K
  )
})

test_that("infeasible and carried failure routes have typed proof evidence", {
  infeasible <- DPprior_target_K(
    J = 10,
    K_interval = list(
      lower = 3L, upper = 10L, type = "equal_tail",
      coverage = 0.8, family = "maxent"
    )
  )
  expect_true(.dpprior_validate_target_v1(infeasible, collect = TRUE)$valid)
  expect_identical(infeasible$status, "infeasible")
  expect_identical(infeasible$normalized, infeasible$used)
  expect_null(infeasible$derivation$normalized_to_used)
  expect_null(infeasible$pmf)
  expect_null(infeasible$implied)
  expect_identical(
    infeasible$verification$method,
    "analytic_interval_infeasibility_certificate"
  )
  certificate <- infeasible$verification$settings$certificate
  expect_identical(
    names(certificate),
    c(
      "version", "method", "kind", "assumptions", "request", "J",
      "support", "interval", "family", "group_masses",
      "group_support_counts", "empty_groups", "requested_mean",
      "feasible_mean_lower", "feasible_mean_upper",
      "feasibility_tolerance", "side", "outside_distance", "certified",
      "source"
    )
  )
  expect_identical(certificate$request, infeasible$request)
  expect_identical(certificate$interval, infeasible$interval)
  expect_identical(certificate$family, infeasible$family)
  expect_identical(
    certificate$kind, "positive_mass_group_has_empty_support"
  )
  expect_identical(certificate$empty_groups, "right")
  expect_identical(
    names(infeasible$verification$components), "infeasibility_certificate"
  )
  expect_null(infeasible$verification$selected_snapshot)
  expect_null(infeasible$verification$verifier_snapshot)
  expect_identical(length(infeasible$attempts), 1L)
  expect_identical(infeasible$attempts[[1L]]$stage, "feasibility")
  expect_identical(
    infeasible$attempts[[1L]]$reason_code,
    "globally_infeasible_by_certificate"
  )
  expect_identical(
    names(infeasible$attempts[[1L]]$unavailable),
    c("candidate_parameters", "candidate_objective")
  )
  expect_silent(.dpprior_schema_validate_plain_record_value(
    infeasible$compatibility$views$target_v0,
    "test.target.compatibility.views.target_v0"
  ))
  condition <- tryCatch(
    .dp_target_K_stop_unusable(infeasible), error = identity
  )
  expect_s3_class(condition, "dpprior_interval_empty_tail_error")
  expect_identical(condition$code, "empty_interval_group")
  expect_identical(condition$result, infeasible)
  expect_null(condition$call)

  forged_compatibility <- infeasible
  forged_compatibility$compatibility$views$target_v0$condition$classes <- c(
    "forged_error", "error", "condition"
  )
  forged_compatibility$compatibility$views$target_v0$condition$fields$code <-
    "forged_code"
  forged_compatibility$compatibility$views$target_v0$condition$fields$message <-
    "forged message"
  expect_true(.dpprior_validate_target_v1(
    forged_compatibility, collect = TRUE
  )$valid)
  canonical_condition <- tryCatch(
    .dp_target_K_stop_unusable(forged_compatibility), error = identity
  )
  expect_s3_class(canonical_condition, "dpprior_interval_empty_tail_error")
  expect_false(inherits(canonical_condition, "forged_error"))
  expect_identical(canonical_condition$code, "empty_interval_group")
  expect_identical(canonical_condition$message, infeasible$message)
  expect_identical(canonical_condition$result, forged_compatibility)

  approximate <- DPprior_target_K(
    J = 20, mu_K = 6.5,
    K_interval = list(
      lower = 3L, upper = 10L, type = "central_mass",
      coverage = 0.8, family = "maxent"
    ),
    root_control = list(max_iterations = 1L)
  )
  expect_true(.dpprior_validate_target_v1(approximate, collect = TRUE)$valid)
  expect_identical(approximate$status, "approximate")
  expect_null(approximate$pmf)
  expect_null(approximate$implied)
  expect_null(approximate$derivation$normalized_to_used)
  expect_false(approximate$verification$passed)
  expect_null(approximate$verification$selected_snapshot)
  expect_true(is.numeric(
    approximate$compatibility$views$target_v0$pmf
  ))
  expect_identical(approximate$attempts[[1L]]$exit_code, 1L)
})

test_that("interval truth controls are closed before construction", {
  interval <- list(
    lower = 3L, upper = 10L, type = "equal_tail",
    coverage = 0.8, family = "maxent"
  )
  explicit <- DPprior_target_K(
    J = 20, K_interval = interval,
    root_control = list(root_tol = 1e-13, max_iterations = 77L)
  )
  expect_identical(
    explicit$computation$used$controls$root_tolerance, 1e-13
  )
  expect_identical(
    explicit$computation$used$controls$max_iterations, 77L
  )
  expect_identical(
    explicit$derivation$normalized_to_used$evidence$root_tolerance,
    1e-13
  )
  expect_error(
    DPprior_target_K(J = 20, K_interval = interval, tolerance = 1e-8),
    class = "dpprior_bounds_error"
  )
  expect_error(
    DPprior_target_K(
      J = 20, K_interval = interval,
      root_control = list(root_tol = 1e-11)
    ),
    class = "dpprior_bounds_error"
  )
  expect_error(
    DPprior_target_K(
      J = 20, K_interval = interval,
      root_control = list(max_iterations = 1001L)
    ),
    class = "dpprior_bounds_error"
  )
})

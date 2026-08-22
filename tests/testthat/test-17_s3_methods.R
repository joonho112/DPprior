# Canonical S3 consumer tests -------------------------------------------------

.s3_17_cache <- new.env(parent = emptyenv())

.s3_17_cached <- function(key, producer) {
  if (!exists(key, envir = .s3_17_cache, inherits = FALSE)) {
    assign(key, producer(), envir = .s3_17_cache)
  }
  unserialize(serialize(
    get(key, envir = .s3_17_cache, inherits = FALSE), NULL, version = 3L
  ))
}

.s3_17_a1 <- function() .s3_17_cached("a1", function() {
  DPprior_fit(
    20L, mu_K = 4, confidence = "medium", method = "A1",
    check_diagnostics = FALSE
  )
})

.s3_17_a2 <- function() .s3_17_cached("a2", function() {
  DPprior_fit(
    20L, 4, 8, method = "A2-MN", M = 80L,
    check_diagnostics = FALSE
  )
})

.s3_17_hard <- function() .s3_17_cached("hard", function() {
  fit <- DPprior_fit(
    50L, 3, 10, method = "A2-MN", M = 80L,
    check_diagnostics = FALSE
  )
  DPprior_dual_hard(
    fit,
    list(
      metric = "wsb_tail", threshold = 0.5,
      relation = "<=", bound = 0.25
    ),
    M = 20L, M_verify = 60L,
    control = list(
      scan_points = 9L, scan_keep = 5L,
      profile_starts = 2L, maxit = 100L
    )
  )
})

.s3_17_soft <- function() .s3_17_cached("soft", function() {
  DPprior_dual_soft(
    .s3_17_a2(),
    list(
      metric = "wsb_tail", relation = "target",
      threshold = 0.5, value = 0.3
    ),
    lambda = 1, M_fit = 80L, M_verify = 160L
  )
})

.s3_17_legacy <- function() .s3_17_cached("legacy", function() {
  suppressWarnings(DPprior_dual(
    .s3_17_a2(),
    list(prob = list(threshold = 0.5, value = 0.3)),
    lambda = 1, M = 80L
  ))
})

.s3_17_no_candidate <- function() .s3_17_cached("no-candidate", function() {
  condition <- tryCatch(
    DPprior_dual_hard(
      .s3_17_a2(),
      list(
        metric = "wsb_tail", threshold = 0.5,
        relation = "<=", bound = 1e-8
      ),
      M = 20L, M_verify = 60L, log_bounds = c(-2, 2),
      control = list(
        scan_points = 7L, scan_keep = 5L,
        profile_starts = 2L, maxit = 40L
      )
    ),
    error = identity
  )
  expect_s3_class(condition, "dpprior_dual_infeasible")
  condition$result
})

.s3_17_expect_canonical <- function(fit, mode) {
  expect_s3_class(fit, "DPprior_fit")
  expect_s3_class(fit, "dpprior_result")
  expect_identical(unclass(fit)[["mode", exact = TRUE]], mode)
  expect_silent(.dpprior_validate_result_v1(fit))
}

.s3_17_with_null_device <- function(expr) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
}

test_that("fixtures exercise the five canonical fit modes", {
  fits <- list(
    a1 = .s3_17_a1(), a2 = .s3_17_a2(), hard = .s3_17_hard(),
    soft = .s3_17_soft(), legacy = .s3_17_legacy()
  )
  modes <- c("a1_proxy", "a2_moment", "dual_hard", "dual_soft", "dual_legacy")
  Map(.s3_17_expect_canonical, fits, modes)
})

test_that("print, summary, and data-frame methods expose canonical authority", {
  for (fit in list(
    .s3_17_a1(), .s3_17_a2(), .s3_17_hard(),
    .s3_17_soft(), .s3_17_legacy()
  )) {
    raw <- unclass(fit)
    printed <- capture.output(visible <- withVisible(print(fit)))
    expect_false(visible$visible)
    expect_identical(visible$value, fit)
    expect_true(any(grepl("Schema: dpprior.result/1", printed, fixed = TRUE)))
    expect_true(any(grepl(raw[["mode", exact = TRUE]], printed, fixed = TRUE)))
    expect_true(any(grepl(raw[["status", exact = TRUE]], printed, fixed = TRUE)))

    summary <- summary(fit, print_output = FALSE)
    expect_s3_class(summary, "summary.DPprior_fit")
    expect_identical(summary$schema, "dpprior.result/1")
    expect_identical(summary$mode, raw[["mode", exact = TRUE]])
    expect_identical(summary$status, raw[["status", exact = TRUE]])
    expect_identical(summary$verified, raw[["verified", exact = TRUE]])
    expect_identical(
      summary$gamma_prior,
      unclass(raw[["parameters", exact = TRUE]])[c("a", "b")]
    )

    frame <- as.data.frame(fit)
    expect_identical(nrow(frame), 1L)
    expect_identical(frame$schema, "dpprior.result/1")
    expect_identical(frame$mode, raw[["mode", exact = TRUE]])
    expect_identical(frame$status, raw[["status", exact = TRUE]])
    expect_identical(frame$verified, raw[["verified", exact = TRUE]])
  }
})

test_that("summary print is status-aware", {
  summary <- summary(.s3_17_a1(), print_output = FALSE)
  output <- capture.output(visible <- withVisible(print(summary)))
  expect_false(visible$visible)
  expect_true(any(grepl("DPprior Prior Elicitation Summary", output)))
  expect_true(any(grepl("approximate", output)))
  expect_true(any(grepl("verified: no", output, fixed = TRUE)))
})

test_that("the S3 plot method passes the canonical object without an adapter", {
  fit <- .s3_17_a2()
  received <- NULL
  testthat::local_mocked_bindings(
    plot_alpha_prior = function(fit, ...) {
      received <<- fit
      structure(list(), class = "captured_plot")
    },
    .package = "DPprior"
  )
  result <- plot(fit, type = "alpha", engine = "base", show = FALSE)
  expect_s3_class(result, "captured_plot")
  expect_identical(received, fit)
})

test_that("auto dispatch selects dual comparison only for current hard and soft", {
  received <- character()
  testthat::local_mocked_bindings(
    plot_prior_dashboard = function(fit, ...) {
      received <<- c(received, paste0("single:", unclass(fit)[["mode"]]))
      invisible(NULL)
    },
    plot_dual_comparison = function(fit_dual, ...) {
      received <<- c(received, paste0("dual:", unclass(fit_dual)[["mode"]]))
      invisible(NULL)
    },
    .package = "DPprior"
  )
  invisible(plot(.s3_17_a1(), type = "auto", engine = "base", show = FALSE))
  invisible(plot(.s3_17_a2(), type = "auto", engine = "base", show = FALSE))
  invisible(plot(.s3_17_hard(), type = "auto", engine = "base", show = FALSE))
  invisible(plot(.s3_17_soft(), type = "auto", engine = "base", show = FALSE))
  invisible(plot(.s3_17_legacy(), type = "auto", engine = "base", show = FALSE))
  expect_identical(received, c(
    "single:a1_proxy", "single:a2_moment", "dual:dual_hard",
    "dual:dual_soft", "single:dual_legacy"
  ))
})

test_that("all S3 plot types preserve base show=FALSE", {
  for (fit in list(.s3_17_a1(), .s3_17_a2(), .s3_17_legacy())) {
    for (type in c("auto", "dashboard", "alpha", "K", "w1")) {
      expect_null(plot(fit, type = type, engine = "base", show = FALSE))
    }
  }
  for (fit in list(.s3_17_hard(), .s3_17_soft())) {
    for (type in c("auto", "dual", "comparison")) {
      expect_null(plot(fit, type = type, engine = "base", show = FALSE))
    }
  }
})

test_that("ggplot S3 routes work with canonical A1, A2, hard, soft, and legacy", {
  skip_if_not_installed("ggplot2")
  for (fit in list(.s3_17_a1(), .s3_17_a2(), .s3_17_legacy())) {
    expect_s3_class(plot(fit, type = "alpha", show = FALSE), "ggplot")
    expect_s3_class(plot(fit, type = "K", show = FALSE), "ggplot")
    expect_s3_class(plot(fit, type = "w1", show = FALSE), "ggplot")
  }
  for (fit in list(.s3_17_hard(), .s3_17_soft())) {
    result <- .s3_17_with_null_device(plot(fit, type = "auto", show = FALSE))
    expect_true(inherits(result, "gtable") || is.list(result))
  }
})

test_that("legacy comparison has typed unavailable lineage", {
  for (type in c("dual", "comparison")) {
    condition <- tryCatch(
      plot(.s3_17_legacy(), type = type, engine = "base", show = FALSE),
      error = identity
    )
    expect_s3_class(condition, "dpprior_s3_unavailable_error")
    expect_s3_class(condition, "dpprior_visualization_data_error")
    expect_identical(condition$code, "canonical_comparison_lineage_unavailable")
  }
})

test_that("non-dual comparison fails typed instead of falling back", {
  condition <- tryCatch(
    plot(.s3_17_a2(), type = "comparison", engine = "base", show = FALSE),
    error = identity
  )
  expect_s3_class(condition, "dpprior_visualization_input_error")
  expect_identical(condition$code, "visualization_not_dual_fit")
})

test_that("candidate-unavailable canonical objects fail before plotting", {
  fit <- .s3_17_no_candidate()
  for (type in c("auto", "dashboard", "alpha", "K", "w1")) {
    condition <- tryCatch(
      plot(fit, type = type, engine = "base", show = FALSE),
      error = identity
    )
    expect_s3_class(condition, "dpprior_s3_unavailable_error")
    expect_s3_class(condition, "dpprior_visualization_data_error")
    expect_identical(condition$code, "canonical_fit_candidate_unavailable")
  }
})

test_that("flat and unknown DPprior_fit records are rejected by the schema gate", {
  flat <- structure(
    list(a = 2, b = 1, J = 20L, status = "converged"),
    class = "DPprior_fit"
  )
  for (operation in list(
    function() print(flat),
    function() summary(flat, print_output = FALSE),
    function() as.data.frame(flat),
    function() plot(flat, engine = "base", show = FALSE)
  )) {
    condition <- tryCatch(operation(), error = identity)
    expect_s3_class(condition, "dpprior_schema_unsupported_error")
    expect_s3_class(condition, "dpprior_serialization_error")
  }
})

test_that("dual detection is canonical mode detection", {
  expect_false(.dpprior_is_dual(.s3_17_a1()))
  expect_false(.dpprior_is_dual(.s3_17_a2()))
  expect_true(.dpprior_is_dual(.s3_17_hard()))
  expect_true(.dpprior_is_dual(.s3_17_soft()))
  expect_true(.dpprior_is_dual(.s3_17_legacy()))
  expect_false(.dpprior_is_dual(list(dual_anchor = list(init = list()))))
})

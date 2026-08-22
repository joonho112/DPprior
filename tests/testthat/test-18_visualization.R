# Canonical visualization consumer tests ------------------------------------

.viz18_cache <- new.env(parent = emptyenv())

.viz18_cached <- function(key, producer) {
  if (!exists(key, envir = .viz18_cache, inherits = FALSE)) {
    assign(key, producer(), envir = .viz18_cache)
  }
  unserialize(serialize(
    get(key, envir = .viz18_cache, inherits = FALSE), NULL, version = 3L
  ))
}

.viz18_a1 <- function() .viz18_cached("a1", function() {
  DPprior_fit(
    20L, mu_K = 4, confidence = "medium", method = "A1",
    check_diagnostics = FALSE
  )
})

.viz18_a2 <- function() .viz18_cached("a2", function() {
  DPprior_fit(
    20L, 4, 8, method = "A2-MN", M = 80L,
    check_diagnostics = FALSE
  )
})

.viz18_hard <- function() .viz18_cached("hard", function() {
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

.viz18_soft <- function() .viz18_cached("soft", function() {
  DPprior_dual_soft(
    .viz18_a2(),
    list(
      metric = "wsb_tail", relation = "target",
      threshold = 0.5, value = 0.3
    ),
    lambda = 1, M_fit = 80L, M_verify = 160L
  )
})

.viz18_legacy <- function() .viz18_cached("legacy", function() {
  suppressWarnings(DPprior_dual(
    .viz18_a2(),
    list(prob = list(threshold = 0.5, value = 0.3)),
    lambda = 1, M = 80L
  ))
})

.viz18_no_candidate <- function() .viz18_cached("no-candidate", function() {
  fit <- DPprior_fit(
    20L, 5, 8, method = "A2-MN", M = 80L,
    check_diagnostics = FALSE
  )
  condition <- tryCatch(
    DPprior_dual_hard(
      fit,
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

.viz18_target <- function() {
  list(
    metric = "wsb_tail", relation = "target",
    threshold = 0.5, value = 0.3
  )
}

.viz18_curve <- function(kind = c("complete", "mixed", "failed")) {
  kind <- match.arg(kind)
  .viz18_cached(paste0("curve-", kind), function() {
    fit_fun <- switch(
      kind,
      complete = NULL,
      mixed = function(fit, target, lambda, ...) {
        if (identical(lambda, 0.7)) stop("injected curve gap")
        DPprior_dual_soft(fit, target, lambda, ...)
      },
      failed = function(...) stop("injected all-point failure")
    )
    control <- if (is.null(fit_fun)) list() else list(.fit_fun = fit_fun)
    compute_tradeoff_curve(
      20L, list(mu_K = 4, var_K = 8),
      lambda_seq = if (identical(kind, "failed")) 0.7 else c(0.4, 0.7, 1),
      M = 80L, M_verify = 160L, target = .viz18_target(),
      control = control
    )
  })
}

.viz18_null_device <- function(expr) {
  grDevices::pdf(file = nullfile())
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
}

.viz18_capture <- function(expr) {
  tryCatch(eval.parent(substitute(expr)), error = identity)
}

test_that("colors and theme retain the public visualization contract", {
  colors <- DPprior_colors()
  expect_identical(
    names(colors),
    c("primary", "secondary", "accent", "ink", "warning", "shade",
      "k_only", "dual")
  )
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", unlist(colors))))
  skip_if_not_installed("ggplot2")
  expect_s3_class(theme_DPprior(), "theme")
})

test_that("one canonical fit view serves A1, A2, hard, soft, and legacy", {
  fits <- list(
    .viz18_a1(), .viz18_a2(), .viz18_hard(),
    .viz18_soft(), .viz18_legacy()
  )
  modes <- c("a1_proxy", "a2_moment", "dual_hard", "dual_soft", "dual_legacy")
  for (index in seq_along(fits)) {
    view <- .dpprior_visualization_fit_view(fits[[index]])
    raw <- unclass(fits[[index]])
    expect_identical(view$core$mode, modes[[index]])
    expect_identical(
      view$parameters,
      .dpprior_visualization_plain_record(raw[["parameters", exact = TRUE]])
    )
    expect_identical(
      view$target$K,
      .dpprior_visualization_plain_record(
        raw[["target", exact = TRUE]][["K", exact = TRUE]]
      )
    )
    expect_identical(
      view$achieved$K,
      .dpprior_visualization_plain_record(
        raw[["achieved", exact = TRUE]][["K", exact = TRUE]]
      )
    )
  }
})

test_that("single-fit ggplot views are canonical for A1, A2, and legacy", {
  skip_if_not_installed("ggplot2")
  for (fit in list(.viz18_a1(), .viz18_a2(), .viz18_legacy())) {
    expect_s3_class(plot_alpha_prior(fit, show = FALSE), "ggplot")
    expect_s3_class(plot_K_prior(fit, show = FALSE), "ggplot")
    expect_s3_class(plot_w1_prior(fit, show = FALSE), "ggplot")
    dashboard <- .viz18_null_device(plot_prior_dashboard(fit, show = FALSE))
    expect_true(inherits(dashboard, "gtable") || is.list(dashboard))
  }
})

test_that("single-fit captions name status, verification, and estimand", {
  skip_if_not_installed("ggplot2")
  a1_alpha <- plot_alpha_prior(.viz18_a1(), show = FALSE)
  expect_match(a1_alpha$labels$caption, "mode=a1_proxy")
  expect_match(a1_alpha$labels$caption, "status=approximate")
  expect_match(a1_alpha$labels$caption, "verified=no")
  expect_match(a1_alpha$labels$caption, "estimand=Gamma concentration")

  a2_K <- plot_K_prior(.viz18_a2(), show = FALSE)
  expect_match(a2_K$labels$caption, "target estimand=")
  expect_match(a2_K$labels$caption, "achieved estimand=")

  legacy_w <- plot_w1_prior(.viz18_legacy(), show = FALSE)
  expect_match(legacy_w$labels$caption, "mode=dual_legacy")
  expect_match(legacy_w$labels$caption, "deprecated legacy dual fit")
})

test_that("base show=FALSE validates but creates no plotting device", {
  before <- grDevices::dev.cur()
  for (fit in list(.viz18_a1(), .viz18_a2(), .viz18_legacy())) {
    expect_null(plot_alpha_prior(fit, engine = "base", show = FALSE))
    expect_null(plot_K_prior(fit, engine = "base", show = FALSE))
    expect_null(plot_w1_prior(fit, engine = "base", show = FALSE))
    expect_null(plot_prior_dashboard(fit, engine = "base", show = FALSE))
  }
  expect_identical(grDevices::dev.cur(), before)
})

test_that("direct parameter plotting remains separate from fit science", {
  skip_if_not_installed("ggplot2")
  expect_s3_class(plot_alpha_prior(a = 2, b = 1.5, show = FALSE), "ggplot")
  expect_s3_class(
    plot_K_prior(J = 20L, a = 2, b = 1.5, show = FALSE), "ggplot"
  )
  expect_s3_class(plot_w1_prior(a = 2, b = 1.5, show = FALSE), "ggplot")
  expect_null(plot_alpha_prior(
    a = 2, b = 1.5, engine = "base", show = FALSE
  ))
})

test_that("canonical summary uses nested target and achieved records", {
  fit <- .viz18_a2()
  raw <- unclass(fit)
  summary <- .dpprior_compute_summary(fit)
  expect_identical(
    summary$target$mu_K,
    raw[["target", exact = TRUE]][["K", exact = TRUE]][[
      "implied", exact = TRUE
    ]][["mean", exact = TRUE]]
  )
  expect_identical(
    summary$achieved$mu_K,
    raw[["achieved", exact = TRUE]][["K", exact = TRUE]][[
      "mean", exact = TRUE
    ]]
  )
  expect_identical(
    summary$achieved$estimand,
    raw[["achieved", exact = TRUE]][["K", exact = TRUE]][[
      "estimand", exact = TRUE
    ]]
  )
})

test_that("compatibility alias tampering is rejected rather than plotted", {
  fit <- .viz18_a2()
  raw <- unclass(fit)
  raw[["a"]] <- raw[["parameters", exact = TRUE]][["a", exact = TRUE]] + 10
  class(raw) <- class(fit)
  condition <- .viz18_capture(plot_alpha_prior(raw, show = FALSE))
  expect_s3_class(condition, "dpprior_schema_error")
  expect_identical(condition$code, "alias_identity")
})

test_that("canonical fit access bypasses poisoned fit accessors", {
  local({
    calls <- 0L
    assign("$.DPprior_fit", function(...) {
      calls <<- calls + 1L
      stop("poison dollar")
    }, envir = environment())
    assign("[[.DPprior_fit", function(...) {
      calls <<- calls + 1L
      stop("poison bracket")
    }, envir = environment())
    expect_s3_class(plot_alpha_prior(.viz18_a2(), show = FALSE), "ggplot")
    expect_s3_class(plot_K_prior(.viz18_a2(), show = FALSE), "ggplot")
    expect_s3_class(plot_w1_prior(.viz18_a2(), show = FALSE), "ggplot")
    expect_identical(calls, 0L)
  })
})

test_that("nested canonical target classes are distilled before scientific reads", {
  fit <- .viz18_a2()
  calls <- 0L
  method_name <- "[[.dpprior_K_target"
  had_method <- exists(method_name, envir = .GlobalEnv, inherits = FALSE)
  old_method <- if (had_method) {
    get(method_name, envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (had_method) {
      assign(method_name, old_method, envir = .GlobalEnv)
    } else if (exists(method_name, envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = method_name, envir = .GlobalEnv)
    }
  }, add = TRUE)
  assign(method_name, function(...) {
    calls <<- calls + 1L
    stop("poison nested target accessor")
  }, envir = .GlobalEnv)

  expect_null(plot_K_prior(fit, engine = "base", show = FALSE))
  expect_null(plot(fit, type = "K", engine = "base", show = FALSE))
  expect_identical(calls, 0L)
})

test_that("hard and soft dual comparisons use authoritative input-fit lineage", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("gtable")
  for (fit in list(.viz18_hard(), .viz18_soft())) {
    expect_null(plot_dual_comparison(
      fit, engine = "base", show = FALSE
    ))
    result <- .viz18_null_device(plot_dual_comparison(fit, show = FALSE))
    expect_true(inherits(result, "gtable") || is.list(result))
    dashboard <- .viz18_null_device(plot_dual_dashboard(fit, show = FALSE))
    expect_true(inherits(dashboard, "gtable") || is.list(dashboard))
  }
})

test_that("dual summary is hard-constraint or soft-tradeoff specific", {
  skip_if_not_installed("ggplot2")
  for (case in list(
    list(fit = .viz18_hard(), required = "Constraint", absent = "Lambda"),
    list(fit = .viz18_soft(), required = "Lambda", absent = "Constraint")
  )) {
    dual <- .dpprior_visualization_dual_view(case$fit)
    baseline <- dual$baseline$parameters
    current <- dual$current$parameters
    J <- dual$current$core$J
    table <- .dpprior_dual_comparison_table(
      dual,
      .dpprior_compute_summary(a = baseline$a, b = baseline$b, J = J),
      .dpprior_compute_summary(a = current$a, b = current$b, J = J)
    )
    expect_true(case$required %in% table$data$Metric)
    expect_false(case$absent %in% table$data$Metric)
  }
})

test_that("explicit comparison baseline must exactly match provenance.input_fit", {
  expect_null(plot_dual_comparison(
    .viz18_soft(), .viz18_a2(), engine = "base", show = FALSE
  ))
  mismatch <- DPprior_fit(
    20L, 4.2, 8, method = "A2-MN", M = 80L,
    check_diagnostics = FALSE
  )
  condition <- .viz18_capture(plot_dual_comparison(
    .viz18_soft(), mismatch, engine = "base", show = FALSE
  ))
  expect_s3_class(condition, "dpprior_visualization_input_error")
  expect_identical(condition$code, "comparison_input_fit_mismatch")
  expect_identical(condition$mismatch, "canonical_input_fit_reference")

  different_verifier <- DPprior_a2_newton(
    20L, 4, 8, M = 80L, M_verify = 240L, verbose = FALSE
  )
  condition <- .viz18_capture(plot_dual_comparison(
    .viz18_soft(), different_verifier, engine = "base", show = FALSE
  ))
  expect_s3_class(condition, "dpprior_visualization_input_error")
  expect_identical(condition$code, "comparison_input_fit_mismatch")
})

test_that("legacy comparison is typed unavailable and never uses compatibility", {
  for (fun in list(plot_dual_comparison, plot_dual_dashboard)) {
    condition <- .viz18_capture(fun(
      .viz18_legacy(), engine = "base", show = FALSE
    ))
    expect_s3_class(condition, "dpprior_s3_unavailable_error")
    expect_s3_class(condition, "dpprior_visualization_data_error")
    expect_identical(condition$code, "canonical_comparison_lineage_unavailable")
  }
})

test_that("candidate NULL is a typed visualization data error", {
  fit <- .viz18_no_candidate()
  for (fun in list(plot_alpha_prior, plot_K_prior, plot_w1_prior,
                   plot_prior_dashboard, plot_dual_comparison)) {
    condition <- .viz18_capture(fun(fit, engine = "base", show = FALSE))
    expect_s3_class(condition, "dpprior_s3_unavailable_error")
    expect_s3_class(condition, "dpprior_visualization_data_error")
    expect_identical(condition$code, "canonical_fit_candidate_unavailable")
    expect_identical(condition$status, "infeasible")
  }
})

test_that("non-dual input to a dual view fails typed", {
  condition <- .viz18_capture(plot_dual_comparison(
    .viz18_a2(), engine = "base", show = FALSE
  ))
  expect_s3_class(condition, "dpprior_visualization_input_error")
  expect_identical(condition$code, "visualization_not_dual_fit")
})

test_that("trade-off plotting requires producer-retained canonical evidence", {
  curve <- .viz18_curve("complete")
  expect_s3_class(curve, "dpprior_tradeoff_curve")
  expect_null(plot_tradeoff_curve(
    curve, engine = "base", show = FALSE
  ))
  skip_if_not_installed("ggplot2")
  for (metric in c("w1_prob_gt_50", "E_w1", "K_loss", "mu_K", "var_K")) {
    expect_s3_class(plot_tradeoff_curve(
      curve, metric = metric, show = FALSE
    ), "ggplot")
  }

  generic <- unclass(curve)
  class(generic) <- "data.frame"
  condition <- .viz18_capture(plot_tradeoff_curve(generic, show = FALSE))
  expect_s3_class(condition, "dpprior_visualization_input_error")
  expect_identical(condition$code, "tradeoff_curve_contract")

  for (field in c("K_loss", "w1_prob_gt_50", "E_w1")) {
    tampered <- curve
    tampered[[field]][[1L]] <- tampered[[field]][[1L]] + 1
    condition <- .viz18_capture(plot_tradeoff_curve(
      tampered, metric = field, show = FALSE
    ))
    expect_s3_class(condition, "dpprior_visualization_input_error")
    expect_identical(condition$code, "tradeoff_curve_science_mismatch")
  }

  tampered <- curve
  tampered$condition_code[[1L]] <- "forged_success_condition"
  condition <- .viz18_capture(plot_tradeoff_curve(tampered, show = FALSE))
  expect_s3_class(condition, "dpprior_visualization_input_error")
  expect_identical(condition$code, "tradeoff_curve_science_mismatch")
})

test_that("every public curve decision column is bound to retained evidence", {
  curve <- .viz18_curve("complete")
  mutators <- list(
    lambda = function(x) {
      x$lambda[[1L]] <- 0.41
      x
    },
    metric = function(x) {
      x$metric[[1L]] <- "wsb_mean"
      x
    },
    relation = function(x) {
      x$relation[[1L]] <- "at_most"
      x
    },
    target_value = function(x) {
      x$target_value[[1L]] <- 0.31
      x
    },
    attempt_count = function(x) {
      x$attempt_count[[1L]] <- x$attempt_count[[1L]] + 1L
      x
    },
    selected_method = function(x) {
      x$selected_method[[1L]] <- "forged-method"
      x
    },
    warm_start_from = function(x) {
      x$warm_start_from[[1L]] <- "forged warm source"
      x
    }
  )
  for (field in names(mutators)) {
    condition <- .viz18_capture(plot_tradeoff_curve(
      mutators[[field]](curve), engine = "base", show = FALSE
    ))
    expect_s3_class(condition, "dpprior_visualization_input_error")
    expect_identical(
      condition$code, "tradeoff_curve_science_mismatch"
    )
  }
})

test_that("curve identities, metadata, and condition presence are exact", {
  curve <- .viz18_curve("complete")

  coordinated_id <- curve
  fits <- attr(coordinated_id, "fits", exact = TRUE)
  conditions <- attr(coordinated_id, "conditions", exact = TRUE)
  forged_id <- paste0(coordinated_id$point_id[[1L]], "|forged")
  coordinated_id$point_id[[1L]] <- forged_id
  names(fits)[[1L]] <- forged_id
  names(conditions)[[1L]] <- forged_id
  attr(coordinated_id, "fits") <- fits
  attr(coordinated_id, "conditions") <- conditions
  condition <- .viz18_capture(plot_tradeoff_curve(
    coordinated_id, engine = "base", show = FALSE
  ))
  expect_identical(condition$code, "tradeoff_curve_science_mismatch")

  missing_metadata <- curve
  attr(missing_metadata, "metadata") <- list()
  condition <- .viz18_capture(plot_tradeoff_curve(
    missing_metadata, engine = "base", show = FALSE
  ))
  expect_identical(condition$code, "tradeoff_curve_retained_evidence")

  metadata_mutators <- list(
    evaluation_order = function(metadata) {
      metadata$evaluation_order <- rev(metadata$evaluation_order)
      metadata
    },
    target_K = function(metadata) {
      metadata$target_K <- list(mu_K = 99, var_K = 99)
      metadata
    },
    target_weight = function(metadata) {
      metadata$target_weight$value <- 0.99
      metadata
    },
    fixed_scales = function(metadata) {
      metadata$fixed_scales$K_mean <- 99
      metadata
    }
  )
  for (field in names(metadata_mutators)) {
    tampered <- curve
    attr(tampered, "metadata") <- metadata_mutators[[field]](
      attr(tampered, "metadata", exact = TRUE)
    )
    condition <- .viz18_capture(plot_tradeoff_curve(
      tampered, engine = "base", show = FALSE
    ))
    expect_identical(
      condition$code, "tradeoff_curve_science_mismatch", info = field
    )
  }

  fabricated <- curve
  fits <- attr(fabricated, "fits", exact = TRUE)
  conditions <- attr(fabricated, "conditions", exact = TRUE)
  conditions[[1L]] <- .dpprior_new_condition(
    "fabricated condition on successful fit",
    c(
      "dpprior_dual_soft_backend_contract_error",
      "dpprior_dual_soft_error", "dpprior_calibration_error",
      "dpprior_error", "error"
    ),
    code = "dual_soft_backend_contract", result = fits[[1L]]
  )
  fields <- .dpprior_soft_condition_fields(conditions[[1L]])
  fabricated$outcome[[1L]] <- "condition_retained"
  fabricated$condition_class[[1L]] <- fields$class
  fabricated$condition_code[[1L]] <- fields$code
  fabricated$condition_message[[1L]] <- fields$message
  attr(fabricated, "conditions") <- conditions
  condition <- .viz18_capture(plot_tradeoff_curve(
    fabricated, engine = "base", show = FALSE
  ))
  expect_identical(condition$code, "tradeoff_curve_science_mismatch")

  swapped <- curve
  fits <- attr(swapped, "fits", exact = TRUE)
  fits[[1L]] <- fits[[2L]]
  attr(swapped, "fits") <- fits
  copied_fields <- c(
    "mode", "status", "usable", "verified", "converged", "outcome",
    "condition_class", "condition_code", "condition_message", "a", "b",
    "mu_K", "var_K", "achieved_weight", "target_residual", "K_loss",
    "weight_loss", "total_loss", "w_loss", "w1_prob_gt_50", "E_w1"
  )
  for (field in copied_fields) {
    swapped[[field]][[1L]] <- swapped[[field]][[2L]]
  }
  condition <- .viz18_capture(plot_tradeoff_curve(
    swapped, engine = "base", show = FALSE
  ))
  expect_identical(condition$code, "tradeoff_curve_science_mismatch")
})

test_that("classed curve leaves are rejected before hostile dispatch", {
  curve <- .viz18_curve("complete")
  calls <- 0L
  specifications <- list(
    numeric_metric = list(field = "K_loss", class = "evil_curve_numeric"),
    numeric_lambda = list(field = "lambda", class = "evil_curve_numeric"),
    logical = list(field = "usable", class = "evil_curve_logical"),
    character = list(field = "metric", class = "evil_curve_character"),
    point_id = list(field = "point_id", class = "evil_curve_character")
  )
  cases <- lapply(names(specifications), function(name) {
    specification <- specifications[[name]]
    tampered <- unclass(curve)
    column <- tampered[[specification$field, exact = TRUE]]
    if (identical(name, "numeric_metric")) column[[1L]] <- 999
    attr(column, "class") <- specification$class
    tampered[[specification$field]] <- column
    attr(tampered, "class") <- c("dpprior_tradeoff_curve", "data.frame")
    tampered
  })
  names(cases) <- names(specifications)

  metadata_character <- unclass(curve)
  metadata <- attr(metadata_character, "metadata", exact = TRUE)
  leaf <- metadata$target_weight$raw_input$metric
  attr(leaf, "class") <- "evil_curve_character"
  metadata$target_weight$raw_input$metric <- leaf
  attr(metadata_character, "metadata") <- metadata
  attr(metadata_character, "class") <- c(
    "dpprior_tradeoff_curve", "data.frame"
  )
  cases$metadata_character <- metadata_character

  metadata_numeric <- unclass(curve)
  metadata <- attr(metadata_numeric, "metadata", exact = TRUE)
  leaf <- metadata$target_K$mu_K
  attr(leaf, "class") <- "evil_curve_numeric"
  metadata$target_K$mu_K <- leaf
  attr(metadata_numeric, "metadata") <- metadata
  attr(metadata_numeric, "class") <- c(
    "dpprior_tradeoff_curve", "data.frame"
  )
  cases$metadata_numeric <- metadata_numeric

  for (attribute in c("fits", "conditions", "metadata")) {
    tampered <- unclass(curve)
    retained <- attr(tampered, attribute, exact = TRUE)
    attr(retained, "class") <- "evil_curve_list"
    attr(tampered, attribute) <- retained
    attr(tampered, "class") <- c("dpprior_tradeoff_curve", "data.frame")
    cases[[paste0("classed_", attribute)]] <- tampered
  }

  condition_message <- unclass(.viz18_curve("failed"))
  conditions <- attr(condition_message, "conditions", exact = TRUE)
  condition_classes <- attr(conditions[[1L]], "class", exact = TRUE)
  condition_raw <- unclass(conditions[[1L]])
  leaf <- condition_raw$message
  attr(leaf, "class") <- "evil_curve_character"
  condition_raw$message <- leaf
  attr(condition_raw, "class") <- condition_classes
  conditions[[1L]] <- condition_raw
  attr(condition_message, "conditions") <- conditions
  attr(condition_message, "class") <- c(
    "dpprior_tradeoff_curve", "data.frame"
  )
  cases$condition_message <- condition_message

  methods <- list(
    "[[.evil_curve_numeric" = function(...) {
      calls <<- calls + 1L
      0
    },
    "is.finite.evil_curve_numeric" = function(...) {
      calls <<- calls + 1L
      TRUE
    },
    "Ops.evil_curve_numeric" = function(...) {
      calls <<- calls + 1L
      TRUE
    },
    "[[.evil_curve_logical" = function(...) {
      calls <<- calls + 1L
      TRUE
    },
    "[[.evil_curve_character" = function(...) {
      calls <<- calls + 1L
      "forged"
    },
    "length.evil_curve_character" = function(...) {
      calls <<- calls + 1L
      1L
    },
    "as.character.evil_curve_character" = function(...) {
      calls <<- calls + 1L
      "forged"
    },
    "xtfrm.evil_curve_character" = function(...) {
      calls <<- calls + 1L
      1
    },
    "dim.evil_curve_character" = function(...) {
      calls <<- calls + 1L
      NULL
    },
    "dim.evil_curve_numeric" = function(...) {
      calls <<- calls + 1L
      NULL
    },
    "dim.evil_curve_logical" = function(...) {
      calls <<- calls + 1L
      NULL
    },
    "length.evil_curve_numeric" = function(...) {
      calls <<- calls + 1L
      1L
    },
    "length.evil_curve_logical" = function(...) {
      calls <<- calls + 1L
      1L
    },
    "[[.evil_curve_list" = function(...) {
      calls <<- calls + 1L
      stop("hostile retained-list accessor")
    },
    "length.evil_curve_list" = function(...) {
      calls <<- calls + 1L
      1L
    },
    "names.evil_curve_list" = function(...) {
      calls <<- calls + 1L
      "forged"
    }
  )
  method_names <- names(methods)
  existed <- vapply(
    method_names, exists, logical(1L), envir = .GlobalEnv, inherits = FALSE
  )
  previous <- lapply(method_names, function(name) {
    if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
      get(name, envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
  })
  on.exit({
    for (index in seq_along(method_names)) {
      name <- method_names[[index]]
      if (existed[[index]]) {
        assign(name, previous[[index]], envir = .GlobalEnv)
      } else if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
        rm(list = name, envir = .GlobalEnv)
      }
    }
  }, add = TRUE)
  for (name in method_names) assign(name, methods[[name]], envir = .GlobalEnv)

  for (name in names(cases)) {
    condition <- .viz18_capture(plot_tradeoff_curve(
      cases[[name]], metric = "K_loss", engine = "base", show = FALSE
    ))
    expect_s3_class(condition, "dpprior_visualization_input_error")
  }
  expect_identical(calls, 0L)
})

test_that("mixed curve failures create visible gaps without imputation", {
  curve <- .viz18_curve("mixed")
  view <- .dpprior_visualization_curve_view(curve, "K_loss")
  expect_identical(view$requested_points, 3L)
  expect_identical(view$available_points, 2L)
  expect_identical(view$unavailable_points, 1L)
  expect_identical(view$line_data$segment, c(1L, 2L))
  expect_null(plot_tradeoff_curve(
    curve, metric = "K_loss", engine = "base", show = FALSE
  ))
  skip_if_not_installed("ggplot2")
  plot <- plot_tradeoff_curve(curve, metric = "K_loss", show = FALSE)
  expect_s3_class(plot, "ggplot")
  expect_match(plot$labels$subtitle, "gaps are not imputed")

  failure_index <- which(curve$lambda == 0.7)
  conditions <- attr(curve, "conditions", exact = TRUE)
  condition_raw <- unclass(conditions[[failure_index]])
  condition_raw[["result"]] <- attr(curve, "fits", exact = TRUE)[[1L]]
  class(condition_raw) <- class(conditions[[failure_index]])
  conditions[[failure_index]] <- condition_raw
  attr(curve, "conditions") <- conditions
  condition <- .viz18_capture(plot_tradeoff_curve(
    curve, metric = "K_loss", show = FALSE
  ))
  expect_s3_class(condition, "dpprior_visualization_input_error")
  expect_identical(condition$code, "tradeoff_curve_science_mismatch")
})

test_that("an all-unavailable curve metric fails typed", {
  condition <- .viz18_capture(plot_tradeoff_curve(
    .viz18_curve("failed"), metric = "K_loss",
    engine = "base", show = FALSE
  ))
  expect_s3_class(condition, "dpprior_s3_unavailable_error")
  expect_s3_class(condition, "dpprior_visualization_data_error")
  expect_identical(condition$code, "tradeoff_metric_unavailable")
})

test_that("trade-off dashboard uses the same validated no-imputation routes", {
  curve <- .viz18_curve("complete")
  expect_null(plot_tradeoff_dashboard(
    curve, engine = "base", show = FALSE
  ))
  tampered <- curve
  tampered$lambda[[1L]] <- 0.41
  condition <- .viz18_capture(plot_tradeoff_dashboard(
    tampered, engine = "base", show = FALSE
  ))
  expect_identical(condition$code, "tradeoff_curve_science_mismatch")
  skip_if_not_installed("ggplot2")
  result <- .viz18_null_device(plot_tradeoff_dashboard(curve, show = FALSE))
  expect_true(inherits(result, "gtable") || is.list(result))
})

test_that("threshold and fit contracts fail closed", {
  expect_error(
    plot_w1_prior(a = 2, b = 1.5, thresholds = c(0.9, 0.5), show = FALSE),
    class = "dpprior_visualization_input_error"
  )
  flat <- structure(list(a = 2, b = 1, J = 20L), class = "DPprior_fit")
  expect_error(
    plot_alpha_prior(flat, show = FALSE),
    class = "dpprior_schema_unsupported_error"
  )
})

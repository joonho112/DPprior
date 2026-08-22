.s3_v1_diagnostics <- function() {
  structure(list(
    J = 20L,
    a = 1.27223464009856,
    b = 0.58075377814199,
    alpha = list(
      mean = 2.19, sd = 1.94, cv = 0.887, median = 1.65,
      quantiles = c(`5%` = 0.19, `50%` = 1.65, `95%` = 6.03)
    ),
    K = list(
      mean = 5, var = 8, sd = sqrt(8), mode = 4L, median = 5L,
      quantiles = c(`5%` = 1L, `50%` = 5L, `95%` = 10L),
      pmf = rep(0.05, 20L)
    ),
    weights = list(
      mean = 0.42, median = 0.343,
      quantiles = c(`5%` = 0.024, `50%` = 0.343, `95%` = 0.996),
      prob_exceeds = c(`0.5` = 0.368, `0.9` = 0.13),
      dominance_risk = "moderate"
    ),
    coclustering = list(
      mean = 0.42, var = 0.073, sd = 0.269,
      interpretation = "Moderate prior co-clustering"
    ),
    warnings = character()
  ), class = "DPprior_diagnostics")
}


.s3_v1_a2_mn <- function() {
  structure(list(
    a = 1.27223464009856,
    b = 0.58075377814199,
    J = 20,
    target = list(mu_K = 5, var_K = 8, type = "moments"),
    method = "A2-MN",
    status = "success",
    converged = TRUE,
    iterations = 7L,
    termination = "residual",
    fit = list(
      mu_K = 5.00000000001361,
      var_K = 7.99999999987637,
      residual = 1.24381729905704e-10
    ),
    diagnostics = list(
      a0 = 4, b0 = 2.99573227355399, tol_F = 1e-8,
      tol_step = 1e-10, M = 80L, fallback_used = FALSE
    ),
    trace = data.frame(iteration = 1L, residual = 1.24381729905704e-10)
  ), class = "DPprior_fit")
}


.s3_v1_a1 <- function() {
  structure(list(
    a = 4,
    b = 2.99573227355399,
    J = 20,
    target = list(mu_K = 5, var_K = 8, type = "moments"),
    method = "A1",
    status = "success",
    scaling = list(cJ = log(20), rule = "legacy"),
    cJ = log(20),
    var_K_used = 8,
    converged = TRUE,
    iterations = 0L,
    fit = NULL,
    diagnostics = NULL,
    trace = NULL
  ), class = "DPprior_fit")
}


.s3_v1_a2_kl <- function() {
  pmf <- rep(1 / 20, 20L)
  support <- seq_along(pmf)
  mean <- sum(support * pmf)
  variance <- sum((support - mean)^2 * pmf)
  structure(list(
    a = 1.5,
    b = 0.7,
    J = 20L,
    target = list(
      type = "chisq", pmf = pmf, mu_K = mean, var_K = variance,
      df = 6.25, scale = 0.8,
      mu_K_discrete = mean, var_K_discrete = variance
    ),
    method = "A2-KL",
    status = "success",
    converged = TRUE,
    iterations = 8L,
    termination = "optim_converged",
    fit = list(mu_K = 5, var_K = 8, kl = 0.02, residual = 0.02),
    diagnostics = list(M = 80L, fallback_used = FALSE),
    trace = data.frame(evaluation = 1L, kl = 0.02)
  ), class = "DPprior_fit")
}


.s3_v1_wrapper <- function(method = "A2-MN",
                            include_diagnostics = FALSE) {
  raw <- unclass(.s3_v1_a2_mn())
  raw$target <- list(
    mu_K = 5, var_K = 8, var_K_used = 8,
    confidence = NULL, type = "moments"
  )
  raw$method <- method
  raw$solver_diagnostics <- raw$diagnostics
  raw$diagnostics <- NULL
  raw <- raw[c(
    "a", "b", "J", "target", "method", "status", "converged",
    "iterations", "termination", "fit", "solver_diagnostics", "trace"
  )]
  if (include_diagnostics) {
    raw$diagnostics <- .s3_v1_diagnostics()
  }
  structure(raw, class = "DPprior_fit")
}


.s3_v1_dual <- function() {
  raw <- unclass(.s3_v1_a2_mn())
  raw$target <- list(
    mu_K = 5, var_K = 8, var_K_used = 8,
    confidence = NULL, type = "moments"
  )
  raw$method <- "dual-anchor"
  raw$solver_diagnostics <- raw$diagnostics
  raw$diagnostics <- NULL
  raw$dual_anchor <- list(
    w1_target = list(prob = list(threshold = 0.5, value = 0.3)),
    lambda = 1,
    loss_type = "relative",
    w1_achieved = list(
      mean = 0.42, prob_gt_50 = 0.368, prob_gt_90 = 0.13
    ),
    K_loss = 0,
    init = list(a = raw$a, b = raw$b),
    note = "lambda = 1 returns the K-only solution"
  )
  raw <- raw[c(
    "a", "b", "J", "target", "method", "status", "converged",
    "iterations", "termination", "fit", "solver_diagnostics", "trace",
    "dual_anchor"
  )]
  structure(raw, class = "DPprior_fit")
}


.s3_upgrade <- function(x, allow_legacy = TRUE) {
  suppressWarnings(DPprior:::upgrade_DPprior_object(
    x, verify = TRUE, allow_legacy = allow_legacy
  ))
}


.s3_migrated_fits <- list(
  direct_a1 = .s3_upgrade(.s3_v1_a1()),
  direct_a2_mn = .s3_upgrade(.s3_v1_a2_mn()),
  direct_a2_kl = .s3_upgrade(.s3_v1_a2_kl()),
  wrapper_a1 = .s3_upgrade(.s3_v1_wrapper("A1")),
  wrapper_a2_mn = .s3_upgrade(.s3_v1_wrapper("A2-MN")),
  wrapper_with_diagnostics = .s3_upgrade(
    .s3_v1_wrapper("A2-MN", include_diagnostics = TRUE)
  ),
  dual = .s3_upgrade(.s3_v1_dual())
)


.s3_public_candidates <- c(
  direct_a1 = TRUE,
  direct_a2_mn = FALSE,
  direct_a2_kl = FALSE,
  wrapper_a1 = TRUE,
  wrapper_a2_mn = FALSE,
  wrapper_with_diagnostics = FALSE,
  dual = TRUE
)


.s3_catch <- function(expr) {
  tryCatch(expr, error = identity)
}


test_that("seven migrated fit shapes have truthful S3 consumer views", {
  for (name in names(.s3_migrated_fits)) {
    fit <- .s3_migrated_fits[[name]]
    output <- capture.output(visible <- withVisible(print(fit)))
    summary <- summary(fit, print_output = FALSE)
    compact <- as.data.frame(fit)
    candidate_available <- unname(.s3_public_candidates[[name]])

    expect_false(visible$visible, info = name)
    expect_identical(visible$value, fit, info = name)
    expect_true(any(grepl("Schema: dpprior.result/1", output, fixed = TRUE)),
                info = name)
    expect_true(any(grepl("Status: approximate", output, fixed = TRUE)),
                info = name)
    expect_true(any(grepl("verified: no", output, fixed = TRUE)), info = name)
    expect_false(any(grepl("SUCCESS", output, fixed = TRUE)), info = name)

    expect_s3_class(summary, "summary.DPprior_fit")
    expect_true(summary$canonical, info = name)
    expect_identical(summary$status, "approximate", info = name)
    expect_false(summary$verified, info = name)
    expect_identical(summary$candidate_available, candidate_available,
                     info = name)
    expect_match(summary$guidance, "refit", ignore.case = TRUE, info = name)

    expect_s3_class(compact, "data.frame")
    expect_identical(nrow(compact), 1L, info = name)
    expect_identical(compact$status, "approximate", info = name)
    expect_false(compact$verified, info = name)
    expect_identical(compact$candidate_available, candidate_available,
                     info = name)
    expect_match(compact$guidance, "refit", ignore.case = TRUE, info = name)

    if (!candidate_available) {
      expect_true(any(grepl("quarantined", output, ignore.case = TRUE)),
                  info = name)
      expect_true(all(is.na(compact[c("a", "b", "achieved_mu_K",
                                      "achieved_var_K")])), info = name)
    } else {
      expect_true(all(is.finite(unlist(compact[c(
        "a", "b", "achieved_mu_K", "achieved_var_K"
      )]))), info = name)
    }
    if (identical(name, "dual")) {
      expect_identical(compact$weight_metric, "wsb_tail")
      expect_true(is.finite(compact$weight_target_value))
      expect_true(is.finite(compact$weight_achieved_value))
      expect_true(any(grepl("Weight target", output, fixed = TRUE)))
    }
  }
})


test_that("seven migrated fits obey the two-engine plot matrix", {
  engines <- c("base", if (requireNamespace("ggplot2", quietly = TRUE)) {
    "ggplot2"
  })
  for (name in names(.s3_migrated_fits)) {
    fit <- .s3_migrated_fits[[name]]
    candidate_available <- unname(.s3_public_candidates[[name]])
    for (engine in engines) {
      value <- .s3_catch(plot(
        fit, type = "alpha", engine = engine, show = FALSE
      ))
      if (candidate_available) {
        expect_false(inherits(value, "error"),
                     info = paste(name, engine))
        if (identical(engine, "ggplot2")) {
          expect_s3_class(value, "ggplot")
        }
      } else {
        expect_s3_class(value, "dpprior_s3_unavailable_error")
        expect_s3_class(value, "dpprior_visualization_data_error")
        expect_s3_class(value, "dpprior_condition")
        expect_identical(value$code, "migrated_fit_candidate_unavailable",
                         info = paste(name, engine))
        expect_match(value$upgrade_action, "Refit", info = name)
        expect_false(inherits(value, "simpleError"), info = name)
      }
    }
  }
})


test_that("retained migrated candidates support safe basic plot routes", {
  fit <- .s3_migrated_fits$direct_a1
  engines <- c("base", if (requireNamespace("ggplot2", quietly = TRUE)) {
    "ggplot2"
  })
  for (engine in engines) {
    for (type in c("alpha", "K", "w1", "dashboard", "auto")) {
      value <- .s3_catch(plot(
        fit, type = type, engine = engine, show = FALSE
      ))
      expect_false(inherits(value, "error"), info = paste(engine, type))
    }
  }

  dual <- .s3_migrated_fits$dual
  for (type in c("dual", "comparison")) {
    condition <- .s3_catch(plot(
      dual, type = type, engine = "base", show = FALSE
    ))
    expect_s3_class(condition, "dpprior_s3_unavailable_error")
    expect_identical(
      condition$code, "migrated_comparison_lineage_unavailable"
    )
    expect_match(
      condition$message, "authoritative provenance.input_fit", fixed = TRUE
    )
  }
  expect_no_error(plot(dual, type = "auto", engine = "base", show = FALSE))
})


test_that("migrated diagnostics consume only fresh canonical components", {
  diagnostic <- .s3_upgrade(.s3_v1_diagnostics(), allow_legacy = FALSE)
  output <- capture.output(visible <- withVisible(print(diagnostic)))
  compact <- summary(diagnostic)

  expect_false(visible$visible)
  expect_identical(visible$value, diagnostic)
  expect_true(any(grepl("Schema: dpprior.result/1", output, fixed = TRUE)))
  expect_true(any(grepl("Status: CONVERGED", output, fixed = TRUE)))
  expect_true(any(grepl("alpha", output, fixed = TRUE)))
  expect_true(any(grepl("K_J", output, fixed = TRUE)))
  expect_true(any(grepl("W_SB is not W_max", output, fixed = TRUE)))
  expect_true(any(grepl("W_max: unavailable", output, fixed = TRUE)))
  expect_true(any(grepl("rho", output, fixed = TRUE)))
  expect_false(any(grepl("moderate|dominance risk:", output,
                         ignore.case = TRUE)))

  expect_s3_class(compact, "data.frame")
  expect_identical(nrow(compact), 1L)
  expect_true(all(is.finite(unlist(compact[c(
    "E_alpha", "CV_alpha", "E_K_J", "SD_K_J", "Mode_K_J",
    "E_W_SB", "P_W_SB_gt_50", "P_W_SB_gt_90", "E_rho"
  )]))))
  expect_false(compact$W_max_available)
  expect_true(is.na(compact$P_W_max_gt_50))
  expect_true(is.na(compact$P_W_max_gt_90))
  expect_match(compact$W_max_guidance, "not reconstructed")
})


test_that("canonical S3 consumers bypass forged extraction methods", {
  fit <- .s3_migrated_fits$direct_a1
  diagnostic <- .s3_upgrade(.s3_v1_diagnostics(), allow_legacy = FALSE)
  bindings <- c(
    "$.DPprior_fit", "[[.DPprior_fit",
    "$.DPprior_diagnostics", "[[.DPprior_diagnostics"
  )
  old <- lapply(bindings, get0, envir = .GlobalEnv, inherits = FALSE)
  names(old) <- bindings
  on.exit({
    for (name in bindings) {
      if (is.null(old[[name]])) {
        if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
          rm(list = name, envir = .GlobalEnv)
        }
      } else {
        assign(name, old[[name]], envir = .GlobalEnv)
      }
    }
  }, add = TRUE)
  poison <- function(...) stop("forged extraction dispatch", call. = FALSE)
  for (name in bindings) assign(name, poison, envir = .GlobalEnv)

  expect_no_error(capture.output(print(fit)))
  expect_no_error(summary(fit, print_output = FALSE))
  expect_no_error(as.data.frame(fit))
  expect_no_error(plot(fit, type = "K", engine = "base", show = FALSE))
  expect_no_error(capture.output(print(diagnostic)))
  expect_no_error(summary(diagnostic))
})


test_that("raw v1 consumers reject before legacy partial-name access", {
  operations <- list(
    function(x) print(x),
    function(x) summary(x, print_output = FALSE),
    function(x) as.data.frame(x),
    function(x) plot(x, type = "alpha", engine = "base", show = FALSE)
  )
  for (operation in operations) {
    condition <- .s3_catch(operation(.s3_v1_a2_mn()))
    expect_s3_class(condition, "dpprior_legacy_object_error")
    expect_s3_class(condition, "dpprior_serialization_error")
    expect_identical(condition$code, "legacy_upgrade_required")
    expect_false(inherits(condition, "simpleError"))
  }
  for (operation in list(print, summary)) {
    condition <- .s3_catch(operation(.s3_v1_diagnostics()))
    expect_s3_class(condition, "dpprior_legacy_object_error")
    expect_identical(condition$code, "legacy_upgrade_required")
    expect_false(inherits(condition, "simpleError"))
  }
})


.s3_current_dual_cache <- new.env(parent = emptyenv())


.s3_current_dual <- function(kind = c("soft", "hard")) {
  kind <- match.arg(kind)
  if (!exists(kind, envir = .s3_current_dual_cache, inherits = FALSE)) {
    fit <- if (identical(kind, "soft")) {
      input <- DPprior_fit(
        20L, 4, 8, method = "A2-MN", M = 80L,
        check_diagnostics = FALSE
      )
      DPprior_dual_soft(
        input,
        list(
          metric = "wsb_tail", relation = "target",
          threshold = 0.5, value = 0.3
        ),
        lambda = 1, M_fit = 80L, M_verify = 160L
      )
    } else {
      input <- DPprior_fit(
        50L, 3, 10, method = "A2-MN", M = 80L,
        check_diagnostics = FALSE
      )
      DPprior_dual_hard(
        input,
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
    }
    assign(kind, fit, envir = .s3_current_dual_cache)
  }
  unserialize(serialize(
    get(kind, envir = .s3_current_dual_cache, inherits = FALSE),
    NULL, version = 3L
  ))
}


test_that("exact canonical hard and soft class vectors remain consumable", {
  for (kind in c("hard", "soft")) {
    fit <- .s3_current_dual(kind)
    expect_silent(.dpprior_validate_result_v1(fit))
    expect_no_error(capture.output(print(fit)))
    expect_s3_class(summary(fit, print_output = FALSE), "summary.DPprior_fit")
    expect_s3_class(as.data.frame(fit), "data.frame")
    expect_no_error(plot(fit, type = "alpha", engine = "base", show = FALSE))
  }

  forged <- .s3_current_dual("soft")
  class(forged) <- c("evil", class(forged))
  condition <- .s3_catch(print(forged))
  expect_s3_class(condition, "dpprior_schema_error")
})


test_that("flat and unknown current-looking records are unsupported schemas", {
  flat <- structure(
    list(
      mode = "dual_soft", method = "dual-soft", J = 20L,
      status = "converged", usable = TRUE, verified = TRUE,
      parameters = list(a = 1.4, b = 0.7)
    ),
    class = c("DPprior_dual_soft", "DPprior_fit")
  )
  unknown <- structure(
    list(a = 1.4, b = 0.7, J = 20L, status = "converged"),
    class = "DPprior_fit"
  )
  operations <- list(
    function(x) print(x),
    function(x) summary(x, print_output = FALSE),
    function(x) as.data.frame(x),
    function(x) plot(x, type = "alpha", engine = "base", show = FALSE)
  )
  for (fit in list(flat, unknown)) {
    for (operation in operations) {
      condition <- .s3_catch(operation(fit))
      expect_s3_class(condition, "dpprior_schema_unsupported_error")
      expect_s3_class(condition, "dpprior_serialization_error")
      expect_identical(condition$code, "unsupported_schema")
    }
  }
})

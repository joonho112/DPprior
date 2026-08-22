# =============================================================================
# Module 17: S3 Methods for DPprior_fit Class
# =============================================================================
#
# This module implements enhanced S3 methods (print, summary, plot) for the
# DPprior_fit class, providing user-friendly output and visualization.
#
# Key Features:
# - print.DPprior_fit(): Concise output with estimand-labelled weight metrics
# - summary.DPprior_fit(): Detailed comparison of target vs achieved fit
# - plot.DPprior_fit(): Flexible visualization with multiple plot types
#
# Author: JoonHo Lee (jlee296@ua.edu)
# Date: December 2025
# Part of: DPprior R Package
# Dependencies: Modules 00 (constants), 08 (weights_w1), 14 (diagnostics),
#               15 (visualization)
# =============================================================================


# =============================================================================
# S3 Method: print.DPprior_fit()
# =============================================================================

.dpprior_wsb_tail_from_diagnostics <- function(weights, threshold) {
  size_biased <- weights$size_biased
  if (!is.null(size_biased$thresholds) &&
      !is.null(size_biased$tail_probability)) {
    index <- which(abs(size_biased$thresholds - threshold) <=
                     16 * .Machine$double.eps * max(1, abs(threshold)))
    if (length(index) > 0L) {
      return(as.numeric(size_biased$tail_probability[index[[1L]]]))
    }
  }
  legacy_name <- paste0("prob_gt_", threshold)
  if (!is.null(weights$prob_exceeds) &&
      legacy_name %in% names(weights$prob_exceeds)) {
    return(as.numeric(weights$prob_exceeds[[legacy_name]]))
  }
  NA_real_
}


# Gate every fit-like S3 consumer before it reads fit fields. Only validated
# dpprior.result/1 fits proceed, through an unclassed exact-name view.
.dpprior_s3_fit_gate <- function(x) {
  validated <- .dpprior_require_schema(
    x, kind = "fit", allow_legacy = FALSE
  )
  list(
    fit = validated,
    raw = unclass(validated)
  )
}


.dpprior_s3_fit_canonical_view <- function(raw) {
  target <- raw[["target", exact = TRUE]]
  target_K <- unclass(target[["K", exact = TRUE]])
  target_weight <- if ("weight" %in% names(target)) {
    unclass(target[["weight", exact = TRUE]])
  } else {
    NULL
  }
  implied <- target_K[["implied", exact = TRUE]]
  request <- target_K[["request", exact = TRUE]]
  parameters <- raw[["parameters", exact = TRUE]]
  achieved <- raw[["achieved", exact = TRUE]]
  achieved_K <- achieved[["K", exact = TRUE]]
  achieved_weight <- achieved[["weight", exact = TRUE]]
  residuals <- raw[["residuals", exact = TRUE]]
  residual_K <- residuals[["K", exact = TRUE]]
  computation <- raw[["computation", exact = TRUE]]
  termination <- computation[["termination", exact = TRUE]]
  diagnostics <- if ("diagnostics" %in% names(raw)) {
    raw[["diagnostics", exact = TRUE]]
  } else {
    NULL
  }
  provenance <- raw[["provenance", exact = TRUE]]
  migration <- provenance[["migration", exact = TRUE]]
  migrated <- !is.null(migration) && identical(
    migration[["adapter", exact = TRUE]], "upgrade_DPprior_object"
  )

  candidate_available <- !is.null(parameters) && !is.null(achieved_K)
  status <- raw[["status", exact = TRUE]]
  usable <- raw[["usable", exact = TRUE]]
  verified <- raw[["verified", exact = TRUE]]
  message <- raw[["message", exact = TRUE]]
  guidance <- if (!candidate_available && migrated) {
    paste(
      message,
      "Public candidate values are unavailable or quarantined; refit with",
      "the current API before plotting or any decision-readiness claim."
    )
  } else if (!candidate_available) {
    paste(
      message,
      "No public candidate is available; inspect the canonical status and",
      "rerun calibration with admissible inputs before plotting."
    )
  } else if (migrated &&
             (!identical(status, "converged") || !isTRUE(verified))) {
    paste(
      message,
      "This retained approximate candidate is not decision-ready; refit with",
      "the current API."
    )
  } else if (!identical(status, "converged") || !isTRUE(verified)) {
    paste(
      message,
      "This canonical candidate is not decision-ready; inspect status and",
      "verification evidence before use."
    )
  } else {
    message
  }

  requested_mean <- request[["mu_K", exact = TRUE]]
  requested_variance <- request[["var_K", exact = TRUE]]
  target_mean <- implied[["mean", exact = TRUE]]
  target_variance <- implied[["variance", exact = TRUE]]
  residual_values <- if (is.null(residual_K)) numeric() else {
    unlist(residual_K[c("mean", "variance")], use.names = FALSE)
  }

  list(
    schema = .DPPRIOR_RESULT_SCHEMA_V1,
    mode = raw[["mode", exact = TRUE]],
    method = raw[["method", exact = TRUE]],
    J = raw[["J", exact = TRUE]],
    status = status,
    usable = usable,
    verified = verified,
    message = message,
    guidance = guidance,
    migrated = migrated,
    candidate_available = candidate_available,
    a = if (is.null(parameters)) NA_real_ else {
      parameters[["a", exact = TRUE]]
    },
    b = if (is.null(parameters)) NA_real_ else {
      parameters[["b", exact = TRUE]]
    },
    target_mean = target_mean,
    target_variance = target_variance,
    requested_mean = if (is.null(requested_mean)) target_mean else {
      requested_mean
    },
    requested_variance = if (is.null(requested_variance)) {
      target_variance
    } else {
      requested_variance
    },
    weight_metric = if (is.null(target_weight)) NA_character_ else {
      target_weight[["metric", exact = TRUE]]
    },
    weight_estimand = if (is.null(target_weight)) NA_character_ else {
      target_weight[["estimand", exact = TRUE]]
    },
    weight_target_value = if (is.null(target_weight)) NA_real_ else {
      target_weight[["value", exact = TRUE]]
    },
    weight_achieved_value = if (is.null(achieved_weight)) NA_real_ else {
      achieved_weight[["value", exact = TRUE]]
    },
    achieved_mean = if (is.null(achieved_K)) NA_real_ else {
      achieved_K[["mean", exact = TRUE]]
    },
    achieved_variance = if (is.null(achieved_K)) NA_real_ else {
      achieved_K[["variance", exact = TRUE]]
    },
    residual = if (length(residual_values) && all(is.finite(residual_values))) {
      max(abs(residual_values))
    } else {
      NA_real_
    },
    M = if (is.null(achieved_K)) NULL else achieved_K[["M", exact = TRUE]],
    iterations = termination[["iterations", exact = TRUE]],
    diagnostics = diagnostics,
    legacy = if ("legacy" %in% names(raw)) {
      raw[["legacy", exact = TRUE]]
    } else {
      NULL
    }
  )
}


.dpprior_s3_abort_fit_unavailable <- function(view, operation, code) {
  message <- switch(
    code,
    migrated_comparison_lineage_unavailable = paste(
      "A migrated dual fit does not retain the K-only baseline lineage",
      "required for", operation, "; no comparison was fabricated. Refit",
      "with the current API."
    ),
    visualization_backend_unavailable = paste(
      "The visualization backend required for", operation,
      "is unavailable. Reinstall DPprior and its plotting dependencies."
    ),
    canonical_comparison_lineage_unavailable = paste(
      "The canonical fit does not expose the baseline lineage required for",
      operation, "; no comparison was fabricated."
    ),
    canonical_fit_candidate_unavailable = paste(
      "The canonical fit has no public candidate values for", operation,
      "; inspect its status/message and rerun calibration before plotting."
    ),
    paste(
      "The migrated fit has no public candidate values for", operation,
      "because legacy optimizer/selection lineage was not retained; refit",
      "with the current API."
    )
  )
  stop(.dpprior_new_condition(
    message = message,
    classes = c(
      "dpprior_s3_unavailable_error", "dpprior_visualization_data_error",
      "dpprior_error", "error"
    ),
    code = code,
    operation = operation,
    schema = view[["schema", exact = TRUE]],
    mode = view[["mode", exact = TRUE]],
    status = view[["status", exact = TRUE]],
    usable = view[["usable", exact = TRUE]],
    verified = view[["verified", exact = TRUE]],
    action = "refit_with_current_API",
    upgrade_action = if (isTRUE(view[["migrated", exact = TRUE]])) {
      paste(
        "Refit with the current API; do not recover candidate values from",
        "compatibility-only migration evidence."
      )
    } else {
      "Inspect canonical status/message and rerun calibration if appropriate."
    }
  ))
}


#' Print Method for DPprior_fit Objects
#'
#' Prints the canonical schema, mode, method, decision fields, target, and
#' available candidate evidence. Migrated or failed objects without a public
#' candidate are labelled explicitly rather than exposing quarantined values.
#'
#' @param x A \code{DPprior_fit} object.
#' @param digits Integer; number of significant digits for display.
#'   Default is 4.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @details The method first validates the complete \code{dpprior.result/1}
#' object. Printed convergence language is derived from \code{status},
#' \code{usable}, and \code{verified}; compatibility aliases and optimizer exit
#' codes are never treated as verification.
#'
#' @examples
#' # Create a fit object
#' fit <- DPprior_fit(J = 50, mu_K = 5, var_K = 8)
#' print(fit)
#'
#' # With custom digits
#' print(fit, digits = 6)
#'
#' @seealso \code{\link{summary.DPprior_fit}}, \code{\link{plot.DPprior_fit}},
#'   \code{\link{DPprior_fit}}
#'
#' @method print DPprior_fit
#' @export
print.DPprior_fit <- function(x, digits = 4, ...) {

  original <- x
  gate <- .dpprior_s3_fit_gate(x)
  view <- .dpprior_s3_fit_canonical_view(gate[["raw", exact = TRUE]])
    cat("DPprior Prior Elicitation Result\n")
    cat(strrep("=", 45), "\n\n")
    cat(sprintf("Schema: %s\n", view[["schema", exact = TRUE]]))
    cat(sprintf("Method: %s (mode: %s)\n",
                view[["method", exact = TRUE]],
                view[["mode", exact = TRUE]]))
    cat(sprintf("Status: %s; usable: %s; verified: %s\n\n",
                view[["status", exact = TRUE]],
                if (isTRUE(view[["usable", exact = TRUE]])) "yes" else "no",
                if (isTRUE(view[["verified", exact = TRUE]])) "yes" else "no"))
    cat(sprintf("Target (J = %d):\n", view[["J", exact = TRUE]]))
    cat(sprintf("  E[K_J]   = %.4f\n",
                view[["target_mean", exact = TRUE]]))
    cat(sprintf("  Var(K_J) = %.4f\n",
                view[["target_variance", exact = TRUE]]))

    if (isTRUE(view[["candidate_available", exact = TRUE]])) {
      if (isTRUE(view[["migrated", exact = TRUE]])) {
        cat("\nRetained candidate (approximate migration evidence):\n")
      } else {
        cat("\nCanonical candidate:\n")
      }
      cat(sprintf("  alpha ~ Gamma(a = %.*f, b = %.*f)\n",
                  digits, view[["a", exact = TRUE]],
                  digits, view[["b", exact = TRUE]]))
      cat(sprintf("  Achieved E[K_J] = %.*f; Var(K_J) = %.*f\n",
                  digits + 2L, view[["achieved_mean", exact = TRUE]],
                  digits + 2L, view[["achieved_variance", exact = TRUE]]))
      if (is.finite(view[["residual", exact = TRUE]])) {
        cat(sprintf("  Maximum absolute moment residual = %.2e\n",
                    view[["residual", exact = TRUE]]))
      }
      if (is.finite(view[["weight_target_value", exact = TRUE]])) {
        cat(sprintf("  Weight target (%s) = %.6g; achieved = %.6g\n",
                    view[["weight_metric", exact = TRUE]],
                    view[["weight_target_value", exact = TRUE]],
                    view[["weight_achieved_value", exact = TRUE]]))
      }
    } else {
      if (isTRUE(view[["migrated", exact = TRUE]])) {
        cat(paste0(
          "\nPublic candidate: unavailable (legacy A2 candidate is ",
          "quarantined).\n"
        ))
      } else {
        cat("\nPublic candidate: unavailable under the canonical status.\n")
      }
    }
    cat(sprintf("\n%s guidance: %s\n",
                if (isTRUE(view[["migrated", exact = TRUE]])) {
                  "Migration"
                } else {
                  "Canonical"
                },
                view[["guidance", exact = TRUE]]))
  invisible(original)
}


# =============================================================================
# S3 Method: summary.DPprior_fit()
# =============================================================================

#' Summary Method for DPprior_fit Objects
#'
#' Produces a status-aware summary from validated canonical fields, including
#' candidate availability, Gamma parameters, target-versus-achieved evidence,
#' and any attached diagnostics.
#'
#' @param object A \code{DPprior_fit} object.
#' @param print_output Logical; if \code{TRUE} (default), prints the summary
#'   to the console. If \code{FALSE}, returns the summary list silently.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class \code{"summary.DPprior_fit"} containing:
#'   \itemize{
#'     \item \code{schema}, \code{mode}, \code{method}, \code{status},
#'       \code{usable}, and \code{verified};
#'     \item \code{candidate_available}, \code{gamma_prior}, and derived
#'       \code{alpha_summary};
#'     \item canonical \code{target}, \code{achieved}, \code{weight}, and
#'       \code{errors} views;
#'     \item \code{message}, \code{guidance}, \code{iterations}, and attached
#'       \code{diagnostics}.
#'   }
#'
#' @details The complete input is schema-validated before summarization. A
#' migrated result without a public candidate remains visibly unavailable;
#' quarantined compatibility evidence is not promoted into the summary.
#'
#' @examples
#' # Create a fit object
#' fit <- DPprior_fit(J = 50, mu_K = 5, var_K = 8, check_diagnostics = TRUE)
#' summary(fit)
#'
#' # Store summary without printing
#' summ <- summary(fit, print_output = FALSE)
#' str(summ)
#'
#' @seealso \code{\link{print.DPprior_fit}}, \code{\link{DPprior_diagnostics}}
#'
#' @method summary DPprior_fit
#' @export
summary.DPprior_fit <- function(object, print_output = TRUE, ...) {

  gate <- .dpprior_s3_fit_gate(object)
  view <- .dpprior_s3_fit_canonical_view(gate[["raw", exact = TRUE]])
    candidate_available <- view[["candidate_available", exact = TRUE]]
    a <- view[["a", exact = TRUE]]
    b <- view[["b", exact = TRUE]]
    target_mean <- view[["target_mean", exact = TRUE]]
    target_variance <- view[["target_variance", exact = TRUE]]
    achieved_mean <- view[["achieved_mean", exact = TRUE]]
    achieved_variance <- view[["achieved_variance", exact = TRUE]]
    mean_error <- if (candidate_available) {
      abs(target_mean - achieved_mean)
    } else {
      NA_real_
    }
    variance_error <- if (candidate_available) {
      abs(target_variance - achieved_variance)
    } else {
      NA_real_
    }
    result <- list(
      schema = view[["schema", exact = TRUE]],
      canonical = TRUE,
      migrated = view[["migrated", exact = TRUE]],
      mode = view[["mode", exact = TRUE]],
      method = view[["method", exact = TRUE]],
      status = view[["status", exact = TRUE]],
      usable = view[["usable", exact = TRUE]],
      verified = view[["verified", exact = TRUE]],
      message = view[["message", exact = TRUE]],
      guidance = view[["guidance", exact = TRUE]],
      candidate_available = candidate_available,
      gamma_prior = list(a = a, b = b),
      alpha_summary = list(
        E_alpha = if (candidate_available) a / b else NA_real_,
        Var_alpha = if (candidate_available) a / b^2 else NA_real_,
        SD_alpha = if (candidate_available) sqrt(a) / b else NA_real_,
        CV_alpha = if (candidate_available) 1 / sqrt(a) else NA_real_
      ),
      target = list(
        mu_K = view[["requested_mean", exact = TRUE]],
        var_K = view[["requested_variance", exact = TRUE]],
        var_K_used = target_variance,
        confidence = NA_character_
      ),
      achieved = list(
        mu_K = achieved_mean,
        var_K = achieved_variance,
        residual = view[["residual", exact = TRUE]]
      ),
      weight = list(
        metric = view[["weight_metric", exact = TRUE]],
        estimand = view[["weight_estimand", exact = TRUE]],
        target = view[["weight_target_value", exact = TRUE]],
        achieved = view[["weight_achieved_value", exact = TRUE]]
      ),
      errors = list(
        mu_K_abs = mean_error,
        var_K_abs = variance_error,
        mu_K_rel_pct = if (is.finite(mean_error) && target_mean > 0) {
          100 * mean_error / target_mean
        } else {
          NA_real_
        },
        var_K_rel_pct = if (is.finite(variance_error) &&
                            target_variance > 0) {
          100 * variance_error / target_variance
        } else {
          NA_real_
        }
      ),
      scaling = list(J = view[["J", exact = TRUE]]),
      converged = identical(view[["status", exact = TRUE]], "converged") &&
        isTRUE(view[["usable", exact = TRUE]]) &&
        isTRUE(view[["verified", exact = TRUE]]),
      iterations = view[["iterations", exact = TRUE]],
      diagnostics = view[["diagnostics", exact = TRUE]],
      trace = NULL,
      dual_anchor = view[["legacy", exact = TRUE]]
    )
    class(result) <- "summary.DPprior_fit"
    if (isTRUE(print_output)) {
      print(result)
    }
  invisible(result)
}


#' Print Method for summary.DPprior_fit
#'
#' @param x A \code{summary.DPprior_fit} object.
#' @param diagnostics Logical; if TRUE, print full diagnostics. Default is FALSE.
#' @param max_trace Integer; maximum number of trace rows to display. Default is 10.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @method print summary.DPprior_fit
#' @export
print.summary.DPprior_fit <- function(x, diagnostics = FALSE, max_trace = 10L, ...) {

  if (isTRUE(x[["canonical", exact = TRUE]])) {
    cat("DPprior Prior Elicitation Summary\n")
    cat(strrep("=", 60), "\n\n")
    cat(sprintf("Schema: %s\n", x[["schema", exact = TRUE]]))
    cat(sprintf("Sample size: J = %d\n",
                x[["scaling", exact = TRUE]][["J", exact = TRUE]]))
    cat(sprintf("Method: %s (mode: %s)\n",
                x[["method", exact = TRUE]], x[["mode", exact = TRUE]]))
    cat(sprintf("Status: %s; usable: %s; verified: %s\n\n",
                x[["status", exact = TRUE]],
                if (isTRUE(x[["usable", exact = TRUE]])) "yes" else "no",
                if (isTRUE(x[["verified", exact = TRUE]])) "yes" else "no"))

    target <- x[["target", exact = TRUE]]
    cat("Canonical K_J target:\n")
    cat(strrep("-", 40), "\n")
    cat(sprintf("  Requested E[K_J] = %.4f\n",
                target[["mu_K", exact = TRUE]]))
    if (is.finite(target[["var_K", exact = TRUE]]) &&
        abs(target[["var_K", exact = TRUE]] -
            target[["var_K_used", exact = TRUE]]) > 1e-10) {
      cat(sprintf("  Requested Var(K_J) = %.4f; canonical used = %.4f\n",
                  target[["var_K", exact = TRUE]],
                  target[["var_K_used", exact = TRUE]]))
    } else {
      cat(sprintf("  Var(K_J) = %.4f\n",
                  target[["var_K_used", exact = TRUE]]))
    }

    if (isTRUE(x[["candidate_available", exact = TRUE]])) {
      gamma <- x[["gamma_prior", exact = TRUE]]
      alpha <- x[["alpha_summary", exact = TRUE]]
      achieved <- x[["achieved", exact = TRUE]]
      errors <- x[["errors", exact = TRUE]]
      if (isTRUE(x[["migrated", exact = TRUE]])) {
        cat("\nRetained candidate (approximate migration evidence):\n")
      } else {
        cat("\nCanonical candidate:\n")
      }
      cat(strrep("-", 40), "\n")
      cat(sprintf("  Gamma(a = %.6f, b = %.6f)\n",
                  gamma[["a", exact = TRUE]], gamma[["b", exact = TRUE]]))
      cat(sprintf("  E[alpha] = %.4f; SD[alpha] = %.4f; CV[alpha] = %.4f\n",
                  alpha[["E_alpha", exact = TRUE]],
                  alpha[["SD_alpha", exact = TRUE]],
                  alpha[["CV_alpha", exact = TRUE]]))
      cat(sprintf("  Achieved E[K_J] = %.4f; Var(K_J) = %.4f\n",
                  achieved[["mu_K", exact = TRUE]],
                  achieved[["var_K", exact = TRUE]]))
      cat(sprintf("  Absolute errors: mean %.3e; variance %.3e\n",
                  errors[["mu_K_abs", exact = TRUE]],
                  errors[["var_K_abs", exact = TRUE]]))
      weight <- x[["weight", exact = TRUE]]
      if (is.finite(weight[["target", exact = TRUE]])) {
        cat(sprintf("  Weight target (%s) = %.6g; achieved = %.6g\n",
                    weight[["metric", exact = TRUE]],
                    weight[["target", exact = TRUE]],
                    weight[["achieved", exact = TRUE]]))
      }
    } else {
      if (isTRUE(x[["migrated", exact = TRUE]])) {
        cat(paste0(
          "\nPublic candidate: unavailable. The legacy A2 parameter pair is ",
          "quarantined and was not exposed as a successful fit.\n"
        ))
      } else {
        cat("\nPublic candidate: unavailable under the canonical status.\n")
      }
    }
    cat(sprintf("\n%s guidance: %s\n",
                if (isTRUE(x[["migrated", exact = TRUE]])) {
                  "Migration"
                } else {
                  "Canonical"
                },
                x[["guidance", exact = TRUE]]))
    nested_diagnostics <- x[["diagnostics", exact = TRUE]]
    if (isTRUE(diagnostics) && !is.null(nested_diagnostics)) {
      cat("\nCanonical nested diagnostics:\n")
      print(nested_diagnostics)
    }
    invisible(x)
  }
}


#' Coerce a canonical DPprior fit to a data frame
#'
#' Returns a one-row, status-aware public view of a canonical
#' \code{dpprior.result/1} fit. The method validates the complete object before
#' extracting any value; it never repairs a result or reads unregistered flat
#' fields as scientific evidence.
#'
#' @param x A canonical \code{DPprior_fit} object.
#' @param row.names Optional row names.
#' @param optional Retained for the data-frame S3 signature; it does not alter
#'   the canonical column contract.
#' @param ... Additional arguments, currently unused.
#'
#' @return A one-row data frame containing schema, mode, method, status,
#'   \code{usable}, \code{verified}, candidate availability, Gamma parameters,
#'   canonical target and achieved K moments, any named weight target evidence,
#'   and status guidance. Parameter columns are \code{NA} when the canonical
#'   result has no public candidate.
#'
#' @examples
#' fit <- DPprior_fit(
#'   J = 50, mu_K = 5, var_K = 8, check_diagnostics = FALSE
#' )
#' as.data.frame(fit)
#'
#' @seealso \code{\link{DPprior_fit}}, \code{\link{summary.DPprior_fit}}
#' @method as.data.frame DPprior_fit
#' @export
as.data.frame.DPprior_fit <- function(x, row.names = NULL,
                                      optional = FALSE, ...) {
  gate <- .dpprior_s3_fit_gate(x)
  view <- .dpprior_s3_fit_canonical_view(gate[["raw", exact = TRUE]])
  candidate_available <- view[["candidate_available", exact = TRUE]]
  data.frame(
    schema = view[["schema", exact = TRUE]],
    migrated = view[["migrated", exact = TRUE]],
    mode = view[["mode", exact = TRUE]],
    method = view[["method", exact = TRUE]],
    status = view[["status", exact = TRUE]],
    usable = view[["usable", exact = TRUE]],
    verified = view[["verified", exact = TRUE]],
    candidate_available = candidate_available,
    a = view[["a", exact = TRUE]],
    b = view[["b", exact = TRUE]],
    J = view[["J", exact = TRUE]],
    mu_K = view[["target_mean", exact = TRUE]],
    var_K = view[["target_variance", exact = TRUE]],
    achieved_mu_K = view[["achieved_mean", exact = TRUE]],
    achieved_var_K = view[["achieved_variance", exact = TRUE]],
    weight_metric = view[["weight_metric", exact = TRUE]],
    weight_target_value = view[["weight_target_value", exact = TRUE]],
    weight_achieved_value = view[["weight_achieved_value", exact = TRUE]],
    mean_alpha = if (candidate_available) {
      view[["a", exact = TRUE]] / view[["b", exact = TRUE]]
    } else {
      NA_real_
    },
    cv_alpha = if (candidate_available) {
      1 / sqrt(view[["a", exact = TRUE]])
    } else {
      NA_real_
    },
    scaling = NA_character_,
    converged = identical(view[["status", exact = TRUE]], "converged") &&
      isTRUE(view[["usable", exact = TRUE]]) &&
      isTRUE(view[["verified", exact = TRUE]]),
    iterations = if (is.null(view[["iterations", exact = TRUE]])) {
      NA_integer_
    } else {
      view[["iterations", exact = TRUE]]
    },
    message = view[["message", exact = TRUE]],
    guidance = view[["guidance", exact = TRUE]],
    stringsAsFactors = FALSE,
    row.names = row.names,
    check.names = FALSE
  )
}


# =============================================================================
# S3 Method: plot.DPprior_fit()
# =============================================================================

#' Plot Method for DPprior_fit Objects
#'
#' Creates visualizations of a prior elicitation result. Multiple plot types
#' are available, including individual distribution plots and comprehensive
#' dashboards.
#'
#' @param x A \code{DPprior_fit} object.
#' @param type Character; the type of plot to create:
#'   \describe{
#'     \item{"auto"}{(Default) Automatically selects the appropriate plot type.
#'       Uses \code{"dual"} for current hard and soft dual fits. A1, A2, and
#'       retained legacy fits use the descriptive single-fit dashboard.}
#'     \item{"dashboard"}{4-panel dashboard showing alpha, K, the first
#'       size-biased weight \eqn{W_{SB}}, and summary.}
#'     \item{"alpha"}{Prior density of the concentration parameter alpha.}
#'     \item{"K"}{Prior PMF of the number of clusters \eqn{K_J}.}
#'     \item{"w1"}{Prior density of the first size-biased DP weight
#'       \eqn{W_{SB}}.}
#'     \item{"dual"}{Authoritative K-only comparison for current hard or soft
#'       dual fits. Legacy fits raise a typed lineage-unavailable condition.}
#'     \item{"comparison"}{Same as "dual".}
#'   }
#' @param engine Character; graphics engine to use:
#'   \code{"ggplot2"} (default) or \code{"base"}.
#' @param ... Additional arguments passed to the underlying plot functions.
#'   Common options include:
#'   \describe{
#'     \item{base_size}{Base font size (default: 11)}
#'     \item{ci_level}{Credible interval level for alpha plot (default: 0.95)}
#'     \item{title}{Optional title for the dashboard}
#'     \item{show}{If TRUE, display the plot; if FALSE, return silently}
#'   }
#'
#' @return Depends on the plot type and engine:
#'   \itemize{
#'     \item For ggplot2: Returns a ggplot object or gtable (for dashboards)
#'     \item For base: Returns invisible(NULL)
#'   }
#'
#' @details
#' The \code{"auto"} type is recommended for most use cases. It automatically
#' uses an authoritative comparison only for current hard or soft calibration.
#' Retained legacy results remain available as status-labelled descriptive
#' single-fit plots; no compatibility initialization is used as comparison
#' science.
#'
#' For dual-anchor fits, the comparison dashboard shows:
#' \itemize{
#'   \item Alpha prior: K-only vs Dual-anchor
#'   \item K distribution comparison
#'   \item \eqn{W_{SB}} distribution comparison with named tail thresholds
#'   \item Summary comparison table
#' }
#'
#' @section Plot Type Details:
#' \describe{
#'   \item{dashboard}{
#'     A 2x2 grid showing:
#'     (A) Alpha prior density with CI
#'     (B) \eqn{K_J} prior PMF with mode and mean
#'     (C) first size-biased weight density with threshold shading
#'     (D) Summary statistics table
#'   }
#'   \item{alpha}{
#'     Gamma(a, b) density with:
#'     - Mean line (dashed)
#'     - Credible interval (shaded region)
#'     - Annotation with moments and CI
#'   }
#'   \item{K}{
#'     Bar plot of \eqn{P(K_J = k)} with:
#'     - Target mean line
#'     - Achieved mean line
#'     - Optional CDF overlay
#'   }
#'   \item{w1}{
#'     Density plot with:
#'     - Explicit threshold-region shading for \eqn{W_{SB}}
#'     - Threshold lines
#'     - Estimand-labelled exceedance probabilities
#'   }
#' }
#'
#' @examples
#' # Create a fit object
#' fit <- DPprior_fit(J = 50, mu_K = 5, var_K = 8)
#'
#' # Auto-detect best plot type
#' plot(fit)
#'
#' # Specific plot types
#' plot(fit, type = "alpha")
#' plot(fit, type = "K")
#' plot(fit, type = "w1")
#' plot(fit, type = "dashboard")
#'
#' # With custom options
#' plot(fit, type = "dashboard", title = "My Prior Analysis")
#'
#' # Current soft dual comparison
#' fit_K <- DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8)
#' fit_dual <- DPprior_dual_soft(
#'   fit_K,
#'   target = list(
#'     metric = "wsb_tail", relation = "target",
#'     threshold = 0.5, value = 0.3
#'   ),
#'   lambda = 0.5
#' )
#' plot(fit_dual)  # Auto-selects dual comparison
#' plot(fit_dual, type = "comparison")  # Explicit
#'
#' # Retained legacy results use a descriptive single-fit dashboard.
#' fit_legacy <- suppressWarnings(DPprior_dual(
#'   fit_K, list(prob = list(threshold = 0.5, value = 0.3))
#' ))
#' plot(fit_legacy)
#' try(plot(fit_legacy, type = "comparison")) # typed lineage unavailable
#'
#' @seealso \code{\link{plot_prior_dashboard}}, \code{\link{plot_alpha_prior}},
#'   \code{\link{plot_K_prior}}, \code{\link{plot_w1_prior}},
#'   \code{\link{plot_dual_comparison}}
#'
#' @method plot DPprior_fit
#' @export
plot.DPprior_fit <- function(x, type = c("auto", "dashboard", "alpha", "K", "w1",
                                         "dual", "comparison"),
                             engine = c("ggplot2", "base"),
                             ...) {
  gate <- .dpprior_s3_fit_gate(x)
  type <- match.arg(type)
  engine <- match.arg(engine)
  view <- .dpprior_s3_fit_canonical_view(gate[["raw", exact = TRUE]])
  if (!isTRUE(view[["candidate_available", exact = TRUE]])) {
    .dpprior_s3_abort_fit_unavailable(
      view, paste0("plot(type='", type, "')"),
      if (isTRUE(view[["migrated", exact = TRUE]])) {
        "migrated_fit_candidate_unavailable"
      } else {
        "canonical_fit_candidate_unavailable"
      }
    )
  }
  mode <- view[["mode", exact = TRUE]]
  if (identical(type, "auto")) {
    type <- if (mode %in% c("dual_hard", "dual_soft")) {
      "dual"
    } else {
      # A1, A2, and retained legacy fits are descriptive single-fit views.
      # Legacy comparison is unavailable because it has no authoritative
      # provenance.input_fit lineage.
      "dashboard"
    }
  }
  fit <- gate[["fit", exact = TRUE]]
  backend <- switch(
    type,
    dashboard = "plot_prior_dashboard",
    alpha = "plot_alpha_prior",
    K = "plot_K_prior",
    w1 = "plot_w1_prior",
    dual = "plot_dual_comparison",
    comparison = "plot_dual_comparison"
  )
  if (!exists(backend, mode = "function")) {
    .dpprior_s3_abort_fit_unavailable(
      view, paste0("plot(type='", type, "')"),
      "visualization_backend_unavailable"
    )
  }
  switch(
    type,
    dashboard = plot_prior_dashboard(fit, engine = engine, ...),
    alpha = plot_alpha_prior(fit = fit, engine = engine, ...),
    K = plot_K_prior(fit = fit, engine = engine, ...),
    w1 = plot_w1_prior(fit = fit, engine = engine, ...),
    dual = plot_dual_comparison(fit_dual = fit, engine = engine, ...),
    comparison = plot_dual_comparison(
      fit_dual = fit, engine = engine, ...
    )
  )
}


# =============================================================================
# Helper Functions for S3 Methods
# =============================================================================

#' Check whether a fit uses a dual-anchor mode
#'
#' Schema-validates a DPprior_fit object and recognizes current hard, current
#' soft, and retained legacy dual modes.
#'
#' @param fit A DPprior_fit object.
#' @return Logical; TRUE if dual-anchor fit with valid structure.
#'
#' @keywords internal
.dpprior_is_dual <- function(fit) {
  if (!inherits(fit, "DPprior_fit")) return(FALSE)
  gate <- .dpprior_s3_fit_gate(fit)
  gate[["raw", exact = TRUE]][["mode", exact = TRUE]] %in% c(
    "dual_hard", "dual_soft", "dual_legacy"
  )
}


# =============================================================================
# Verification Function
# =============================================================================

#' Verify S3 Methods Module
#'
#' Runs comprehensive verification tests for the S3 methods module.
#'
#' @param verbose Logical; if TRUE, print detailed test output.
#'
#' @return Invisibly returns TRUE if all tests pass.
#'
#' @examples
#' \dontrun{
#' verify_s3_methods()
#'
#' }
#' @keywords internal
verify_s3_methods <- function(verbose = TRUE) {
  fit <- DPprior_fit(
    20L, mu_K = 4, var_K = 8, method = "A2-MN", M = 80L,
    check_diagnostics = FALSE
  )

  checks <- list(
    schema = tryCatch({
      .dpprior_require_schema(fit, kind = "fit", allow_legacy = FALSE)
      TRUE
    }, error = function(error) FALSE),
    print = tryCatch({
      output <- capture.output(print(fit))
      any(grepl("Schema: dpprior.result/1", output, fixed = TRUE))
    }, error = function(error) FALSE),
    summary = tryCatch({
      result <- summary(fit, print_output = FALSE)
      output <- capture.output(print(result))
      inherits(result, "summary.DPprior_fit") &&
        identical(result[["schema", exact = TRUE]], .DPPRIOR_RESULT_SCHEMA_V1) &&
        any(grepl("Canonical K_J target", output, fixed = TRUE))
    }, error = function(error) FALSE),
    data_frame = tryCatch({
      result <- as.data.frame(fit)
      identical(nrow(result), 1L) &&
        identical(result[["schema", exact = TRUE]], .DPPRIOR_RESULT_SCHEMA_V1)
    }, error = function(error) FALSE),
    mode = tryCatch(!.dpprior_is_dual(fit), error = function(error) FALSE),
    plot = tryCatch({
      plot(fit, type = "alpha", engine = "base", show = FALSE)
      TRUE
    }, error = function(error) FALSE)
  )

  passed <- vapply(checks, isTRUE, logical(1L))
  if (isTRUE(verbose)) {
    cat(strrep("=", 60), "\n")
    cat("Module 17: canonical S3 methods verification\n")
    cat(strrep("=", 60), "\n")
    for (name in names(passed)) {
      cat(sprintf("  %-12s %s\n", name, if (passed[[name]]) "PASS" else "FAIL"))
    }
  }
  invisible(all(passed))
}

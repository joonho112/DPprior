#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats density optim qgamma
#' @importFrom utils capture.output
## usethis namespace: end
NULL

# Global variables used in ggplot2 aes() and data.table/dplyr operations
# These are column names referenced non-standardly and are not actual global variables
utils::globalVariables(c(

  # ggplot2 aesthetics
  "density", "k", "pmf", "lambda", "mu_K", "var_K",
  "Method", "Metric", "Value", "Component", "Weight",
  "K_only", "Dual",


  # Internal function references
  "exact_K_pmf",

  # rlang .data pronoun
  ".data"
))

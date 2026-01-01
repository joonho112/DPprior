# Find Threshold J for A1 Adequacy

Determines the minimum sample size J for which the A1 approximation
achieves target accuracy in moment matching.

## Usage

``` r
find_a1_threshold_J(
  a,
  b,
  target_error = 0.05,
  target_var_error = NULL,
  J_min = 10,
  J_max = 500,
  step = 10
)
```

## Arguments

- a, b:

  Numeric; Gamma hyperparameters.

- target_error:

  Numeric; target relative error for mean (default: 5%).

- target_var_error:

  Numeric; target relative error for variance (default: 2 \*
  target_error).

- J_min, J_max:

  Integer; search range for J.

- step:

  Integer; step size for search.

## Value

Integer; minimum J achieving target accuracy, or NA if not found.

## Details

Searches over J values to find the smallest J where:

- Mean relative error \< target_error

- Variance relative error \< target_var_error

Note: For many parameter combinations, especially with high
\\E\[\alpha\]\\, the A1 approximation may never achieve low errors
within practical J ranges. In such cases, A2 refinement is recommended.

## See also

[`a1_moment_error`](https://joonho112.github.io/DPprior/reference/a1_moment_error.md),
[`DPprior_error_bounds`](https://joonho112.github.io/DPprior/reference/DPprior_error_bounds.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Find threshold for 5% mean error
find_a1_threshold_J(a = 1, b = 2)

# For higher E[alpha], threshold may not exist
find_a1_threshold_J(a = 2, b = 1)  # Likely returns NA
} # }
```

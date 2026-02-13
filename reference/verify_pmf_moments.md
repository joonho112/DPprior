# Verify PMF Consistency with Moments

Verifies that moments computed from the PMF match the closed-form
digamma/trigamma formulas.

## Usage

``` r
verify_pmf_moments(J, alpha, logS, tol = 1e-08, verbose = TRUE)
```

## Arguments

- J:

  Integer; sample size.

- alpha:

  Numeric; DP concentration parameter.

- logS:

  Matrix; pre-computed log-Stirling matrix.

- tol:

  Numeric; tolerance for comparison (default: 1e-8).

- verbose:

  Logical; if `TRUE`, print results.

## Value

Logical; `TRUE` if verification passes.

## Examples

``` r
if (FALSE) { # \dontrun{
logS <- compute_log_stirling(50)
verify_pmf_moments(50, 2.0, logS)

} # }
```

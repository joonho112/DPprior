# Verify A2-MN Moment Matching

Tests that the A2-MN solver achieves exact moment matching.

## Usage

``` r
verify_a2_moment_matching(J, mu_K, var_K, tol = 1e-06, verbose = TRUE)
```

## Arguments

- J:

  Integer; sample size.

- mu_K:

  Numeric; target mean.

- var_K:

  Numeric; target variance.

- tol:

  Numeric; tolerance for verification.

- verbose:

  Logical; if TRUE, print results.

## Value

Logical; TRUE if verification passes.

## Examples

``` r
if (FALSE) { # \dontrun{
verify_a2_moment_matching(J = 50, mu_K = 5, var_K = 8)

} # }
```

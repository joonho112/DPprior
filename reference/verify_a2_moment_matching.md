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
verify_a2_moment_matching(J = 50, mu_K = 5, var_K = 8)
#> A2-MN Moment Matching Verification
#> -------------------------------------------------- 
#> Target: E[K]=5.0000, Var(K)=8.0000
#> Achieved: E[K]=4.9999999992, Var(K)=8.0000000076
#> Mean error: 8.31e-10
#> Var error: 7.55e-09
#> Termination: residual
#> Status: PASS
```

# Verify Derivative via Finite Difference

Verifies the analytic derivative of the conditional mean against finite
difference approximation.

## Usage

``` r
verify_derivative(J, alpha, eps = 1e-06, tol = 1e-05, verbose = TRUE)
```

## Arguments

- J:

  Sample size (integer \>= 1).

- alpha:

  Concentration parameter (positive scalar).

- eps:

  Finite difference step size (default: 1e-6).

- tol:

  Tolerance for comparison (default: 1e-5).

- verbose:

  Logical; if TRUE, print results.

## Value

Logical; TRUE if derivative matches finite difference.

## Examples

``` r
verify_derivative(50, 2.0)
#> Derivative verification (J=50, alpha=2.00):
#>   Analytic:    2.2677787792
#>   Finite diff: 2.2677787790
#>   Error:       2.13e-10 [PASS]
```

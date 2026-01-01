# Verify Conditional Variance Formula

Verifies that `Var(rho|alpha) = E(rho^2|alpha) - E(rho|alpha)^2`.

## Usage

``` r
verify_rho_conditional_variance(alpha, tol = 1e-12)
```

## Arguments

- alpha:

  Numeric; concentration parameter (must be positive).

- tol:

  Numeric; tolerance for comparison. Default is 1e-12.

## Value

Logical; TRUE if formula holds.

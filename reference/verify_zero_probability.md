# Verify Zero Probability at K=0

Verifies that P(K_J = 0 \| alpha) = 0.

## Usage

``` r
verify_zero_probability(J, alpha, logS, tol = 1e-15, verbose = TRUE)
```

## Arguments

- J:

  Integer; sample size.

- alpha:

  Numeric; DP concentration parameter.

- logS:

  Matrix; pre-computed log-Stirling matrix.

- tol:

  Numeric; tolerance (default: 1e-15).

- verbose:

  Logical; if `TRUE`, print results.

## Value

Logical; `TRUE` if verification passes.

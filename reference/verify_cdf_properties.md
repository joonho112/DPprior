# Verify CDF Properties

Verifies that the CDF is non-decreasing and ends at 1.

## Usage

``` r
verify_cdf_properties(J, alpha, logS, tol = 1e-10, verbose = TRUE)
```

## Arguments

- J:

  Integer; sample size.

- alpha:

  Numeric; DP concentration parameter.

- logS:

  Matrix; pre-computed log-Stirling matrix.

- tol:

  Numeric; tolerance (default: 1e-10).

- verbose:

  Logical; if `TRUE`, print results.

## Value

Logical; `TRUE` if verification passes.

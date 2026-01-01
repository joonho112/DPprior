# Validate Conditional Moments Computation

Validates polygamma formulas against direct summation.

## Usage

``` r
validate_moments_conditional(J, alpha, tol = 1e-10, verbose = TRUE)
```

## Arguments

- J:

  Sample size.

- alpha:

  Concentration parameter.

- tol:

  Tolerance (default: 1e-10).

- verbose:

  Print results if TRUE.

## Value

Logical; TRUE if validation passes.

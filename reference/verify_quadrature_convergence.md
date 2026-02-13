# Verify Quadrature Convergence for Marginal Moments

Tests that marginal moments converge as the number of quadrature nodes
increases.

## Usage

``` r
verify_quadrature_convergence(
  J,
  a,
  b,
  M_values = c(20, 40, 60, 80, 100, 120),
  verbose = TRUE
)
```

## Arguments

- J:

  Integer; sample size.

- a:

  Numeric; shape parameter.

- b:

  Numeric; rate parameter.

- M_values:

  Integer vector; numbers of quadrature nodes to test.

- verbose:

  Logical; if `TRUE`, print detailed results.

## Value

Data frame with convergence results.

## Examples

``` r
if (FALSE) { # \dontrun{
verify_quadrature_convergence(50, 1.5, 0.5)

} # }
```

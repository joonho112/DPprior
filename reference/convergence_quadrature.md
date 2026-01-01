# Convergence Diagnostic for Quadrature

Examines how quadrature accuracy improves with increasing number of
nodes.

## Usage

``` r
convergence_quadrature(
  f,
  a,
  b,
  M_values = c(10, 20, 50, 80, 100),
  true_value = NULL
)
```

## Arguments

- f:

  Function to integrate.

- a:

  Numeric; shape parameter of Gamma distribution.

- b:

  Numeric; rate parameter of Gamma distribution.

- M_values:

  Integer vector; number of nodes to test.

- true_value:

  Numeric; known true value (optional).

## Value

A data frame with columns: M, estimate, change, relative_change.

## Examples

``` r
# Check convergence for E[alpha]
convergence_quadrature(identity, 2.5, 1.5,
                       M_values = c(10, 20, 50, 80, 100),
                       true_value = 2.5/1.5)
#>     M estimate        change relative_change        error relative_error
#> 1  10 1.666667            NA              NA 4.440892e-16   2.664535e-16
#> 2  20 1.666667 -1.332268e-15    7.993606e-16 8.881784e-16   5.329071e-16
#> 3  50 1.666667  0.000000e+00    0.000000e+00 8.881784e-16   5.329071e-16
#> 4  80 1.666667  4.440892e-16    2.664535e-16 4.440892e-16   2.664535e-16
#> 5 100 1.666667  4.440892e-16    2.664535e-16 0.000000e+00   0.000000e+00
```

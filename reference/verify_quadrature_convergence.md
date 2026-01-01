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
verify_quadrature_convergence(50, 1.5, 0.5)
#> Quadrature Convergence (J=50, a=1.50, b=0.50):
#>    M     mean      var  mean_change   var_change
#>   20 8.355497 22.76871           NA           NA
#>   40 8.355487 22.76895 1.052945e-05 2.394687e-04
#>   60 8.355487 22.76895 5.804876e-08 1.314309e-06
#>   80 8.355487 22.76895 1.056184e-09 2.390260e-08
#>  100 8.355487 22.76895 3.567457e-11 8.072334e-10
#>  120 8.355487 22.76895 1.790568e-12 4.058620e-11
#> 
#> Mean: CONVERGED to machine precision
```

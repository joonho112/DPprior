# Verify Quadrature Convergence

Checks that the marginal PMF converges as the number of quadrature nodes
increases.

## Usage

``` r
verify_pmf_marginal_convergence(
  J,
  a,
  b,
  logS,
  M_values = c(20L, 40L, 80L, 120L),
  verbose = TRUE
)
```

## Arguments

- J:

  Integer; sample size.

- a:

  Numeric; shape parameter of Gamma prior.

- b:

  Numeric; rate parameter of Gamma prior.

- logS:

  Matrix; pre-computed log-Stirling matrix.

- M_values:

  Integer vector; quadrature node counts to test.

- verbose:

  Logical; if TRUE, print results.

## Value

Data frame with convergence results.

## Examples

``` r
logS <- compute_log_stirling(50)
verify_pmf_marginal_convergence(50, 1.5, 0.5, logS)
#> Quadrature Convergence (J=50, a=1.50, b=0.50):
#>      M         mean          var     mean_chg      var_chg       L1_chg
#>     20   8.35549735  22.76870932     0.00e+00     0.00e+00     0.00e+00
#>     40   8.35548682  22.76894879     1.05e-05     2.39e-04     4.09e-03
#>     80   8.35548676  22.76895012     5.91e-08     1.34e-06     7.89e-05
#>    120   8.35548676  22.76895013     3.75e-11     8.48e-10     8.62e-08
```

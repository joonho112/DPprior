# Verify Marginal PMF Properties

Verifies that the marginal PMF satisfies basic probability properties.

## Usage

``` r
verify_pmf_marginal_properties(
  J,
  a,
  b,
  logS,
  M = .QUAD_NODES_DEFAULT,
  tol = 1e-10,
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

- M:

  Integer; number of quadrature nodes.

- tol:

  Numeric; tolerance for comparisons.

- verbose:

  Logical; if TRUE, print results.

## Value

Logical; TRUE if all verifications pass.

## Examples

``` r
logS <- compute_log_stirling(50)
verify_pmf_marginal_properties(50, 1.5, 0.5, logS)
#> Marginal PMF Properties (J=50, a=1.50, b=0.50):
#>   Sum = 1:          PASS (sum = 1.000000000000)
#>   P(K=0) = 0:       PASS (P(K=0) = 0.00e+00)
#>   Non-negative:     PASS (min = 0.00e+00)
#>   CDF monotonic:    PASS
#>   CDF[J] = 1:       PASS (CDF[J] = 1.000000000000)
#>   Overall: PASS
```

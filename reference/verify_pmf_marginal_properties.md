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
if (FALSE) { # \dontrun{
logS <- compute_log_stirling(50)
verify_pmf_marginal_properties(50, 1.5, 0.5, logS)

} # }
```

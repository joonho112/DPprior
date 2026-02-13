# Verify Moments Consistency

Verifies that moments computed from the marginal PMF match those from
[`exact_K_moments()`](https://joonho112.github.io/DPprior/reference/exact_K_moments.md).

## Usage

``` r
verify_pmf_marginal_moments(
  J,
  a,
  b,
  logS,
  M = .QUAD_NODES_DEFAULT,
  tol = 1e-06,
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

Logical; TRUE if moments match within tolerance.

## Examples

``` r
if (FALSE) { # \dontrun{
logS <- compute_log_stirling(50)
verify_pmf_marginal_moments(50, 1.5, 0.5, logS)

} # }
```

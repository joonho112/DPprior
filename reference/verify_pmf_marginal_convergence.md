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
if (FALSE) { # \dontrun{
logS <- compute_log_stirling(50)
verify_pmf_marginal_convergence(50, 1.5, 0.5, logS)

} # }
```

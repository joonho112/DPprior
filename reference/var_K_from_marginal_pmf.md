# Variance of Marginal K from PMF

Computes \\Var(K_J \mid a, b)\\ from the marginal PMF.

## Usage

``` r
var_K_from_marginal_pmf(J, a, b, logS, M = .QUAD_NODES_DEFAULT)
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

## Value

Numeric; marginal variance.

## Examples

``` r
if (FALSE) { # \dontrun{
logS <- compute_log_stirling(50)
var_K_from_marginal_pmf(50, 1.5, 0.5, logS)

} # }
```

# Variance of K from PMF

Computes \\Var(K_J \mid \alpha)\\ by summing over the PMF. This is
primarily for verification against the closed-form trigamma formula.

## Usage

``` r
var_K_from_pmf(J, alpha, logS)
```

## Arguments

- J:

  Integer; sample size.

- alpha:

  Numeric; DP concentration parameter.

- logS:

  Matrix; pre-computed log-Stirling matrix.

## Value

Numeric; the conditional variance \\Var(K_J \mid \alpha)\\.

## Examples

``` r
logS <- compute_log_stirling(50)

# Should match var_K_given_alpha(50, 2.0)
var_K_from_pmf(50, 2.0, logS)
#> [1] 4.535558
var_K_given_alpha(50, 2.0)
#> [1] 4.535558
```

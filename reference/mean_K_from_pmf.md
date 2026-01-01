# Mean of K from PMF

Computes \\E\[K_J \mid \alpha\]\\ by summing over the PMF. This is
primarily for verification against the closed-form digamma formula.

## Usage

``` r
mean_K_from_pmf(J, alpha, logS)
```

## Arguments

- J:

  Integer; sample size.

- alpha:

  Numeric; DP concentration parameter.

- logS:

  Matrix; pre-computed log-Stirling matrix.

## Value

Numeric; the conditional mean \\E\[K_J \mid \alpha\]\\.

## Examples

``` r
logS <- compute_log_stirling(50)

# Should match mean_K_given_alpha(50, 2.0)
mean_K_from_pmf(50, 2.0, logS)
#> [1] 7.037626
mean_K_given_alpha(50, 2.0)
#> [1] 7.037626
```

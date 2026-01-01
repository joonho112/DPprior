# Mode of K Given Alpha

Computes the mode (most likely value) of \\K_J \mid \alpha\\.

## Usage

``` r
mode_K_given_alpha(J, alpha, logS)
```

## Arguments

- J:

  Integer; sample size.

- alpha:

  Numeric; DP concentration parameter.

- logS:

  Matrix; pre-computed log-Stirling matrix.

## Value

Integer; the value \\k\\ that maximizes \\P(K_J = k \mid \alpha)\\.

## Examples

``` r
logS <- compute_log_stirling(50)
mode_K_given_alpha(50, 2.0, logS)
#> [1] 7
```

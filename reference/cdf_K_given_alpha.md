# CDF of K Given Alpha

Computes the cumulative distribution function \\P(K_J \leq k \mid
\alpha)\\ for \\k = 0, 1, \ldots, J\\.

## Usage

``` r
cdf_K_given_alpha(J, alpha, logS)
```

## Arguments

- J:

  Integer; sample size.

- alpha:

  Numeric; DP concentration parameter.

- logS:

  Matrix; pre-computed log-Stirling matrix.

## Value

Numeric vector of length \\J+1\\ containing \\P(K_J \leq k \mid
\alpha)\\ for \\k = 0, 1, \ldots, J\\.

## Details

The CDF satisfies:

- \\F(0) = 0\\ (since \\P(K_J = 0) = 0\\)

- \\F(J) = 1\\

- \\F(k)\\ is non-decreasing in \\k\\

## Examples

``` r
logS <- compute_log_stirling(50)
cdf <- cdf_K_given_alpha(50, 2.0, logS)

# Verify CDF ends at 1
cdf[51]  # Should be 1
#> [1] 1

# P(K <= 5)
cdf[6]
#> [1] 0.2418684
```

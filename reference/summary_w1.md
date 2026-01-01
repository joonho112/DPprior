# Summary Statistics for w₁ Distribution

Computes comprehensive summary statistics for the w₁ distribution.

## Usage

``` r
summary_w1(
  a,
  b,
  probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
  M = .QUAD_NODES_DEFAULT
)
```

## Arguments

- a:

  Numeric; shape parameter of the Gamma prior on α (a \> 0).

- b:

  Numeric; rate parameter of the Gamma prior on α (b \> 0).

- probs:

  Numeric vector; quantile probabilities. Default is
  `c(0.05, 0.25, 0.5, 0.75, 0.95)`.

- M:

  Integer; number of quadrature nodes for mean/variance. Default is 80.

## Value

A list of class "w1_summary" containing:

- mean:

  E(w₁)

- var:

  Var(w₁)

- sd:

  SD(w₁) = sqrt(Var(w₁))

- median:

  Median of w₁

- quantiles:

  Named vector of quantiles

- prob_gt_50:

  P(w₁ \> 0.5), dominance indicator

- prob_gt_90:

  P(w₁ \> 0.9), extreme dominance indicator

- params:

  List of input parameters (a, b)

## See also

[`cdf_w1`](https://joonho112.github.io/DPprior/reference/cdf_w1.md),
[`quantile_w1`](https://joonho112.github.io/DPprior/reference/quantile_w1.md),
[`mean_w1`](https://joonho112.github.io/DPprior/reference/mean_w1.md)

## Examples

``` r
# Standard summary
summary_w1(a = 2, b = 1)
#> w1 Distribution Summary
#> ============================================= 
#> 
#> Gamma prior: alpha ~ Gamma(2.0000, 1.0000)
#> E[alpha] = 2.0000, CV(alpha) = 70.71%
#> 
#> Location and Scale:
#> ------------------------------ 
#>   Mean:   0.4037
#>   Median: 0.3391
#>   SD:     0.2995
#> 
#> Quantiles:
#> ------------------------------ 
#>   q5: 0.0256
#>   q25: 0.1433
#>   q50: 0.3391
#>   q75: 0.6321
#>   q95: 0.9689 
#> 
#> Dominance Risk:
#> ------------------------------ 
#>   P(w1 > 0.5): 0.3488
#>   P(w1 > 0.9): 0.0917

# Lee et al. DP-inform prior
summary_w1(a = 1.6, b = 1.22)
#> w1 Distribution Summary
#> ============================================= 
#> 
#> Gamma prior: alpha ~ Gamma(1.6000, 1.2200)
#> E[alpha] = 1.3115, CV(alpha) = 79.06%
#> 
#> Location and Scale:
#> ------------------------------ 
#>   Mean:   0.5084
#>   Median: 0.4839
#>   SD:     0.3244
#> 
#> Quantiles:
#> ------------------------------ 
#>   q5: 0.0390
#>   q25: 0.2136
#>   q50: 0.4839
#>   q75: 0.8139
#>   q95: 0.9988 
#> 
#> Dominance Risk:
#> ------------------------------ 
#>   P(w1 > 0.5): 0.4868
#>   P(w1 > 0.9): 0.1833
```

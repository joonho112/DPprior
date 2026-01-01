# Summary Statistics for rho Distribution

Computes comprehensive summary statistics for the co-clustering
probability rho under the hierarchical prior alpha ~ Gamma(a, b).

## Usage

``` r
summary_rho(a, b, M = .QUAD_NODES_DEFAULT)
```

## Arguments

- a:

  Numeric; shape parameter of the Gamma prior on alpha (a \> 0).

- b:

  Numeric; rate parameter of the Gamma prior on alpha (b \> 0).

- M:

  Integer; number of quadrature nodes. Default is 80.

## Value

A list of class "rho_summary" containing:

- mean:

  E(rho \| a, b)

- var:

  Var(rho \| a, b)

- sd:

  SD(rho \| a, b) = sqrt(Var)

- cv:

  Coefficient of variation SD/mean

- interpretation:

  Qualitative interpretation of co-clustering level

- params:

  List of input parameters (a, b)

- alpha_prior:

  Summary of the alpha prior (mean, sd, cv)

- conditional_at_alpha_mean:

  Conditional moments evaluated at E(alpha)

## Details

The co-clustering probability rho indicates how likely two randomly
chosen observations are to belong to the same cluster a priori.

**Interpretation guidelines:**

- E(rho) \> 0.5: High co-clustering; most pairs expected in same cluster

- E(rho) in (0.2, 0.5): Moderate co-clustering

- E(rho) \< 0.2: Low co-clustering; fragmented prior

The `conditional_at_alpha_mean` component provides a "plug-in" estimate
for comparison: what the moments would be if alpha were fixed at its
prior mean.

## See also

[`mean_rho`](https://joonho112.github.io/DPprior/reference/mean_rho.md),
[`var_rho`](https://joonho112.github.io/DPprior/reference/var_rho.md),
[`summary_w1`](https://joonho112.github.io/DPprior/reference/summary_w1.md)

## Examples

``` r
summary_rho(a = 2, b = 1)
#> Co-Clustering Probability (rho) Summary
#> ================================================== 
#> 
#> Gamma prior: alpha ~ Gamma(2.0000, 1.0000)
#> E[alpha] = 2.0000, SD(alpha) = 1.4142, CV(alpha) = 70.7%
#> 
#> Marginal distribution of rho:
#> ----------------------------------- 
#>   Mean:   0.4037
#>   SD:     0.2397
#>   CV:     59.4%
#> 
#> Conditional at E[alpha] = 2.0000 (plug-in):
#> ----------------------------------- 
#>   E[rho | E[alpha]]:   0.3333
#>   Var(rho | E[alpha]): 0.0222
#> 
#> Interpretation:
#> ----------------------------------- 
#>   Moderate co-clustering
summary_rho(a = 1.6, b = 1.22)
#> Co-Clustering Probability (rho) Summary
#> ================================================== 
#> 
#> Gamma prior: alpha ~ Gamma(1.6000, 1.2200)
#> E[alpha] = 1.3115, SD(alpha) = 1.0368, CV(alpha) = 79.1%
#> 
#> Marginal distribution of rho:
#> ----------------------------------- 
#>   Mean:   0.5084
#>   SD:     0.2664
#>   CV:     52.4%
#> 
#> Conditional at E[alpha] = 1.3115 (plug-in):
#> ----------------------------------- 
#>   E[rho | E[alpha]]:   0.4326
#>   Var(rho | E[alpha]): 0.0344
#> 
#> Interpretation:
#> ----------------------------------- 
#>   High co-clustering: most pairs expected in same cluster
```

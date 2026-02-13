# Weight Distribution Diagnostics (w1)

Computes comprehensive diagnostics for the first stick-breaking weight.
This is the KEY diagnostic for concerns about unintended prior behavior
(Lee, 2026, Section 4).

## Usage

``` r
compute_weight_diagnostics(
  a,
  b,
  thresholds = c(0.3, 0.5, 0.7, 0.9),
  M = .QUAD_NODES_DEFAULT
)
```

## Arguments

- a:

  Numeric; shape parameter of the Gamma prior (\> 0).

- b:

  Numeric; rate parameter of the Gamma prior (\> 0).

- thresholds:

  Numeric vector; thresholds for computing P(w1 \> t). Default is c(0.3,
  0.5, 0.7, 0.9). All values must be in (0, 1).

- M:

  Integer; number of quadrature nodes for moment computation (default:
  80).

## Value

A list with components: mean (expected value of the first weight),
median, quantiles (named vector), prob_exceeds (named vector of
exceedance probabilities for each threshold), and dominance_risk
(character: "low", "moderate", or "high").

## Details

The weight w1 is the first stick-breaking weight in GEM order
(size-biased permutation), representing the asymptotic cluster share of
a randomly chosen unit.

Dominance risk classification:

- "low": P(w1 \> 0.5) \< 0.2

- "moderate": 0.2 \<= P(w1 \> 0.5) \< 0.4

- "high": P(w1 \> 0.5) \>= 0.4

## References

Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet
Process Mixtures. *arXiv preprint* arXiv:2602.06301.

## See also

Other diagnostics:
[`DPprior_diagnostics()`](https://joonho112.github.io/DPprior/reference/DPprior_diagnostics.md),
[`DPprior_error_bounds()`](https://joonho112.github.io/DPprior/reference/DPprior_error_bounds.md),
[`compute_error_landscape()`](https://joonho112.github.io/DPprior/reference/compute_error_landscape.md),
[`dual_anchor_diagnostics()`](https://joonho112.github.io/DPprior/reference/dual_anchor_diagnostics.md)

## Examples

``` r
# Lee et al. DP-inform prior (known high dominance)
compute_weight_diagnostics(1.60, 1.22)
#> $mean
#> [1] 0.508368
#> 
#> $median
#> [1] 0.4839219
#> 
#> $quantiles
#>         q5        q25        q50        q75        q95 
#> 0.03896534 0.21361987 0.48392192 0.81393615 0.99878645 
#> 
#> $prob_exceeds
#> prob_gt_0.3 prob_gt_0.5 prob_gt_0.7 prob_gt_0.9 
#>   0.6634195   0.4868311   0.3333737   0.1833147 
#> 
#> $dominance_risk
#> [1] "high"
#> 

# Lower dominance case
compute_weight_diagnostics(5, 1)
#> $mean
#> [1] 0.1915145
#> 
#> $median
#> [1] 0.138171
#> 
#> $quantiles
#>         q5        q25        q50        q75        q95 
#> 0.01025848 0.05750422 0.13817096 0.27349354 0.55981677 
#> 
#> $prob_exceeds
#> prob_gt_0.3 prob_gt_0.5 prob_gt_0.7 prob_gt_0.9 
#> 0.217581005 0.071866491 0.019229538 0.002545247 
#> 
#> $dominance_risk
#> [1] "low"
#> 
```

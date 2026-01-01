# Get Quadrature Information Summary

Returns summary information about the quadrature nodes and weights for a
given Gamma(a, b) distribution.

## Usage

``` r
summary_quadrature(a, b, M = .QUAD_NODES_DEFAULT)
```

## Arguments

- a:

  Numeric; shape parameter of Gamma distribution.

- b:

  Numeric; rate parameter of Gamma distribution.

- M:

  Integer; number of quadrature nodes.

## Value

A list with components:

- `n_nodes`:

  Number of quadrature nodes.

- `alpha_range`:

  Range of alpha nodes (min, max).

- `weight_range`:

  Range of normalized weights (min, max).

- `gamma_mean`:

  Theoretical mean of Gamma(a, b).

- `gamma_sd`:

  Theoretical SD of Gamma(a, b).

- `coverage`:

  Approximate coverage in terms of SD from mean.

## Examples

``` r
summary_quadrature(2.5, 1.5, M = 80)
#> $n_nodes
#> [1] 80
#> 
#> $alpha_range
#>          min          max 
#>   0.04141985 199.91970270 
#> 
#> $weight_range
#>       min       max 
#> 0.0000000 0.1413702 
#> 
#> $gamma_mean
#> [1] 1.666667
#> 
#> $gamma_sd
#> [1] 1.054093
#> 
#> $coverage
#>   lower_sd   upper_sd 
#>   1.541845 188.079344 
#> 
```

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
if (FALSE) { # \dontrun{
summary_quadrature(2.5, 1.5, M = 80)

} # }
```

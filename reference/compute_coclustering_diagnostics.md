# Co-Clustering Diagnostics (rho)

Computes diagnostics for rho = sum of w_h squared, the prior
co-clustering probability.

## Usage

``` r
compute_coclustering_diagnostics(a, b, M = .QUAD_NODES_DEFAULT)
```

## Arguments

- a:

  Numeric; shape parameter of the Gamma prior (\> 0).

- b:

  Numeric; rate parameter of the Gamma prior (\> 0).

- M:

  Integer; number of quadrature nodes (default: 80).

## Value

A list with components: mean (expected co-clustering probability), var
(variance), sd, and interpretation (qualitative description of
co-clustering level).

## Details

Key identity from RN-06: the expected co-clustering probability equals
the expected first stick-breaking weight.

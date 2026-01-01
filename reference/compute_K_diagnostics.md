# K Distribution Diagnostics

Computes summary statistics for the marginal distribution of K_J under
alpha distributed as Gamma(a, b).

## Usage

``` r
compute_K_diagnostics(J, a, b, M = .QUAD_NODES_DEFAULT)
```

## Arguments

- J:

  Integer; sample size (positive integer \>= 1).

- a:

  Numeric; shape parameter of the Gamma prior (\> 0).

- b:

  Numeric; rate parameter of the Gamma prior (\> 0).

- M:

  Integer; number of Gauss-Laguerre quadrature nodes (default: 80).

## Value

A list with components: mean (expected value of K), var (variance), sd,
mode, median, quantiles (named integer vector), and pmf (full PMF vector
for k = 1, ..., J).

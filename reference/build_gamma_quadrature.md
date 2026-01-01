# Build Quadrature for Gamma(a, b) Integration

Transforms standard Gauss-Laguerre quadrature for integration against a
Gamma(a, b) distribution with shape `a` and rate `b`.

## Usage

``` r
build_gamma_quadrature(a, b, M = .QUAD_NODES_DEFAULT)
```

## Arguments

- a:

  Numeric; shape parameter of Gamma distribution (must be \> 0).

- b:

  Numeric; rate parameter of Gamma distribution (must be \> 0).

- M:

  Integer; number of quadrature nodes (default: 80).

## Value

A list with components:

- `a`:

  Shape parameter.

- `b`:

  Rate parameter.

- `alpha_nodes`:

  Numeric vector of transformed nodes \\\alpha_m = x_m / b\\ on the
  \\\alpha\\ scale.

- `weights_normalized`:

  Numeric vector of normalized weights that sum to 1.

## Details

For \\\alpha \sim \text{Gamma}(a, b)\\: \$\$E\[g(\alpha)\] =
\frac{1}{\Gamma(a)} \int_0^\infty g(x/b) x^{a-1} e^{-x} dx\$\$

Using generalized Laguerre quadrature with parameter
\\\alpha\_{\text{param}} = a - 1\\: \$\$E\[g(\alpha)\] \approx
\sum\_{m=1}^M \tilde{w}\_m g(\alpha_m)\$\$

where \\\alpha_m = x_m / b\\ and \\\tilde{w}\_m\\ are normalized weights
summing to 1.

The weights are normalized in log-space for numerical stability, which
ensures monotone convergence as M increases.

## See also

[`gauss_laguerre_nodes`](https://joonho112.github.io/DPprior/reference/gauss_laguerre_nodes.md)
for raw quadrature computation,
[`integrate_gamma`](https://joonho112.github.io/DPprior/reference/integrate_gamma.md)
for high-level expectation computation

## Examples

``` r
# Build quadrature for Gamma(2.5, 1.5)
quad <- build_gamma_quadrature(2.5, 1.5)

# Check weights sum to 1
sum(quad$weights_normalized)
#> [1] 1

# Check mean: E[alpha] should be a/b = 2.5/1.5
sum(quad$weights_normalized * quad$alpha_nodes)
#> [1] 1.666667
```

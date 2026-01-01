# Integrate Function Against Gamma Distribution

High-level interface for computing \\E\[f(\alpha)\]\\ where \\\alpha
\sim \text{Gamma}(a, b)\\.

## Usage

``` r
integrate_gamma(f, a, b, M = .QUAD_NODES_DEFAULT)
```

## Arguments

- f:

  Function to integrate; must accept a single numeric argument and
  return a single numeric value.

- a:

  Numeric; shape parameter of Gamma distribution (must be \> 0).

- b:

  Numeric; rate parameter of Gamma distribution (must be \> 0).

- M:

  Integer; number of quadrature nodes (default: 80).

## Value

Numeric; approximation of \\E\[f(\alpha)\]\\.

## Details

Uses Gauss-Laguerre quadrature to approximate: \$\$E\[f(\alpha)\] =
\int_0^\infty f(\alpha) \frac{b^a}{\Gamma(a)} \alpha^{a-1} e^{-b\alpha}
d\alpha\$\$

The approximation is: \$\$E\[f(\alpha)\] \approx \sum\_{m=1}^M w_m
f(\alpha_m)\$\$

where \\\alpha_m\\ are quadrature nodes and \\w_m\\ are normalized
weights.

For polynomial integrands of degree up to \\2M - 1\\, the quadrature is
exact. For other smooth functions, accuracy improves rapidly with \\M\\.

## See also

[`build_gamma_quadrature`](https://joonho112.github.io/DPprior/reference/build_gamma_quadrature.md)
for the underlying quadrature

## Examples

``` r
# E[alpha] for Gamma(2.5, 1.5) should be 2.5/1.5
integrate_gamma(identity, 2.5, 1.5)
#> [1] 1.666667

# E[alpha^2] for Gamma(2.5, 1.5) should be 2.5*3.5/1.5^2
integrate_gamma(function(x) x^2, 2.5, 1.5)
#> [1] 3.888889

# More complex function
integrate_gamma(function(x) log(x + 1), 2.5, 1.5)
#> [1] 0.9102673
```

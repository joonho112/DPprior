# Density of w₁

Computes the probability density p(w₁ = x \| a, b).

## Usage

``` r
density_w1(x, a, b, log = FALSE)
```

## Arguments

- x:

  Numeric vector; evaluation points.

- a:

  Numeric; shape parameter of the Gamma prior on α (a \> 0).

- b:

  Numeric; rate parameter of the Gamma prior on α (b \> 0).

- log:

  Logical; if `TRUE`, returns log-density. Default is `FALSE`.

## Value

Numeric vector of density (or log-density if `log = TRUE`) values.
Returns 0 (or -Inf on log scale) for x outside (0, 1).

## Details

The marginal density of w₁ is: \$\$p(w_1 \| a, b) = \frac{a \cdot
b^a}{(1-w_1) \cdot \[b - \log(1-w_1)\]^{a+1}}\$\$

**Important:** For small values of a (a \< 1), the density has
significant mass concentrated very close to x = 1.

## See also

[`cdf_w1`](https://joonho112.github.io/DPprior/reference/cdf_w1.md),
[`quantile_w1`](https://joonho112.github.io/DPprior/reference/quantile_w1.md)

## Examples

``` r
# Density at several points
x <- seq(0.1, 0.9, by = 0.1)
density_w1(x, a = 2, b = 1)
#> [1] 1.6454157 1.3661794 1.1442068 0.9665754 0.8240923 0.7105355 0.6227160
#> [8] 0.5628065 0.5552236

# Log-density for numerical stability
density_w1(0.5, a = 2, b = 1, log = TRUE)
#> [1] -0.1934727
```

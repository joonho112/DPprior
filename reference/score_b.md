# Score Function with Respect to Rate Parameter b

Computes the score function \\s_b(\alpha) = \partial/\partial b \log
g\_{a,b}(\alpha)\\ for the Gamma(a, b) distribution.

## Usage

``` r
score_b(alpha, a, b)
```

## Arguments

- alpha:

  Numeric vector; points at which to evaluate.

- a:

  Numeric scalar; shape parameter of the Gamma distribution (\> 0).

- b:

  Numeric scalar; rate parameter of the Gamma distribution (\> 0).

## Value

Numeric vector of the same length as `alpha`.

## Details

For the Gamma(shape = a, rate = b) distribution, the score function with
respect to `b` is: \$\$s_b(\alpha) = \frac{a}{b} - \alpha.\$\$

A fundamental property of score functions is that their expectation is
zero: \$\$E\_{\alpha \sim g\_{a,b}}\[s_b(\alpha)\] = 0.\$\$

Unlike `score_a`, this function is linear in \\\alpha\\, so its
expectation converges very quickly with quadrature.

## References

RN-04, Section 4.2: Jacobian via score identities

## See also

[`score_a`](https://joonho112.github.io/DPprior/reference/score_a.md)
for the score with respect to a

## Examples

``` r
# Evaluate score at several points
alpha_vals <- c(0.5, 1.0, 2.0, 5.0)
score_b(alpha_vals, a = 2.0, b = 1.0)
#> [1]  1.5  1.0  0.0 -3.0
```

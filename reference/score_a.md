# Score Function with Respect to Shape Parameter a

Computes the score function \\s_a(\alpha) = \partial/\partial a \log
g\_{a,b}(\alpha)\\ for the Gamma(a, b) distribution.

## Usage

``` r
score_a(alpha, a, b)
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

For the Gamma(shape = a, rate = b) distribution with density
\$\$g\_{a,b}(\alpha) = \frac{b^a}{\Gamma(a)} \alpha^{a-1}
e^{-b\alpha},\$\$ the score function with respect to `a` is:
\$\$s_a(\alpha) = \log b - \psi(a) + \log \alpha,\$\$ where \\\psi\\ is
the digamma function.

A fundamental property of score functions is that their expectation is
zero: \$\$E\_{\alpha \sim g\_{a,b}}\[s_a(\alpha)\] = 0.\$\$

**Numerical Note:** The `log(alpha)` term causes slower quadrature
convergence for score-weighted integrands compared to the moments
themselves. Use higher M (e.g., 120-200) for verification purposes.

## References

Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet
Process Mixtures. *arXiv preprint* arXiv:2602.06301.

## See also

[`score_b`](https://joonho112.github.io/DPprior/reference/score_b.md)
for the score with respect to b

## Examples

``` r
# Evaluate score at several points
alpha_vals <- c(0.5, 1.0, 2.0, 5.0)
score_a(alpha_vals, a = 2.0, b = 1.0)
#> [1] -1.1159315 -0.4227843  0.2703628  1.1866536
```

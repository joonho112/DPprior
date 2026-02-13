# Marginal Variance of rho

Computes Var(rho \| a, b) when alpha ~ Gamma(a, b) (shape-rate).

## Usage

``` r
var_rho(a, b, M = .QUAD_NODES_DEFAULT)
```

## Arguments

- a:

  Numeric; shape parameter of the Gamma prior on alpha (a \> 0).

- b:

  Numeric; rate parameter of the Gamma prior on alpha (b \> 0).

- M:

  Integer; number of quadrature nodes. Default is 80.

## Value

Numeric; Var(rho \| a, b).

## Details

Uses the law of total variance: \$\$Var(\rho \| a, b) = E\[Var(\rho \|
\alpha)\] + Var(E\[\rho \| \alpha\])\$\$

where:

- `Var(rho | alpha) = 2*alpha / ((1+alpha)^2*(2+alpha)*(3+alpha))`

- `E(rho | alpha) = 1/(1+alpha)`

**Note:** Unlike E(rho), Var(rho) != Var(w1) in general, because the
conditional variances differ.

## References

Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet
Process Mixtures. *arXiv preprint* arXiv:2602.06301.

## See also

[`mean_rho`](https://joonho112.github.io/DPprior/reference/mean_rho.md),
[`cv_rho`](https://joonho112.github.io/DPprior/reference/cv_rho.md),
[`var_w1`](https://joonho112.github.io/DPprior/reference/var_w1.md)

Other co_clustering:
[`cv_rho()`](https://joonho112.github.io/DPprior/reference/cv_rho.md),
[`mean_rho()`](https://joonho112.github.io/DPprior/reference/mean_rho.md),
[`mean_rho_given_alpha()`](https://joonho112.github.io/DPprior/reference/mean_rho_given_alpha.md),
[`var_rho_given_alpha()`](https://joonho112.github.io/DPprior/reference/var_rho_given_alpha.md)

## Examples

``` r
var_rho(a = 2, b = 1)
#> [1] 0.05744825
```

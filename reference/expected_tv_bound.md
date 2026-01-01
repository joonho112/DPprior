# Expected TV Bound Under Gamma Prior

Integrates the conditional TV bound over \\\alpha \sim \text{Gamma}(a,
b)\\ to obtain the marginal error bound.

## Usage

``` r
expected_tv_bound(J, a, b, cJ = log(J), M = .QUAD_NODES_DEFAULT)
```

## Arguments

- J:

  Integer; sample size.

- a, b:

  Numeric; Gamma hyperparameters.

- cJ:

  Numeric; scaling constant (default: log(J)).

- M:

  Integer; number of quadrature nodes.

## Value

Numeric; \\E\[d\_{TV} \text{ bound} \| a, b\]\\.

## Details

From Corollary 1 of RN-05, the TV error between the exact prior
predictive \\p(S_J \| a, b)\\ and the A1 shifted NegBin proxy is bounded
by: \$\$d\_{TV}(P^{\text{exact}}, Q^{A1}) \le E\_{\alpha \sim
\Gamma(a,b)}\[B\_{\text{Pois}} + B\_{\text{lin}}\]\$\$

This follows from the mixture contraction property of TV distance.

## See also

[`compute_total_tv_bound`](https://joonho112.github.io/DPprior/reference/compute_total_tv_bound.md),
[`integrate_gamma`](https://joonho112.github.io/DPprior/reference/integrate_gamma.md)

## Examples

``` r
# Marginal TV bound for J=100, Gamma(1, 1)
expected_tv_bound(J = 100, a = 1, b = 1)
#> [1] 0.2465068
```

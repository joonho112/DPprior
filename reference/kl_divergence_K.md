# KL Divergence Between Target and Induced K_J PMFs

Computes the KL divergence \\D\_{KL}(p^\* \\ p\_{a,b})\\ between a
target PMF and the induced marginal PMF of \\K_J\\ under \\\alpha \sim
Gamma(a, b)\\.

## Usage

``` r
kl_divergence_K(target_pmf, a, b, J, M = .QUAD_NODES_DEFAULT)
```

## Arguments

- target_pmf:

  Numeric vector; target PMF for K_J. Can have length J (support
  k=1,...,J) or J+1 (support k=0,...,J, where k=0 is ignored).

- a:

  Numeric; shape parameter of Gamma hyperprior (a \> 0).

- b:

  Numeric; rate parameter of Gamma hyperprior (b \> 0).

- J:

  Integer; sample size.

- M:

  Integer; number of quadrature nodes. Default: 80.

## Value

Numeric scalar; \\D\_{KL}(p^\* \\ p\_{a,b})\\ (non-negative).

## See also

[`kl_divergence_pmf`](https://joonho112.github.io/DPprior/reference/kl_divergence_pmf.md),
[`DPprior_a2_kl`](https://joonho112.github.io/DPprior/reference/DPprior_a2_kl.md)

## Examples

``` r
J <- 50
target <- rep(1/J, J)  # uniform target over k=1,...,J
kl_divergence_K(target, a = 2, b = 1, J = J)
#> [1] 10.49288
```

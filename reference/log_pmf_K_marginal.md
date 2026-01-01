# Log Marginal PMF of K_J under Gamma Hyperprior

Computes \\\log P(K_J = k \mid a, b)\\ for \\k = 0, 1, \ldots, J\\ using
log-space mixing for numerical stability.

## Usage

``` r
log_pmf_K_marginal(J, a, b, logS, M = .QUAD_NODES_DEFAULT)
```

## Arguments

- J:

  Integer; sample size (must be \>= 1).

- a:

  Numeric; shape parameter of Gamma prior (must be \> 0).

- b:

  Numeric; rate parameter of Gamma prior (must be \> 0).

- logS:

  Matrix; pre-computed log-Stirling matrix from
  [`compute_log_stirling`](https://joonho112.github.io/DPprior/reference/compute_log_stirling.md).

- M:

  Integer; number of quadrature nodes (default: `.QUAD_NODES_DEFAULT`).

## Value

Numeric vector of length \\J+1\\ containing log-probabilities for \\k =
0, 1, \ldots, J\\. Entry `[1]` corresponds to \\k=0\\ and is always
`-Inf`.

## Details

This routine normalizes \\P(K_J = \cdot \mid \alpha_m)\\ at each
quadrature node before mixing, and then mixes in log-space via
`logsumexp_vec`: \$\$\log p_k \approx \log\sum_m \exp\\\log w_m + \log
p\_{k\mid m}\\.\$\$

The log-space computation is essential for numerical stability when:

- J is large (tail probabilities become very small)

- Alpha values span a wide range (extreme quadrature nodes)

- Parameters lead to concentrated distributions

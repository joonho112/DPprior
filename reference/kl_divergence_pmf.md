# KL Divergence Between Two PMFs

Computes the Kullback-Leibler divergence \\D\_{KL}(p \\ q)\\ between two
probability mass functions.

## Usage

``` r
kl_divergence_pmf(p, q, eps = 1e-15)
```

## Arguments

- p:

  Numeric vector; target PMF (reference distribution).

- q:

  Numeric vector; comparison PMF.

- eps:

  Numeric; small value to prevent `log(0)`. Default: 1e-15.

## Value

Numeric scalar; the KL divergence (non-negative).

## Details

The KL divergence is defined as: \$\$D\_{KL}(p \\ q) = \sum_k p(k)
\log\frac{p(k)}{q(k)}\$\$

Only indices where \\p(k) \> \epsilon\\ are included in the sum.

**Properties:**

- \\D\_{KL}(p \\ q) \geq 0\\ with equality iff \\p = q\\

- Not symmetric: \\D\_{KL}(p \\ q) \neq D\_{KL}(q \\ p)\\

## See also

[`kl_divergence_K`](https://joonho112.github.io/DPprior/reference/kl_divergence_K.md)
for KL divergence with induced PMF

## Examples

``` r
p <- c(0.2, 0.5, 0.3)
kl_divergence_pmf(p, p)  # 0
#> [1] 0

q <- c(0.3, 0.4, 0.3)
kl_divergence_pmf(p, q)
#> [1] 0.03047875
```

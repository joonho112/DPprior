# Mean-Linearization Error Bound

Computes an upper bound on the TV distance between two Poisson
distributions, \\\text{Poisson}(\lambda_J(\alpha))\\ and
\\\text{Poisson}(\alpha c_J)\\, using the Poisson-Poisson KL divergence
together with Pinsker's inequality.

## Usage

``` r
compute_linearization_bound(J, alpha, cJ = log(J))
```

## Arguments

- J:

  Integer; sample size.

- alpha:

  Numeric; concentration parameter (vectorized).

- cJ:

  Numeric; scaling constant (default: log(J)).

## Value

Numeric; upper bound via Pinsker's inequality.

## Details

Let \\\lambda = \lambda_J(\alpha)\\ (exact shifted mean) and \\\lambda'
= \alpha c_J\\ (A1 approximate mean).

The KL divergence is: \$\$KL(\text{Poisson}(\lambda) \|\|
\text{Poisson}(\lambda')) = \lambda \log(\lambda/\lambda') + \lambda' -
\lambda\$\$

By Pinsker's inequality: \$\$d\_{TV}(\text{Poisson}(\lambda),
\text{Poisson}(\lambda')) \le \sqrt{KL/2}\$\$

Numerical safeguards handle edge cases where \\\lambda\\ or \\c_J\\ is
zero.

## References

RN-05, Lemma 1: Poisson-Poisson TV bound via KL + Pinsker

## See also

[`compute_poissonization_bound`](https://joonho112.github.io/DPprior/reference/compute_poissonization_bound.md),
[`compute_total_tv_bound`](https://joonho112.github.io/DPprior/reference/compute_total_tv_bound.md)

## Examples

``` r
# Linearization bound for J=50, alpha=1
compute_linearization_bound(J = 50, alpha = 1)
#> [1] 0.1062796

# Effect of J on linearization bound (should decrease)
sapply(c(25, 50, 100, 200), function(J)
  compute_linearization_bound(J, alpha = 2))
#> [1] 0.3579663 0.3328078 0.3098929 0.2899216
```

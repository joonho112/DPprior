# Total TV Error Bound (Conditional)

Computes the combined conditional TV bound using the triangle
inequality: \$\$d\_{TV}(K_J \| \alpha, 1 + \text{NegBin}) \le
B\_{\text{Pois}} + B\_{\text{lin}}\$\$

## Usage

``` r
compute_total_tv_bound(J, alpha, cJ = log(J))
```

## Arguments

- J:

  Integer; sample size.

- alpha:

  Numeric; concentration parameter (vectorized).

- cJ:

  Numeric; scaling constant (default: log(J)).

## Value

Numeric; upper bound on total TV error (capped at 1).

## Details

The result is capped at 1 since TV distance is bounded by 1.

From RN-05, Theorem 1, the total conditional TV error decomposes as:

1.  Poissonization error: \\S_J \| \alpha\\ vs
    \\\text{Poisson}(\lambda_J(\alpha))\\

2.  Linearization error: \\\text{Poisson}(\lambda_J(\alpha))\\ vs
    \\\text{Poisson}(\alpha c_J)\\

## See also

[`compute_poissonization_bound`](https://joonho112.github.io/DPprior/reference/compute_poissonization_bound.md),
[`compute_linearization_bound`](https://joonho112.github.io/DPprior/reference/compute_linearization_bound.md),
[`expected_tv_bound`](https://joonho112.github.io/DPprior/reference/expected_tv_bound.md)

## Examples

``` r
# Total bound at alpha = E[alpha] under Gamma(2, 1)
compute_total_tv_bound(J = 50, alpha = 2)
#> [1] 0.5809985

# Vectorized
compute_total_tv_bound(J = 50, alpha = c(0.5, 1, 2, 5))
#> [1] 0.1075542 0.2795305 0.5809985 1.0000000
```

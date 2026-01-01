# Conditional Mean of rho Given Alpha

Computes the conditional mean E(rho \| alpha) for the co-clustering
probability rho = sum_h w_h^2 under a Dirichlet Process.

## Usage

``` r
mean_rho_given_alpha(alpha)
```

## Arguments

- alpha:

  Numeric vector; concentration parameter(s) (must be positive).

## Value

Numeric vector; `E(rho | alpha) = 1/(1+alpha)`.

## Details

The co-clustering probability rho = sum(w_h^2) over h \>= 1 has
conditional mean: \$\$E\[\rho \| \alpha\] = \frac{1}{1 + \alpha}\$\$

This equals `E(w1 | alpha)` since w1 ~ Beta(1, alpha) has mean
1/(1+alpha).

**Interpretation:**

- alpha -\> 0: E(rho\|alpha) -\> 1 (all observations in one cluster)

- alpha -\> Inf: E(rho\|alpha) -\> 0 (infinitely many small clusters)

- alpha = 1: E(rho\|alpha) = 0.5 (moderate clustering)

## References

Lee, J. (2025). RN-06: Dual-Anchor Design II, Section 3.2.

## See also

[`var_rho_given_alpha`](https://joonho112.github.io/DPprior/reference/var_rho_given_alpha.md),
[`mean_rho`](https://joonho112.github.io/DPprior/reference/mean_rho.md)

## Examples

``` r
mean_rho_given_alpha(1.0)
#> [1] 0.5
mean_rho_given_alpha(c(0.5, 1, 2, 5, 10))
#> [1] 0.66666667 0.50000000 0.33333333 0.16666667 0.09090909
```

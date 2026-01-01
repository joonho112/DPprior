# Conditional Variance of rho Given Alpha

Computes the conditional variance Var(rho \| alpha).

## Usage

``` r
var_rho_given_alpha(alpha)
```

## Arguments

- alpha:

  Numeric vector; concentration parameter(s) (must be positive).

## Value

Numeric vector; Var(rho \| alpha).

## Details

The conditional variance is: \$\$Var(\rho \| \alpha) =
\frac{2\alpha}{(1+\alpha)^2(2+\alpha)(3+\alpha)}\$\$

This is derived from the GEM recursion: rho = V^2 + (1-V)^2 \* rho'
where V ~ Beta(1, alpha) and rho' is an independent copy of rho.

**Properties:**

- Var(rho\|alpha) = 0 when alpha -\> 0 (degenerate at rho = 1)

- Var(rho\|alpha) -\> 0 when alpha -\> Inf (degenerate at rho = 0)

- Maximum variance occurs at intermediate alpha

## References

Lee, J. (2025). RN-06: Dual-Anchor Design II, Appendix A.2.

## See also

[`mean_rho_given_alpha`](https://joonho112.github.io/DPprior/reference/mean_rho_given_alpha.md),
[`var_rho`](https://joonho112.github.io/DPprior/reference/var_rho.md)

## Examples

``` r
var_rho_given_alpha(2)
#> [1] 0.02222222
var_rho_given_alpha(c(0.5, 1, 2, 5, 10))
#> [1] 0.050793651 0.041666667 0.022222222 0.004960317 0.001059547
```

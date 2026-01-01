# Derivative of Conditional Variance w.r.t. Alpha

Computes \\\frac{d}{d\alpha} \mathrm{Var}(K_J \mid \alpha)\\.

## Usage

``` r
dvar_dalpha(J, alpha)
```

## Arguments

- J:

  Sample size (integer \>= 1).

- alpha:

  Concentration parameter (positive numeric, vectorized).

## Value

Numeric vector of derivatives (same length as `alpha`).

## Details

Uses the summation form: \$\$\frac{d}{d\alpha} v_J(\alpha) =
\sum\_{r=1}^{J-1} \frac{r(r - \alpha)}{(\alpha + r)^3}\$\$

This derivative can be positive or negative depending on \\\alpha\\.

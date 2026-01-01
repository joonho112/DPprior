# Conditional Variance of K Given Alpha (Enhanced)

Computes \\Var(K_J \| \alpha)\\ with proper handling of very small alpha
values.

## Usage

``` r
var_K_given_alpha_safe(J, alpha)
```

## Arguments

- J:

  Integer; sample size.

- alpha:

  Numeric vector; concentration parameter values.

## Value

Numeric vector of conditional variances.

## Details

For very small alpha (\< 1e-12), returns 0.0 as the limiting value.

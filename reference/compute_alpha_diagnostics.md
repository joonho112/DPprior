# Alpha Distribution Diagnostics

Computes summary statistics for alpha distributed as Gamma(a, b).

## Usage

``` r
compute_alpha_diagnostics(a, b)
```

## Arguments

- a:

  Numeric; shape parameter of the Gamma prior (\> 0).

- b:

  Numeric; rate parameter of the Gamma prior (\> 0).

## Value

A list with components: mean (expected value
E[scales::alpha](https://scales.r-lib.org/reference/alpha.html) = a/b),
sd (standard deviation), cv (coefficient of variation = 1/sqrt(a)),
median, and quantiles (named vector with q5, q25, q50, q75, q95).

## Details

The coefficient of variation depends only on the shape parameter a,
making it a useful summary of prior uncertainty regardless of the mean.

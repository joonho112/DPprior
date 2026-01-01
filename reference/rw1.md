# Random Generation from w₁ Distribution

Generates random samples from the w₁ distribution by first sampling α ~
Gamma(a, b), then w₁ \| α ~ Beta(1, α).

## Usage

``` r
rw1(n, a, b)
```

## Arguments

- n:

  Integer; number of samples to generate.

- a:

  Numeric; shape parameter of the Gamma prior on α (a \> 0).

- b:

  Numeric; rate parameter of the Gamma prior on α (b \> 0).

## Value

Numeric vector of length n; random samples from the w₁ distribution.

## Details

This function uses the hierarchical representation:

1.  α ~ Gamma(a, b)

2.  w₁ \| α ~ Beta(1, α)

Useful for Monte Carlo validation of the closed-form functions.

## Examples

``` r
# Generate samples
set.seed(42)
samples <- rw1(10000, a = 2, b = 1)

# Compare empirical vs theoretical mean
mean(samples)          # ~0.404
#> [1] 0.4023996
mean_w1(a = 2, b = 1)  # 0.4037
#> [1] 0.4036526
```

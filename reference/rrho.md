# Random Generation from rho Distribution

Generates random samples from the rho = sum w_h^2 distribution by
stick-breaking simulation.

## Usage

``` r
rrho(n, a, b, n_sticks = 500L)
```

## Arguments

- n:

  Integer; number of samples to generate.

- a:

  Numeric; shape parameter of the Gamma prior on alpha (a \> 0).

- b:

  Numeric; rate parameter of the Gamma prior on alpha (b \> 0).

- n_sticks:

  Integer; number of sticks for truncation. Default is 500.

## Value

Numeric vector of length n; random samples from the rho distribution.

## Details

Uses the hierarchical representation:

1.  alpha ~ Gamma(a, b)

2.  v_h \| alpha ~ Beta(1, alpha) independently for h = 1, ..., n_sticks

3.  w_1 = v_1, w_h = v_h \* prod(1 - v_l) for l \< h

4.  rho = sum w_h^2

Useful for Monte Carlo validation of analytical formulas.

## See also

[`rw1`](https://joonho112.github.io/DPprior/reference/rw1.md) for w1
random generation

## Examples

``` r
set.seed(42)
rho_samples <- rrho(1000, a = 2, b = 1)
mean(rho_samples)
#> [1] 0.4096141
mean_rho(a = 2, b = 1)
#> [1] 0.4036526
```

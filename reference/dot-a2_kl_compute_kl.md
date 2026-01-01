# Compute KL Divergence from Target to Induced PMF

Compute KL Divergence from Target to Induced PMF

## Usage

``` r
.a2_kl_compute_kl(target_pmf, induced_pmf, eps = 1e-15)
```

## Arguments

- target_pmf:

  Numeric vector; normalized PMF (length J).

- induced_pmf:

  Numeric vector; normalized induced PMF (length J).

- eps:

  Numeric; small constant to avoid log(0).

## Value

Numeric; KL divergence (non-negative).

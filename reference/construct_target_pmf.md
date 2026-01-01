# Construct Target PMF from User Specification

Creates a normalized target PMF from either a user-provided PMF vector
or target moment specification.

## Usage

``` r
construct_target_pmf(J, target)
```

## Arguments

- J:

  Integer; sample size.

- target:

  Either:

  - Numeric vector of length J or J+1: direct PMF specification

  - Named list with `mu_K` and `var_K`: construct from moments

## Value

A list with components:

- `pmf`:

  Numeric vector of length J; normalized target PMF

- `mu_K`:

  Target mean

- `var_K`:

  Target variance

- `df`:

  (if moments provided) Chi-square degrees of freedom

- `scale`:

  (if moments provided) Chi-square scale parameter

## Details

**Direct PMF specification:** If `target` is a numeric vector of length
J, it is treated as the PMF for k = 1, ..., J. If length J+1, the k=0
entry is dropped.

**Moment specification:** If `target` is a list with `mu_K` and `var_K`,
a discretized chi-square distribution matching these moments is
constructed using: \$\$\text{scale} = \sigma^2_K / (2\mu_K), \quad
\text{df} = 2\mu_K^2 / \sigma^2_K\$\$

## See also

[`discretize_chisq`](https://joonho112.github.io/DPprior/reference/discretize_chisq.md),
[`DPprior_a2_kl`](https://joonho112.github.io/DPprior/reference/DPprior_a2_kl.md)

## Examples

``` r
# Direct PMF
result <- construct_target_pmf(50, rep(1, 50))  # Uniform on 1:50

# Moment specification
result <- construct_target_pmf(50, list(mu_K = 5, var_K = 8))
result$mu_K  # 5
#> [1] 5
result$df    # 6.25
#> [1] 6.25
```

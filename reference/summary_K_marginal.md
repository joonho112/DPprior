# Summary Statistics for Marginal K Distribution

Computes comprehensive summary statistics for the marginal distribution
of \\K_J\\ under a Gamma prior on \\\alpha\\.

## Usage

``` r
summary_K_marginal(
  J,
  a,
  b,
  logS,
  M = .QUAD_NODES_DEFAULT,
  probs = c(0.05, 0.25, 0.5, 0.75, 0.95)
)
```

## Arguments

- J:

  Integer; sample size (positive integer \>= 1).

- a:

  Numeric; shape parameter of Gamma prior (\> 0).

- b:

  Numeric; rate parameter of Gamma prior (\> 0).

- logS:

  Matrix; pre-computed log-Stirling matrix.

- M:

  Integer; number of quadrature nodes (default: 80).

- probs:

  Numeric vector; probability levels for quantiles (default:
  `c(0.05, 0.25, 0.5, 0.75, 0.95)`).

## Value

A list with components:

- `J`:

  Sample size

- `a`:

  Gamma shape parameter

- `b`:

  Gamma rate parameter

- `mean`:

  Mean \\E\[K_J \mid a, b\]\\

- `var`:

  Variance \\Var(K_J \mid a, b)\\

- `sd`:

  Standard deviation

- `cv`:

  Coefficient of variation (sd/mean)

- `mode`:

  Mode (most likely value)

- `median`:

  Median (50th percentile)

- `quantiles`:

  Named integer vector of quantiles at `probs`

- `pmf`:

  Full PMF vector

- `cdf`:

  Full CDF vector

## Details

This function provides a complete summary of the marginal distribution,
combining PMF-based and CDF-based statistics. The mean and variance
computed from the PMF should match
[`exact_K_moments()`](https://joonho112.github.io/DPprior/reference/exact_K_moments.md)
within numerical tolerance.

The `probs` argument allows customization of which quantiles to report,
making it flexible for different reporting needs.

## See also

[`pmf_K_marginal`](https://joonho112.github.io/DPprior/reference/pmf_K_marginal.md),
[`exact_K_moments`](https://joonho112.github.io/DPprior/reference/exact_K_moments.md)

Other marginal_K:
[`K_moments()`](https://joonho112.github.io/DPprior/reference/K_moments.md),
[`cdf_K_marginal()`](https://joonho112.github.io/DPprior/reference/cdf_K_marginal.md),
[`exact_K_moments()`](https://joonho112.github.io/DPprior/reference/exact_K_moments.md),
[`mode_K_marginal()`](https://joonho112.github.io/DPprior/reference/mode_K_marginal.md),
[`pmf_K_marginal()`](https://joonho112.github.io/DPprior/reference/pmf_K_marginal.md),
[`quantile_K_marginal()`](https://joonho112.github.io/DPprior/reference/quantile_K_marginal.md)

## Examples

``` r
logS <- compute_log_stirling(50)
summary <- summary_K_marginal(50, 1.5, 0.5, logS)

# View main statistics
summary$mean
#> [1] 8.355487
summary$var
#> [1] 22.76895
summary$mode
#> [1] 6
summary$quantiles
#>  q5 q25 q50 q75 q95 
#>   2   5   8  11  17 

# Custom quantiles
summary2 <- summary_K_marginal(50, 1.5, 0.5, logS,
                               probs = c(0.025, 0.5, 0.975))
summary2$quantiles
#>  q2 q50 q98 
#>   1   8  19 

# Compare with exact moments
exact <- exact_K_moments(50, 1.5, 0.5)
c(summary$mean - exact$mean, summary$var - exact$var)
#> [1] -7.105427e-15 -4.263256e-14
```

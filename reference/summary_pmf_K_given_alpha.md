# Summary of Conditional PMF

Returns a comprehensive summary of the conditional distribution \\K_J
\mid \alpha\\.

## Usage

``` r
summary_pmf_K_given_alpha(J, alpha, logS)
```

## Arguments

- J:

  Integer; sample size.

- alpha:

  Numeric; DP concentration parameter.

- logS:

  Matrix; pre-computed log-Stirling matrix.

## Value

A list with components:

- `J`:

  Sample size

- `alpha`:

  Concentration parameter

- `mean`:

  Conditional mean \\E\[K_J \mid \alpha\]\\

- `var`:

  Conditional variance \\Var(K_J \mid \alpha)\\

- `sd`:

  Conditional standard deviation

- `mode`:

  Most likely value of K

- `median`:

  Median of K

- `quantiles`:

  25th, 50th, and 75th percentiles

- `pmf`:

  Full PMF vector

- `cdf`:

  Full CDF vector

## Examples

``` r
if (FALSE) { # \dontrun{
logS <- compute_log_stirling(50)
summary <- summary_pmf_K_given_alpha(50, 2.0, logS)
print(summary)

} # }
```

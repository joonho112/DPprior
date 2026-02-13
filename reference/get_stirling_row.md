# Get Stirling Numbers for a Fixed J

Returns a vector of \\\log\|s(J,k)\|\\ for \\k = 1, \ldots, J\\.

## Usage

``` r
get_stirling_row(J, logS)
```

## Arguments

- J:

  Sample size (positive integer).

- logS:

  Pre-computed log-Stirling matrix from
  [`compute_log_stirling`](https://joonho112.github.io/DPprior/reference/compute_log_stirling.md).

## Value

Numeric vector of length J containing \\\log\|s(J,k)\|\\ for \\k = 1,
\ldots, J\\.

## Examples

``` r
if (FALSE) { # \dontrun{
logS <- compute_log_stirling(10)
log_s_row <- get_stirling_row(5, logS)
exp(log_s_row)

} # }
```

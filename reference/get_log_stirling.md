# Get Single Log-Stirling Value with Bounds Checking

Safe accessor for the log-Stirling matrix with automatic bounds
checking. Returns `-Inf` for invalid indices (k \> J, k \< 1).

## Usage

``` r
get_log_stirling(J, k, logS)
```

## Arguments

- J:

  Sample size (non-negative integer).

- k:

  Number of clusters (integer).

- logS:

  Pre-computed log-Stirling matrix from
  [`compute_log_stirling`](https://joonho112.github.io/DPprior/reference/compute_log_stirling.md).

## Value

The value \\\log\|s(J,k)\|\\, or `-Inf` if k \> J or k \< 1.

## See also

[`compute_log_stirling`](https://joonho112.github.io/DPprior/reference/compute_log_stirling.md)
for matrix computation

## Examples

``` r
if (FALSE) { # \dontrun{
logS <- compute_log_stirling(10)

# Valid access
get_log_stirling(4, 2, logS)

# Invalid access returns -Inf
get_log_stirling(4, 5, logS)
get_log_stirling(4, 0, logS)

} # }
```

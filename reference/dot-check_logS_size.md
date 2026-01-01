# Validate Precomputed Log-Stirling Matrix Size

Checks that the log-Stirling matrix is properly formatted and large
enough for the given sample size J.

## Usage

``` r
.check_logS_size(J, logS)
```

## Arguments

- J:

  Sample size (integer \>= 1).

- logS:

  Pre-computed log-Stirling matrix from
  [`compute_log_stirling()`](https://joonho112.github.io/DPprior/reference/compute_log_stirling.md).

## Value

Invisible `TRUE` if validation passes.

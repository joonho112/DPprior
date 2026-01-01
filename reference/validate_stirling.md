# Validate Stirling Number Computation

Validates a log-Stirling matrix against known reference values. Useful
for testing and verification.

## Usage

``` r
validate_stirling(logS, verbose = TRUE)
```

## Arguments

- logS:

  Pre-computed log-Stirling matrix from
  [`compute_log_stirling`](https://joonho112.github.io/DPprior/reference/compute_log_stirling.md).

- verbose:

  Logical; if `TRUE`, print validation results.

## Value

Logical indicating whether all validations passed.

## Details

Checks against known values:

- \\\|s(4,2)\| = 11\\

- \\\|s(5,3)\| = 35\\

- \\\|s(10,5)\| = 269325\\

## Examples

``` r
logS <- compute_log_stirling(10)
validate_stirling(logS)
#> PASS: |s(4,2)| = 11 (expected 11)
#> PASS: |s(5,3)| = 35 (expected 35)
#> PASS: |s(6,3)| = 225 (expected 225)
#> PASS: |s(10,5)| = 269325 (expected 269325)
#> [1] TRUE
```

# Verify Row Sum Identity for Stirling Numbers

Verifies that the row sum of Stirling numbers equals J! for each row.
This is the fundamental identity: \\\sum\_{k=1}^{J} \|s(J,k)\| = J!\\

## Usage

``` r
verify_stirling_row_sum(
  logS,
  J_values = 2:10,
  tolerance = 1e-10,
  verbose = TRUE
)
```

## Arguments

- logS:

  Pre-computed log-Stirling matrix from
  [`compute_log_stirling`](https://joonho112.github.io/DPprior/reference/compute_log_stirling.md).

- J_values:

  Vector of J values to verify (default: 2:10).

- tolerance:

  Numerical tolerance for comparison (default: 1e-10).

- verbose:

  Logical; if `TRUE`, print verification results.

## Value

Logical indicating whether all verifications passed.

## Details

The unsigned Stirling numbers of the first kind satisfy:
\$\$\sum\_{k=1}^{J} \|s(J,k)\| = J!\$\$

This identity follows from the fact that \\\|s(J,k)\|\\ counts
permutations of J elements with exactly k cycles, and the total number
of permutations is J!.

## Examples

``` r
if (FALSE) { # \dontrun{
logS <- compute_log_stirling(15)
verify_stirling_row_sum(logS, J_values = 2:10)

} # }
```

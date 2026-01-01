# Assert Valid Cluster Count k

Validates that k is a valid cluster count for given sample size J.

## Usage

``` r
assert_valid_k(k, J)
```

## Arguments

- k:

  Cluster count to validate.

- J:

  Sample size (k must be in 1:J).

## Value

Invisible `TRUE` if validation passes.

# Numerically Stable Softmax

Computes the softmax transformation of a numeric vector, returning a
probability vector that sums to 1.

## Usage

``` r
softmax(x)
```

## Arguments

- x:

  Numeric vector of log-odds or arbitrary real values.

## Value

Numeric vector of probabilities summing to 1.

## Details

The softmax function is defined as: \$\$p_i = \frac{\exp(x_i)}{\sum_j
\exp(x_j)}\$\$

This implementation subtracts the maximum value before exponentiating to
ensure numerical stability for extreme inputs.

Special cases:

- If `x` contains `Inf` values, the probability mass is split uniformly
  across all `Inf` entries.

- Empty vector throws an error.

## Examples

``` r
if (FALSE) { # \dontrun{
softmax(c(1, 2, 3))
sum(softmax(c(1, 2, 3)))

# Works with extreme values
softmax(c(1000, 1001, 1002))

# Inf handling
softmax(c(1, Inf, Inf))

} # }
```

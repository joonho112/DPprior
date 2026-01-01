# Vectorized Log-Sum-Exp

Computes `log(sum(exp(x)))` for a numeric vector in a numerically stable
way.

## Usage

``` r
logsumexp_vec(x)
```

## Arguments

- x:

  Numeric vector of log-scale values.

## Value

Scalar value equal to `log(sum(exp(x)))`.

## Details

Subtracts the maximum before exponentiating to prevent overflow:
\$\$\log\sum_i \exp(x_i) = \max_i x_i + \log\sum_i \exp(x_i - \max_i
x_i)\$\$

Special cases:

- If all entries are `-Inf`, returns `-Inf`.

- If any entry is `Inf`, returns `Inf`.

- Empty vector throws an error.

## See also

[`logsumexp`](https://joonho112.github.io/DPprior/reference/logsumexp.md)
for binary operation

## Examples

``` r
# Sum of equal values
logsumexp_vec(c(0, 0, 0, 0))
#> [1] 1.386294

# Extreme values
logsumexp_vec(c(1000, 1000, 1000))
#> [1] 1001.099

# All -Inf
logsumexp_vec(c(-Inf, -Inf))
#> [1] -Inf
```

# Numerically Stable Log-Sum-Exp (Binary)

Computes `log(exp(a) + exp(b))` in a numerically stable way, avoiding
overflow and underflow.

## Usage

``` r
logsumexp(a, b)
```

## Arguments

- a:

  Numeric vector of log-scale values.

- b:

  Numeric vector of log-scale values (recycled to match length of `a`).

## Value

Numeric vector of `log(exp(a) + exp(b))`.

## Details

Uses the identity: \$\$\log(\exp(a) + \exp(b)) = \max(a,b) + \log(1 +
\exp(-\|a-b\|))\$\$

This formulation ensures numerical stability even for extreme values
(e.g., `a = 1000` or `a = -1000`).

Special cases:

- If both inputs are `-Inf`, returns `-Inf`.

- If either input is `Inf`, returns `Inf`.

## See also

[`logsumexp_vec`](https://joonho112.github.io/DPprior/reference/logsumexp_vec.md)
for vector input

## Examples

``` r
# Standard case
logsumexp(log(2), log(3))
#> [1] 1.609438

# Extreme values that would overflow with naive implementation
logsumexp(1000, 1000)
#> [1] 1000.693

# Edge cases with Inf
logsumexp(-Inf, -Inf)
#> [1] -Inf
```

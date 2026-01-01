# Verify Score Function Zero Expectation Property

Verifies the fundamental property that \\E\[s\_\theta(\alpha)\] = 0\\
for both score functions.

## Usage

``` r
verify_score_expectation(a, b, M = .QUAD_NODES_VERIFICATION, verbose = TRUE)
```

## Arguments

- a:

  Numeric; shape parameter.

- b:

  Numeric; rate parameter.

- M:

  Integer; number of quadrature nodes.

- verbose:

  Logical; if TRUE, print results.

## Value

A named list with components:

- `E_score_a`:

  Expectation of \\s_a\\

- `E_score_b`:

  Expectation of \\s_b\\

## Details

This is a fundamental property of score functions. Due to quadrature
approximation error, the computed expectations may not be exactly zero.

**Expected behavior:**

- `E[s_b]` should be very close to zero (typically \< 1e-14) because s_b
  is linear in alpha.

- `E[s_a]` may show larger errors (up to 1e-2 for small a) due to the
  log(alpha) term causing slower quadrature convergence.

For rigorous verification, use adaptive integration (e.g., in Python
with scipy.integrate.quad) which provides independent ground truth.

## Examples

``` r
verify_score_expectation(a = 2.0, b = 1.0, verbose = TRUE)
#> Score Expectation Verification (a=2.00, b=1.00, M=200)
#>   E[s_a(alpha)] = 1.448463e-05 (should be ~= 0)
#>   E[s_b(alpha)] = -1.197109e-15 (should be ~= 0)
```

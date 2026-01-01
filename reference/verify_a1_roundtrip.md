# Verify A1 Mapping via Round-Trip

Tests the A1 mapping by computing the forward model (NegBin moments)
from the derived \\(a, b)\\ parameters and comparing to targets.

## Usage

``` r
verify_a1_roundtrip(fit, tol = 1e-08, verbose = TRUE)
```

## Arguments

- fit:

  A `DPprior_fit` object from `DPprior_a1`.

- tol:

  Numeric; tolerance for relative error comparison.

- verbose:

  Logical; if TRUE, print verification details.

## Value

Logical; TRUE if round-trip succeeds within tolerance.

## Details

Under the A1 NegBin approximation: \$\$K_J - 1 \sim \text{NegBin}(a,
p)\$\$ where \\p = b/(b + c_J)\\.

## Examples

``` r
fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)
verify_a1_roundtrip(fit)
#> A1 Round-Trip Verification
#> ---------------------------------------- 
#> mu_K: target = 5.000000, recovered = 5.000000, rel_err = 0.00e+00
#> var_K: target = 8.000000, recovered = 8.000000, rel_err = 0.00e+00
#> Result: PASS
```

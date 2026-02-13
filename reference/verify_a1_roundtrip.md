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
if (FALSE) { # \dontrun{
fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)
verify_a1_roundtrip(fit)

} # }
```

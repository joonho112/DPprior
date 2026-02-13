# Compare Exact Moments to NegBin Approximation

Compares exact marginal moments (via quadrature) to the Negative
Binomial approximation from the A1 method for error analysis.

## Usage

``` r
compare_to_negbin(J, a, b, M = .QUAD_NODES_DEFAULT)
```

## Arguments

- J:

  Integer; sample size.

- a:

  Numeric; shape parameter of Gamma prior.

- b:

  Numeric; rate parameter of Gamma prior.

- M:

  Integer; number of quadrature nodes (default: 80).

## Value

A list with components:

- `exact`:

  List with exact mean and var

- `negbin`:

  List with NegBin approximation mean and var

- `abs_error`:

  Absolute errors (negbin - exact)

- `rel_error`:

  Relative errors

## Details

The NegBin approximation (A1 method from Lee, 2026, Section 3.1)
assumes: \$\$K_J - 1 \| \alpha \approx \text{Poisson}(\alpha \cdot
c_J)\$\$

where \\c_J = \log(J)\\. With \\\alpha \sim \text{Gamma}(a, b)\\:
\$\$E\[K_J\] \approx 1 + (a/b) \cdot c_J\$\$ \$\$Var(K_J) \approx m
\cdot (1 + m/a), \quad m = (a/b) \cdot c_J\$\$

This comparison helps diagnose when the A1 approximation is insufficient
and exact A2 moment matching is needed.

## Examples

``` r
if (FALSE) { # \dontrun{
# Large approximation error for small J
compare_to_negbin(50, 1.5, 0.5)

# Error decreases with J
compare_to_negbin(300, 1.5, 0.5)

} # }
```

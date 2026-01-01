# Compute Scaling Constant for A1 Mapping

Computes the scaling constant \\c_J\\ used in the A1 closed-form
mapping. Three variants are supported based on asymptotic
approximations.

## Usage

``` r
compute_scaling_constant(
  J,
  scaling = c("log", "harmonic", "digamma"),
  mu_K = NULL
)
```

## Arguments

- J:

  Integer; number of items/sites (must be \>= 2).

- scaling:

  Character; one of "log", "harmonic", or "digamma".

- mu_K:

  Numeric; target mean of K (required for "digamma" scaling).

## Value

Numeric scalar; the scaling constant \\c_J\\.

## Details

The scaling constant appears in the Poisson proxy: \$\$K_J - 1 \mid
\alpha \approx \text{Poisson}(\alpha \cdot c_J)\$\$

Available variants:

- log:

  \\c_J = \log(J)\\, the asymptotic leading term (default)

- harmonic:

  \\c_J = H\_{J-1} = \psi(J) + \gamma\\, improves accuracy for
  small/moderate J

- digamma:

  \\c_J = \psi(\tilde{\alpha} + J) - \psi(\tilde{\alpha})\\ where
  \\\tilde{\alpha} = (\mu_K - 1)/\log(J)\\, a local correction

## See also

[`DPprior_a1`](https://joonho112.github.io/DPprior/reference/DPprior_a1.md)
for the main elicitation function

## Examples

``` r
# Default log scaling
compute_scaling_constant(50, "log")
#> [1] 3.912023

# Harmonic scaling (better for moderate J)
compute_scaling_constant(50, "harmonic")
#> [1] 4.479205

# Digamma scaling (requires mu_K)
compute_scaling_constant(50, "digamma", mu_K = 5)
#> [1] 4.463254
```

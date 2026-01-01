# A1 Closed-Form Prior Elicitation

Maps target beliefs about the number of clusters \\(\mu_K, \sigma^2_K)\\
to Gamma hyperprior parameters \\(a, b)\\ using the A1 closed-form
approximation based on Negative Binomial moment matching.

## Usage

``` r
DPprior_a1(
  J,
  mu_K,
  var_K,
  scaling = c("log", "harmonic", "digamma"),
  epsilon = .TOL_PROJECTION_BUFFER
)
```

## Arguments

- J:

  Integer; number of items/sites (must be \>= 2).

- mu_K:

  Numeric; target prior mean of \\K_J\\ (must be \> 1 and \<= J).

- var_K:

  Numeric; target prior variance of \\K_J\\ (must be \> 0).

- scaling:

  Character; scaling constant method: "log" (default), "harmonic", or
  "digamma".

- epsilon:

  Numeric; buffer for feasibility projection. Default is
  `.TOL_PROJECTION_BUFFER` (1e-6).

## Value

An S3 object of class `DPprior_fit` with components:

- a:

  Shape parameter of the Gamma prior

- b:

  Rate parameter of the Gamma prior

- J:

  Sample size used

- target:

  List with target moments and type

- method:

  "A1" indicating closed-form method

- status:

  "success" or "projected" if boundary adjustment needed

- scaling:

  Scaling method used

- cJ:

  Scaling constant value

- var_K_used:

  Actual variance used (may differ if projected)

- converged:

  Always TRUE for A1 (for A2 compatibility)

- iterations:

  Always 0L for A1 (for A2 compatibility)

- fit:

  NULL (placeholder for A2 refinement)

- diagnostics:

  NULL (placeholder for diagnostics)

- trace:

  NULL (placeholder for optimization trace)

## Details

### Theory (RN-03)

The A1 method uses a shifted Negative Binomial approximation: \$\$K_J -
1 \mid \alpha \approx \text{Poisson}(\alpha \cdot c_J)\$\$

With \\\alpha \sim \text{Gamma}(a, b)\\, the marginal becomes: \$\$K_J -
1 \approx \text{NegBin}(a, b/(b + c_J))\$\$

### Inverse Formulas (Theorem 1)

Let \\m = \mu_K - 1\\ (shifted mean) and \\D = \sigma^2_K - m\\. If \\D
\> 0\\ (overdispersion): \$\$a = m^2 / D, \quad b = m \cdot c_J / D\$\$

If \\D \leq 0\\ (infeasible), the variance is projected to the boundary.

### Feasibility

The NegBin model requires overdispersion: \\\sigma^2_K \> \mu_K - 1\\.
High-confidence specifications (low variance) may violate this
constraint under the A1 proxy, even though they may be feasible under
the exact DP.

## See also

[`vif_to_variance`](https://joonho112.github.io/DPprior/reference/vif_to_variance.md)
for VIF conversion,
[`confidence_to_vif`](https://joonho112.github.io/DPprior/reference/confidence_to_vif.md)
for confidence mapping,
[`print.DPprior_fit`](https://joonho112.github.io/DPprior/reference/print.DPprior_fit.md)
for print method

## Examples

``` r
# Basic usage with moment targets
fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)
print(fit)
#> DPprior Prior Elicitation Result
#> ============================================= 
#> 
#> Gamma Hyperprior: α ~ Gamma(a = 4.0000, b = 3.9120)
#>   E[α] = 1.022, SD[α] = 0.511
#> 
#> Target (J = 50):
#>   E[K_J]   = 5.00
#>   Var(K_J) = 8.00
#> 
#> Method: A1 (0 iterations)

# Using VIF specification
fit <- DPprior_a1(J = 50, mu_K = 5, var_K = vif_to_variance(5, 2))

# Using confidence-based specification
vif <- confidence_to_vif("medium")
fit <- DPprior_a1(J = 50, mu_K = 5, var_K = vif_to_variance(5, vif))

# Infeasible variance (triggers projection)
fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 3)  # var < mu-1 = 4
#> Warning: var_K <= mu_K - 1: projected to feasible boundary

# Compare scaling methods
fit_log <- DPprior_a1(50, 5, 8, scaling = "log")
fit_harm <- DPprior_a1(50, 5, 8, scaling = "harmonic")
```

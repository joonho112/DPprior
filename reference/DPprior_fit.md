# Fit a Gamma Hyperprior for DP Concentration Parameter

Main entry point for DPprior prior elicitation. Takes user-specified
cluster count expectations and returns calibrated Gamma(a, b)
parameters.

## Usage

``` r
DPprior_fit(
  J,
  mu_K,
  var_K = NULL,
  confidence = c("medium", "low", "high"),
  method = c("A2-MN", "A1", "A2-KL"),
  target_pmf = NULL,
  check_diagnostics = TRUE,
  warn_dominance = TRUE,
  M = .QUAD_NODES_DEFAULT,
  verbose = FALSE,
  ...
)
```

## Arguments

- J:

  Integer; sample size (number of sites/units). Must be positive and not
  exceed the maximum supported value (default: 500).

- mu_K:

  Numeric; target expected number of clusters `E[K_J]`. Must be in the
  range (1, J\]. Note: mu_K = 1 is trivial (single cluster).

- var_K:

  Numeric; target variance of `K_J`. If NULL, computed from the
  `confidence` argument. Must satisfy var_K \<= (J-1)^2/4 (maximum
  possible variance for K in {1,...,J}).

- confidence:

  Character; alternative to var_K for specifying uncertainty. One of:

  - "low": VIF = 5.0 (very high uncertainty, wide prior)

  - "medium": VIF = 2.5 (moderate uncertainty)

  - "high": VIF = 1.5 (low uncertainty, concentrated prior)

  Only used if `var_K` is NULL.

- method:

  Character; calibration method. One of:

  - "A2-MN" (default): Exact moment matching via Newton's method.
    Recommended for accurate calibration.

  - "A1": Fast closed-form approximation based on shifted NegBin. Good
    for initial exploration or large J. Note: may project var_K to
    feasible region if var_K \< mu_K - 1.

  - "A2-KL": KL divergence minimization. Supports both moment target
    (default) and custom target PMF via `target_pmf`.

- target_pmf:

  Numeric vector; optional custom target PMF for A2-KL. If provided,
  A2-KL uses PMF matching mode. Length must equal J.

- check_diagnostics:

  Logical; if TRUE (default), compute comprehensive diagnostics
  including weight distribution analysis.

- warn_dominance:

  Logical; if TRUE (default), issue a warning if the prior implies high
  probability of a dominant cluster (P(w1 \> 0.5) \> 40%).

- M:

  Integer; number of quadrature nodes for numerical integration. Default
  is 80, which provides good accuracy for most cases.

- verbose:

  Logical; if TRUE, print progress messages during calibration.

- ...:

  Additional arguments passed to method-specific functions. Note: only
  matching formal arguments are forwarded to backends.

## Value

An S3 object of class "DPprior_fit" containing:

- a:

  Numeric; shape parameter of the calibrated Gamma hyperprior.

- b:

  Numeric; rate parameter of the calibrated Gamma hyperprior.

- J:

  Integer; sample size used for calibration.

- target:

  List; target specification including mu_K, var_K, confidence level (if
  used), and target type.

- method:

  Character; calibration method used.

- converged:

  Logical; whether the calibration converged.

- iterations:

  Integer; number of iterations (for iterative methods).

- fit:

  List; achieved moments (mu_K, var_K) and residual.

- diagnostics:

  List; comprehensive RN-07 diagnostics (if requested).

- solver_diagnostics:

  List; backend solver details (if applicable).

- trace:

  Data frame; iteration history (for iterative methods).

## Details

This function provides a convenient interface for specifying beliefs
about the number of clusters and uncertainty, then calibrates a Gamma
hyperprior for the Dirichlet Process concentration parameter alpha.

### Variance Constraints

The target variance must satisfy two constraints:

1.  Upper bound: `var_K <= (J-1)^2/4` (maximum possible variance for a
    distribution on {1,...,J})

2.  Lower bound (A1 only): `var_K >= mu_K - 1` (NegBin feasibility). For
    A1, infeasible variance is projected to the boundary with a warning.
    A2-MN does not have this constraint.

### Diagnostics

When `check_diagnostics = TRUE`, the function computes:

- Alpha distribution: `E[alpha]`, `Var(alpha)`, `CV(alpha)`

- K distribution: achieved moments, mode, quantiles

- Weight distribution: `E[w_1]`, `P(w_1 > 0.5)`, dominance risk

- Co-clustering: `E[rho]`, interpretation

Backend solver details are preserved in `solver_diagnostics`.

## References

Lee, J. (2025). Prior Elicitation for the Dirichlet Process
Concentration Parameter in Low-Information Settings. *Working Paper*.

RN-03: Closed-Form Mapping for the Gamma Hyperprior. RN-04: Small-J
Correction via Newton Refinement. RN-07: Unintended Prior Diagnostic.

## See also

[`DPprior_a1`](https://joonho112.github.io/DPprior/reference/DPprior_a1.md)
for A1 closed-form approximation,
[`DPprior_a2_newton`](https://joonho112.github.io/DPprior/reference/DPprior_a2_newton.md)
for A2 Newton refinement,
[`DPprior_a2_kl`](https://joonho112.github.io/DPprior/reference/DPprior_a2_kl.md)
for A2-KL divergence minimization,
[`DPprior_diagnostics`](https://joonho112.github.io/DPprior/reference/DPprior_diagnostics.md)
for detailed diagnostics,
[`DPprior_dual`](https://joonho112.github.io/DPprior/reference/DPprior_dual.md)
for dual-anchor calibration with weight constraints.

## Examples

``` r
# Basic usage with target moments
fit <- DPprior_fit(J = 50, mu_K = 5, var_K = 8)
#> Warning: HIGH DOMINANCE RISK: P(w1 > 0.5) = 48.1% exceeds 40%.
#>   This may indicate unintended prior behavior (RN-07).
#>   Consider using DPprior_dual() for weight-constrained elicitation.
#>   See ?DPprior_diagnostics for interpretation.
print(fit)
#> DPprior Prior Elicitation Result
#> ============================================= 
#> 
#> Gamma Hyperprior: α ~ Gamma(a = 2.0361, b = 1.6051)
#>   E[α] = 1.269, SD[α] = 0.889
#> 
#> Target (J = 50):
#>   E[K_J]   = 5.00
#>   Var(K_J) = 8.00
#> 
#> Achieved:
#>   E[K_J] = 5.000000, Var(K_J) = 8.000000
#>   Residual = 7.60e-09
#> 
#> Method: A2-MN (6 iterations)
#> 
#> Dominance Risk: HIGH ✘ (P(w₁>0.5) = 48%)

# Using confidence level instead of explicit variance
fit_medium <- DPprior_fit(J = 50, mu_K = 5, confidence = "medium")
#> Warning: HIGH DOMINANCE RISK: P(w1 > 0.5) = 49.7% exceeds 40%.
#>   This may indicate unintended prior behavior (RN-07).
#>   Consider using DPprior_dual() for weight-constrained elicitation.
#>   See ?DPprior_diagnostics for interpretation.
fit_low <- DPprior_fit(J = 50, mu_K = 5, confidence = "low")
#> Warning: HIGH DOMINANCE RISK: P(w1 > 0.5) = 56.3% exceeds 40%.
#>   This may indicate unintended prior behavior (RN-07).
#>   Consider using DPprior_dual() for weight-constrained elicitation.
#>   See ?DPprior_diagnostics for interpretation.

# Quick approximation for exploration
fit_quick <- DPprior_fit(J = 50, mu_K = 5, var_K = 8, method = "A1")
#> Warning: HIGH DOMINANCE RISK: P(w1 > 0.5) = 52.1% exceeds 40%.
#>   This may indicate unintended prior behavior (RN-07).
#>   Consider using DPprior_dual() for weight-constrained elicitation.
#>   See ?DPprior_diagnostics for interpretation.

# Verbose output for debugging
fit_verbose <- DPprior_fit(J = 50, mu_K = 5, var_K = 8, verbose = TRUE)
#> Using A2-MN Newton refinement
#> A2-MN Newton Solver
#> Target: E[K]=5.0000, Var(K)=8.0000
#> A1 initialization: a0=4.000000, b0=3.912023
#> -------------------------------------------------------------------------------- 
#> Iter |          a |          b |       E[K] |     Var(K) |      ||F|| |     step |     det(J)
#> -------------------------------------------------------------------------------- 
#>    1 |   4.000000 |   3.912023 |   4.461351 |   4.783136 |   3.26e+00 |   1.0000 |  -5.30e+00
#>    2 |   1.178650 |   0.911969 |   4.909046 |  10.854537 |   2.86e+00 |   1.0000 |  -2.16e+01
#>    3 |   1.844384 |   1.455254 |   4.974913 |   8.399473 |   4.00e-01 |   1.0000 |  -1.53e+01
#>    4 |   2.029223 |   1.599680 |   4.999187 |   8.013243 |   1.33e-02 |   1.0000 |  -1.43e+01
#>    5 |   2.036082 |   1.605046 |   4.999999 |   8.000021 |   2.08e-05 |   1.0000 |  -1.43e+01
#>    6 |   2.036093 |   1.605054 |   5.000000 |   8.000000 |   7.60e-09 |      --- |  -1.43e+01
#> 
#> Converged: ||F|| = 7.60e-09 < 1.00e-08
#> Warning: HIGH DOMINANCE RISK: P(w1 > 0.5) = 48.1% exceeds 40%.
#>   This may indicate unintended prior behavior (RN-07).
#>   Consider using DPprior_dual() for weight-constrained elicitation.
#>   See ?DPprior_diagnostics for interpretation.

# Access results
fit$a  # Gamma shape
#> [1] 2.036093
fit$b  # Gamma rate
#> [1] 1.605054
fit$fit$mu_K  # Achieved mean
#> [1] 5
fit$diagnostics$weights$dominance_risk  # Dominance risk level
#> [1] "high"
```

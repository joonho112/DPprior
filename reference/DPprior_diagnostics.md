# Comprehensive Prior Diagnostics

Computes a full diagnostic report for a fitted DPprior object,
implementing the "unintended prior" checks from RN-07.

## Usage

``` r
DPprior_diagnostics(fit, thresholds = c(0.5, 0.9))
```

## Arguments

- fit:

  A DPprior_fit object from any calibration method (e.g., DPprior_a1,
  DPprior_a2_newton, DPprior_dual). Must contain fields a, b, and J.
  Optionally can contain M for quadrature nodes.

- thresholds:

  Numeric vector; thresholds for weight dominance checks (default:
  c(0.5, 0.9)).

## Value

An S3 object of class "DPprior_diagnostics" with components: J (sample
size), a and b (Gamma parameters), alpha (alpha distribution summary), K
(K distribution summary), weights (w1 distribution summary with
dominance risk), coclustering (rho summary), and warnings (character
vector of diagnostic warnings).

## Details

Warnings are issued if:

- P(w1 \> 0.5) \> 0.4: "HIGH DOMINANCE RISK"

- P(w1 \> 0.9) \> 0.15: "NEAR-DEGENERATE RISK"

## References

Lee, J. (2025). RN-07: Unintended Prior Diagnostic.

## Examples

``` r
fit <- DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8)
diag <- DPprior_diagnostics(fit)
print(diag)
#> DPprior Comprehensive Diagnostics
#> ============================================================ 
#> 
#> Prior: alpha ~ Gamma(2.0361, 1.6051) for J = 50
#> 
#> alpha Distribution:
#> ---------------------------------------- 
#>   E[alpha] = 1.269, CV(alpha) = 0.701, Median = 1.068
#>   90% CI: [0.230, 2.992]
#> 
#> K_J Distribution:
#> ---------------------------------------- 
#>   E[K] = 5.00, SD(K) = 2.83, Mode = 3
#>   Median = 5, IQR = [3, 7]
#> 
#> w1 Distribution (Size-Biased First Weight):
#> ---------------------------------------- 
#>   E[w1] = 0.501, Median = 0.478
#>   P(w1 > 0.5) = 48.1% (dominance risk: HIGH)
#>   P(w1 > 0.9) = 16.3%
#> 
#> Co-Clustering (rho = sum w_h^2):
#> ---------------------------------------- 
#>   E[rho] = 0.501 (High prior co-clustering: most unit pairs expected in same cluster)
#> 
#> WARNINGS:
#> ---------------------------------------- 
#>   * HIGH DOMINANCE RISK: P(w1 > 0.5) = 48.1% exceeds 40%
#>   * NEAR-DEGENERATE RISK: P(w1 > 0.9) = 16.3% exceeds 15%
#> 
#>   Consider using DPprior_dual() for weight-constrained elicitation.
```

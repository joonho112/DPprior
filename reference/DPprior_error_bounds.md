# Compute A1 Approximation Error Bounds

Comprehensive error analysis for the A1 large-J approximation.
Implements the error quantification framework from Lee (2026, Section
3.3).

## Usage

``` r
DPprior_error_bounds(J, a, b, cJ = log(J), M = .QUAD_NODES_DEFAULT)
```

## Arguments

- J:

  Integer; sample size.

- a, b:

  Numeric; Gamma hyperparameters (shape, rate).

- cJ:

  Numeric; scaling constant (default: log(J)).

- M:

  Integer; number of quadrature nodes (default: 80).

## Value

An S3 object of class "DPprior_error_bounds" with components:

- J, a, b, cJ:

  Input parameters

- moment_errors:

  List of moment error metrics from `a1_moment_error`

- tv_bounds:

  List with conditional and marginal TV bounds

- recommendation:

  "A1_sufficient" or "A2_recommended"

- threshold_J:

  Estimated J threshold for A1 adequacy (or NA)

## Details

The recommendation is based on:

- A1 sufficient if: mean relative error \< 5\\

- A2 recommended otherwise

This function provides:

1.  Moment errors: Exact vs A1 approximation for mean and variance

2.  TV bounds: Conditional bounds at multiple alpha values, plus
    marginal bound

3.  Recommendation: Whether to use A1 or refine with A2

4.  Threshold: Minimum J for A1 adequacy with current prior

## References

Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet
Process Mixtures. *arXiv preprint* arXiv:2602.06301.

## See also

[`a1_moment_error`](https://joonho112.github.io/DPprior/reference/a1_moment_error.md),
[`compute_total_tv_bound`](https://joonho112.github.io/DPprior/reference/compute_total_tv_bound.md),
[`DPprior_a1`](https://joonho112.github.io/DPprior/reference/DPprior_a1.md),
[`DPprior_a2_newton`](https://joonho112.github.io/DPprior/reference/DPprior_a2_newton.md)

Other diagnostics:
[`DPprior_diagnostics()`](https://joonho112.github.io/DPprior/reference/DPprior_diagnostics.md),
[`compute_error_landscape()`](https://joonho112.github.io/DPprior/reference/compute_error_landscape.md),
[`compute_weight_diagnostics()`](https://joonho112.github.io/DPprior/reference/compute_weight_diagnostics.md),
[`dual_anchor_diagnostics()`](https://joonho112.github.io/DPprior/reference/dual_anchor_diagnostics.md)

## Examples

``` r
# Check A1 adequacy for J=50, typical prior
bounds <- DPprior_error_bounds(J = 50, a = 1.6, b = 1.2)
print(bounds)
#> DPprior A1 Approximation Error Analysis
#> ================================================== 
#> 
#> Sample size J = 50, c_J = 3.9120
#> Gamma prior: alpha ~ Gamma(1.6000, 1.2000) [shape-rate]
#> E[alpha] = 1.3333, CV(alpha) = 0.7906
#> 
#> Moment Errors (A1 vs Exact):
#> --------------------------------------------- 
#>   E[K_J]:   exact =   5.0973, A1 =   6.2160, error =  21.95%
#>   Var(K_J): exact =   9.5500, A1 =  22.2204, error = 132.67%
#> 
#> TV Bounds:
#> --------------------------------------------- 
#>   Marginal E[d_TV] <= 0.3642
#>   At E[alpha] = 1.33:
#>     Poissonization: 0.2043 (raw: 0.9128)
#>     Linearization:  0.1804
#>     Total:          0.3846
#> 
#> Recommendation: A2_recommended
#>   (A1 may not achieve < 5%% mean error for J <= 500 with this prior)

# For larger J, A1 becomes more adequate
bounds_200 <- DPprior_error_bounds(J = 200, a = 1.6, b = 1.2)
print(bounds_200)
#> DPprior A1 Approximation Error Analysis
#> ================================================== 
#> 
#> Sample size J = 200, c_J = 5.2983
#> Gamma prior: alpha ~ Gamma(1.6000, 1.2000) [shape-rate]
#> E[alpha] = 1.3333, CV(alpha) = 0.7906
#> 
#> Moment Errors (A1 vs Exact):
#> --------------------------------------------- 
#>   E[K_J]:   exact =   6.9134, A1 =   8.0644, error =  16.65%
#>   Var(K_J): exact =  20.3210, A1 =  38.2557, error =  88.26%
#> 
#> TV Bounds:
#> --------------------------------------------- 
#>   Marginal E[d_TV] <= 0.2995
#>   At E[alpha] = 1.33:
#>     Poissonization: 0.1500 (raw: 0.9389)
#>     Linearization:  0.1571
#>     Total:          0.3071
#> 
#> Recommendation: A2_recommended
#>   (A1 may not achieve < 5%% mean error for J <= 500 with this prior)
```

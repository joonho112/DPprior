# Print Method for DPprior_fit Objects

Displays a concise, informative summary of a prior elicitation result,
including the Gamma hyperprior specification, target vs achieved fit,
and dominance risk assessment.

## Usage

``` r
# S3 method for class 'DPprior_fit'
print(x, digits = 4, ...)
```

## Arguments

- x:

  A `DPprior_fit` object.

- digits:

  Integer; number of significant digits for display. Default is 4.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns `x`.

## Details

The output includes:

- Gamma hyperprior parameters (a, b) with moments
  E[scales::alpha](https://scales.r-lib.org/reference/alpha.html) and
  SD[scales::alpha](https://scales.r-lib.org/reference/alpha.html)

- Target specification (J, \\E\[K_J\]\\, \\Var(K_J)\\)

- Achieved fit with residual error

- Method used and iteration count

- Quick dominance risk summary (if diagnostics available)

## Dominance Risk

If diagnostics are computed, the dominance risk is displayed as:

- **LOW**: P(w1 \> 0.5) \< 20\\

- **MODERATE**: 20\\

- **HIGH**: P(w1 \> 0.5) \>= 40\\

## See also

[`summary.DPprior_fit`](https://joonho112.github.io/DPprior/reference/summary.DPprior_fit.md),
[`plot.DPprior_fit`](https://joonho112.github.io/DPprior/reference/plot.DPprior_fit.md),
[`DPprior_fit`](https://joonho112.github.io/DPprior/reference/DPprior_fit.md)

## Examples

``` r
# Create a fit object
fit <- DPprior_fit(J = 50, mu_K = 5, var_K = 8)
#> Warning: HIGH DOMINANCE RISK: P(w1 > 0.5) = 48.1% exceeds 40%.
#>   This may indicate unintended prior behavior (Lee, 2026).
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

# With custom digits
print(fit, digits = 6)
#> DPprior Prior Elicitation Result
#> ============================================= 
#> 
#> Gamma Hyperprior: α ~ Gamma(a = 2.036093, b = 1.605054)
#>   E[α] = 1.269, SD[α] = 0.889
#> 
#> Target (J = 50):
#>   E[K_J]   = 5.00
#>   Var(K_J) = 8.00
#> 
#> Achieved:
#>   E[K_J] = 5.00000000, Var(K_J) = 8.00000001
#>   Residual = 7.60e-09
#> 
#> Method: A2-MN (6 iterations)
#> 
#> Dominance Risk: HIGH ✘ (P(w₁>0.5) = 48%)
```

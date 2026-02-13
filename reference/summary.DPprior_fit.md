# Summary Method for DPprior_fit Objects

Produces a comprehensive summary of a prior elicitation result,
including detailed parameter information, target vs achieved comparison,
and full diagnostic statistics.

## Usage

``` r
# S3 method for class 'DPprior_fit'
summary(object, print_output = TRUE, ...)
```

## Arguments

- object:

  A `DPprior_fit` object.

- print_output:

  Logical; if `TRUE` (default), prints the summary to the console. If
  `FALSE`, returns the summary list silently.

- ...:

  Additional arguments (currently unused).

## Value

An object of class `"summary.DPprior_fit"` containing:

- `method`: Calibration method used

- `status`: Convergence status

- `gamma_prior`: List with (a, b) parameters

- `alpha_summary`: List with
  E[scales::alpha](https://scales.r-lib.org/reference/alpha.html),
  Var[scales::alpha](https://scales.r-lib.org/reference/alpha.html),
  CV[scales::alpha](https://scales.r-lib.org/reference/alpha.html)

- `target`: Target specification

- `achieved`: Achieved fit

- `errors`: Absolute and relative fitting errors

- `converged`: Logical convergence flag

- `iterations`: Number of iterations

- `diagnostics`: Full diagnostic information (if available)

## Details

The printed summary includes:

1.  Basic information (J, method)

2.  Gamma hyperprior parameters with derived statistics

3.  Target vs Achieved comparison table with error metrics

4.  Full diagnostics for alpha, K, and w1 distributions (if computed)

5.  Iteration trace (first 10 rows, if available)

## See also

[`print.DPprior_fit`](https://joonho112.github.io/DPprior/reference/print.DPprior_fit.md),
[`DPprior_diagnostics`](https://joonho112.github.io/DPprior/reference/DPprior_diagnostics.md)

## Examples

``` r
# Create a fit object
fit <- DPprior_fit(J = 50, mu_K = 5, var_K = 8, check_diagnostics = TRUE)
#> Warning: HIGH DOMINANCE RISK: P(w1 > 0.5) = 48.1% exceeds 40%.
#>   This may indicate unintended prior behavior (Lee, 2026).
#>   Consider using DPprior_dual() for weight-constrained elicitation.
#>   See ?DPprior_diagnostics for interpretation.
summary(fit)
#> DPprior Prior Elicitation Summary
#> ============================================================ 
#> 
#> Sample size: J = 50
#> Method: A2-MN
#> Status: success
#> 
#> Gamma Hyperprior:
#> ---------------------------------------- 
#>   Shape (a) = 2.036093
#>   Rate (b)  = 1.605054
#>   E[α] = 1.2686, SD[α] = 0.8890, CV[α] = 0.7008
#> 
#> Target vs Achieved:
#> ---------------------------------------- 
#>                         Target     Achieved        Error
#>   E[K_J]                5.0000       5.0000     8.31e-10
#>   Var(K_J)              8.0000       8.0000     7.55e-09
#> 
#> Diagnostics:
#> ---------------------------------------- 
#>   Alpha:
#>     Mean = 1.2686, SD = 0.8890, CV = 0.7008
#>     90% CI: [0.230, 2.992]
#>   K_J:
#>     Mean = 5.00, SD = 2.83, Mode = 3
#>   First weight (w₁):
#>     Mean = 0.5014
#>     P(w₁ > 0.5) = 48.1%
#>     P(w₁ > 0.9) = 16.3%
#>     Dominance risk: HIGH
#> 
#> Iteration Trace:
#> ---------------------------------------- 
#>   iter        a         b       M1         V     residual step   det_Jlog
#> 1    1 4.000000 3.9120230 4.461351  4.783136 3.261649e+00    1  -5.300969
#> 2    2 1.178650 0.9119694 4.909046 10.854537 2.855986e+00    1 -21.553612
#> 3    3 1.844384 1.4552538 4.974913  8.399473 4.002603e-01    1 -15.313512
#> 4    4 2.029223 1.5996801 4.999187  8.013243 1.326844e-02    1 -14.298307
#> 5    5 2.036082 1.6050455 4.999999  8.000021 2.078492e-05    1 -14.263052
#> 6    6 2.036093 1.6050541 5.000000  8.000000 7.600363e-09   NA -14.262997
#> 
#> Diagnostics: dominance risk = HIGH (use diagnostics=TRUE for full report)
#> 

# Store summary without printing
summ <- summary(fit, print_output = FALSE)
str(summ)
#> List of 13
#>  $ method       : chr "A2-MN"
#>  $ status       : chr "success"
#>  $ gamma_prior  :List of 2
#>   ..$ a: num 2.04
#>   ..$ b: num 1.61
#>  $ alpha_summary:List of 4
#>   ..$ E_alpha  : num 1.27
#>   ..$ Var_alpha: num 0.79
#>   ..$ SD_alpha : num 0.889
#>   ..$ CV_alpha : num 0.701
#>  $ target       :List of 4
#>   ..$ mu_K      : num 5
#>   ..$ var_K     : num 8
#>   ..$ var_K_used: num 8
#>   ..$ confidence: chr NA
#>  $ achieved     :List of 3
#>   ..$ mu_K    : num 5
#>   ..$ var_K   : num 8
#>   ..$ residual: num 7.6e-09
#>  $ errors       :List of 4
#>   ..$ mu_K_abs     : num 8.31e-10
#>   ..$ var_K_abs    : num 7.55e-09
#>   ..$ mu_K_rel_pct : num 1.66e-08
#>   ..$ var_K_rel_pct: num 9.44e-08
#>  $ scaling      :List of 1
#>   ..$ J: num 50
#>  $ converged    : logi TRUE
#>  $ iterations   : int 6
#>  $ diagnostics  :List of 8
#>   ..$ J           : int 50
#>   ..$ a           : num 2.04
#>   ..$ b           : num 1.61
#>   ..$ alpha       :List of 5
#>   .. ..$ mean     : num 1.27
#>   .. ..$ sd       : num 0.889
#>   .. ..$ cv       : num 0.701
#>   .. ..$ median   : num 1.07
#>   .. ..$ quantiles: Named num [1:5] 0.23 0.615 1.068 1.706 2.992
#>   .. .. ..- attr(*, "names")= chr [1:5] "q5" "q25" "q50" "q75" ...
#>   ..$ K           :List of 7
#>   .. ..$ mean     : num 5
#>   .. ..$ var      : num 8
#>   .. ..$ sd       : num 2.83
#>   .. ..$ mode     : int 3
#>   .. ..$ median   : int 5
#>   .. ..$ quantiles: Named int [1:5] 1 3 5 7 10
#>   .. .. ..- attr(*, "names")= chr [1:5] "q5" "q25" "q50" "q75" ...
#>   .. ..$ pmf      : num [1:50] 0.0746 0.1247 0.1473 0.1472 0.132 ...
#>   ..$ weights     :List of 5
#>   .. ..$ mean          : num 0.501
#>   .. ..$ median        : num 0.478
#>   .. ..$ quantiles     : Named num [1:5] 0.0401 0.2162 0.4784 0.7911 0.9954
#>   .. .. ..- attr(*, "names")= chr [1:5] "q5" "q25" "q50" "q75" ...
#>   .. ..$ prob_exceeds  : Named num [1:2] 0.481 0.163
#>   .. .. ..- attr(*, "names")= chr [1:2] "prob_gt_0.5" "prob_gt_0.9"
#>   .. ..$ dominance_risk: chr "high"
#>   ..$ coclustering:List of 4
#>   .. ..$ mean          : num 0.501
#>   .. ..$ var           : num 0.065
#>   .. ..$ sd            : num 0.255
#>   .. ..$ interpretation: chr "High prior co-clustering: most unit pairs expected in same cluster"
#>   ..$ warnings    : chr [1:2] "HIGH DOMINANCE RISK: P(w1 > 0.5) = 48.1% exceeds 40%" "NEAR-DEGENERATE RISK: P(w1 > 0.9) = 16.3% exceeds 15%"
#>   ..- attr(*, "class")= chr "DPprior_diagnostics"
#>  $ trace        :'data.frame':   6 obs. of  8 variables:
#>   ..$ iter    : int [1:6] 1 2 3 4 5 6
#>   ..$ a       : num [1:6] 4 1.18 1.84 2.03 2.04 ...
#>   ..$ b       : num [1:6] 3.912 0.912 1.455 1.6 1.605 ...
#>   ..$ M1      : num [1:6] 4.46 4.91 4.97 5 5 ...
#>   ..$ V       : num [1:6] 4.78 10.85 8.4 8.01 8 ...
#>   ..$ residual: num [1:6] 3.26 2.86 4.00e-01 1.33e-02 2.08e-05 ...
#>   ..$ step    : num [1:6] 1 1 1 1 1 NA
#>   ..$ det_Jlog: num [1:6] -5.3 -21.6 -15.3 -14.3 -14.3 ...
#>  $ dual_anchor  : NULL
#>  - attr(*, "class")= chr "summary.DPprior_fit"
```

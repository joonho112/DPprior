# Verify Dual-Anchor Module

Verify Dual-Anchor Module

## Usage

``` r
verify_dual_anchor(verbose = TRUE)
```

## Arguments

- verbose:

  Logical; print detailed output.

## Value

Invisible TRUE if all tests pass.

## Examples

``` r
verify_dual_anchor(verbose = TRUE)
#> === Dual-Anchor Verification Tests ===
#> 
#> Test 1: ABSOLUTE loss barely moves (expected)
#>   Before: 0.4815, After: 0.4791, Reduction: 0.5%
#>   PASS (confirmed ABSOLUTE loss is inadequate)
#> 
#> Test 2: RELATIVE loss gives meaningful reduction
#>   Before: 0.4815, After: 0.4379, Reduction: 9.0%
#>   PASS
#> 
#> Test 3: ADAPTIVE loss is more aggressive
#>   Before: 0.4815, After: 0.3427, Reduction: 28.8%
#>   PASS
#> 
#> Test 4: lambda = 1 recovers K-only
#>   PASS
#> 
#> Test 5: lambda = 0 achieves weight target
#>   Target: 0.30, Achieved: 0.3000
#>   PASS
#> 
#> Test 6: Trade-off curve
#>   lambda w1_prob_gt_50     mu_K     var_K
#> 1    0.0     0.3000000 7.116329 11.903303
#> 2    0.5     0.4379077 5.380060  7.930247
#> 3    1.0     0.4814780 5.000000  8.000000
#>   PASS
#> 
#> === All Tests Passed ===
#> 
#> Summary:
#>   ABSOLUTE (broken): 0.5% reduction
#>   RELATIVE (fixed):  9.0% reduction
#>   ADAPTIVE (best):   28.8% reduction
```

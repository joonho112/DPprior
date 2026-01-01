# Verify KL Divergence Properties

Runs verification tests on KL divergence computations.

## Usage

``` r
verify_kl_divergence(verbose = TRUE)
```

## Arguments

- verbose:

  Logical; if TRUE, print detailed results.

## Value

Logical; TRUE if all tests pass.

## Examples

``` r
verify_kl_divergence()
#> ============================================================
#> KL Divergence Verification
#> ============================================================
#> 
#> Test 1: KL(p || p) = 0
#>   KL = 0.00e+00 [PASS]
#> 
#> Test 2: KL(p || q) >= 0
#>   KL(p || q) = 0.0305 [PASS]
#> 
#> Test 3: KL = 0 for matching (a, b)
#>   KL(p_{2.0,1.0} || p_{2.0,1.0}) = 0.00e+00 [PASS]
#> 
#> Test 4: KL > 0 for non-matching (a, b)
#>   KL = 0.3593 [PASS]
#> 
#> ============================================================
#> Overall: ALL TESTS PASSED
#> ============================================================
```

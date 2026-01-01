# Verify S3 Methods Module

Runs comprehensive verification tests for the S3 methods module.

## Usage

``` r
verify_s3_methods(verbose = TRUE)
```

## Arguments

- verbose:

  Logical; if TRUE, print detailed test output.

## Value

Invisibly returns TRUE if all tests pass.

## Examples

``` r
verify_s3_methods()
#> ============================================================ 
#> Module 17: S3 Methods - Verification Suite
#> ============================================================ 
#> 
#> [Test 1] Create mock DPprior_fit object
#> -------------------------------------------------- 
#>   Mock object created: PASS
#> 
#> [Test 2] print.DPprior_fit()
#> -------------------------------------------------- 
#>   print() works: PASS
#> 
#> [Test 3] summary.DPprior_fit()
#> -------------------------------------------------- 
#>   summary() works: PASS
#> 
#> [Test 4] print.summary.DPprior_fit()
#> -------------------------------------------------- 
#>   print(summary()) works: PASS
#> 
#> [Test 5] .dpprior_is_dual() helper (stricter check)
#> -------------------------------------------------- 
#>   .dpprior_is_dual() works (strict): PASS
#> 
#> [Test 6] var_K projection display
#> -------------------------------------------------- 
#>   var_K projection shown: PASS
#> 
#> [Test 7] summary with diagnostics parameter
#> -------------------------------------------------- 
#>   diagnostics parameter works: PASS
#> 
#> [Test 8] Output without diagnostics
#> -------------------------------------------------- 
#>   Works without diagnostics: PASS
#> 
#> ============================================================
#> All S3 methods tests PASSED!
#> ============================================================
```

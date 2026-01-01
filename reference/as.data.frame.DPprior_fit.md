# Coerce DPprior_fit to Data Frame

Coerce DPprior_fit to Data Frame

## Usage

``` r
# S3 method for class 'DPprior_fit'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  A `DPprior_fit` object.

- row.names:

  Optional row names.

- optional:

  Logical; if TRUE, avoid setting names.

- ...:

  Additional arguments (ignored).

## Value

A data frame with one row containing the fit results.

## Examples

``` r
fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)
as.data.frame(fit)
#>   method  status a        b  J mu_K var_K mean_alpha cv_alpha scaling converged
#> 1     A1 success 4 3.912023 50    5     8   1.022489      0.5     log      TRUE
#>   iterations
#> 1          0
```

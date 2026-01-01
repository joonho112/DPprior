# Discretize Chi-Square to K_J Support

Converts a (possibly scaled) chi-square distribution into a discrete PMF
on \\\\1, \dots, J\\\\ using continuity-corrected binning: \$\$p(k) =
P(k - 0.5 \< X \le k + 0.5), \quad k = 1, \dots, J\$\$ followed by
renormalization.

## Usage

``` r
discretize_chisq(J, df, scale = 1)
```

## Arguments

- J:

  Integer; maximum value (support upper bound).

- df:

  Numeric; degrees of freedom.

- scale:

  Numeric; scale parameter (default 1). If `scale != 1`, assumes \\X =
  scale \cdot \chi^2\_{df}\\.

## Value

Numeric vector of length J; PMF on \\\\1, \dots, J\\\\.

## Details

For a scaled chi-square distribution \\Y = scale \cdot X\\ where \\X
\sim \chi^2\_{df}\\:

- \\E\[Y\] = scale \cdot df\\

- \\Var\[Y\] = scale^2 \cdot 2 \cdot df\\

**Matching target moments:** Given target mean \\\mu_K\\ and variance
\\\sigma^2_K\\:

- \\scale = \sigma^2_K / (2 \mu_K)\\

- \\df = 2 \mu_K^2 / \sigma^2_K\\

## See also

[`DPprior_a2_kl`](https://joonho112.github.io/DPprior/reference/DPprior_a2_kl.md)

## Examples

``` r
# Chi-square with target moments mu=5, var=8
mu_K <- 5
var_K <- 8
scale <- var_K / (2 * mu_K)
df <- 2 * mu_K^2 / var_K
pmf <- discretize_chisq(50, df = df, scale = scale)

# Verify moments
k_vals <- 1:50
sum(k_vals * pmf)  # ~5
#> [1] 5.01509
```

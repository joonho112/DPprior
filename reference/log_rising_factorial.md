# Log Rising Factorial (Pochhammer Symbol)

Computes the logarithm of the rising factorial (Pochhammer symbol)
\\(\alpha)\_J = \alpha(\alpha+1)\cdots(\alpha+J-1)\\.

## Usage

``` r
log_rising_factorial(alpha, J)
```

## Arguments

- alpha:

  Numeric; concentration parameter (must be positive scalar).

- J:

  Integer; number of terms (must be \>= 1).

## Value

Numeric; \\\log(\alpha)\_J\\.

## Details

Uses the identity: \$\$(\alpha)\_J =
\frac{\Gamma(\alpha+J)}{\Gamma(\alpha)}\$\$

Therefore: \$\$\log(\alpha)\_J = \log\Gamma(\alpha+J) -
\log\Gamma(\alpha)\$\$

This is numerically stable for all \\\alpha \> 0\\ and \\J \geq 1\\.

## See also

[`lgamma`](https://rdrr.io/r/base/Special.html) for the log-gamma
function

## Examples

``` r
if (FALSE) { # \dontrun{
# (2)_3 = 2 * 3 * 4 = 24
exp(log_rising_factorial(2, 3))

# (1)_5 = 5!
exp(log_rising_factorial(1, 5))
factorial(5)

} # }
```

# Derivative of Conditional Mean w.r.t. Alpha

Computes \\\frac{d}{d\alpha} \mathbb{E}\[K_J \mid \alpha\]\\. This is
used internally when building Jacobians for Newton-type solvers.

## Usage

``` r
dmean_dalpha(J, alpha)
```

## Arguments

- J:

  Sample size (integer \>= 1).

- alpha:

  Concentration parameter (positive numeric, vectorized).

## Value

Numeric vector of derivatives (same length as `alpha`).

## Details

Closed form: \$\$\frac{d}{d\alpha}\mu_J(\alpha) =
\\\psi(\alpha+J)-\psi(\alpha)\\ +
\alpha\\\psi_1(\alpha+J)-\psi_1(\alpha)\\\$\$

This derivative is always positive for \\\alpha \> 0\\, confirming that
\\E\[K_J \| \alpha\]\\ is strictly increasing in \\\alpha\\.

## Examples

``` r
if (FALSE) { # \dontrun{
dmean_dalpha(50, 2.0)

# Verify with finite difference
J <- 50; alpha <- 2.0; eps <- 1e-6
fd <- (mean_K_given_alpha(J, alpha + eps) -
       mean_K_given_alpha(J, alpha - eps)) / (2 * eps)
abs(fd - dmean_dalpha(J, alpha)) < 1e-5  # TRUE

} # }
```

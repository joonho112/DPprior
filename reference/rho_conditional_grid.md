# Compute rho Conditional Moments on alpha Grid

Computes conditional mean and variance of rho on a grid of alpha values.
Useful for visualization of how rho varies with alpha.

## Usage

``` r
rho_conditional_grid(alpha_grid = seq(0.1, 10, length.out = 100))
```

## Arguments

- alpha_grid:

  Numeric vector; grid of alpha values. Default is
  `seq(0.1, 10, length.out = 100)`.

## Value

A data frame with columns:

- alpha:

  Grid points

- mean:

  E(rho \| alpha)

- var:

  Var(rho \| alpha)

- sd:

  SD(rho \| alpha)

## Examples

``` r
df <- rho_conditional_grid()
plot(df$alpha, df$mean, type = "l",
     xlab = expression(alpha), ylab = expression(E(rho)))

```

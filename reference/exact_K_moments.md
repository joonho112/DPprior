# Exact Marginal Moments of K_J under Gamma Prior

Computes the exact marginal mean \\E\[K_J\]\\ and variance \\Var(K_J)\\
when the DP concentration parameter follows a Gamma(a, b) prior.

## Usage

``` r
exact_K_moments(J, a, b, M = .QUAD_NODES_DEFAULT)
```

## Arguments

- J:

  Integer; sample size (positive integer \>= 1).

- a:

  Numeric; shape parameter of Gamma prior (\> 0).

- b:

  Numeric; rate parameter of Gamma prior (\> 0).

- M:

  Integer; number of quadrature nodes (default: 80).

## Value

A named list with components:

- `mean`:

  Marginal mean \\E\[K_J \| a, b\]\\

- `var`:

  Marginal variance \\Var(K_J \| a, b)\\

- `sd`:

  Marginal standard deviation

- `cv`:

  Coefficient of variation (sd/mean)

## Details

Uses Gauss-Laguerre quadrature to numerically evaluate: \$\$M_1(a,b) =
E\_{\alpha \sim \Gamma(a,b)}\[\mu_J(\alpha)\]\$\$ \$\$V(a,b) =
E\[v_J(\alpha)\] + E\[\mu_J(\alpha)^2\] - M_1^2\$\$

where \\\mu_J(\alpha)\\ and \\v_J(\alpha)\\ are the conditional mean and
variance from Module 03.

The Law of Total Variance decomposes the marginal variance into:

- Within-alpha variance: \\E\[v_J(\alpha)\]\\

- Between-alpha variance: \\Var(\mu_J(\alpha))\\

**Key properties:**

- The mean is bounded: \\1 \leq E\[K_J\] \leq J\\

- Despite conditional underdispersion (\\v_J(\alpha) \<
  \mu_J(\alpha)\\), the marginal distribution is typically overdispersed

- Marginal variance exceeds conditional variance at \\\alpha =
  E\[\alpha\]\\

## References

Antoniak, C. E. (1974). Mixtures of Dirichlet Processes. *The Annals of
Statistics*, 2(6), 1152-1174.

## See also

[`K_moments`](https://joonho112.github.io/DPprior/reference/K_moments.md)
for convenience wrapper,
[`mean_K_given_alpha`](https://joonho112.github.io/DPprior/reference/mean_K_given_alpha.md),
[`var_K_given_alpha`](https://joonho112.github.io/DPprior/reference/var_K_given_alpha.md)

## Examples

``` r
# Example from RN-01: J=50, Gamma(1.5, 0.5) prior
result <- exact_K_moments(50, 1.5, 0.5)
print(result)
#> $mean
#> [1] 8.355487
#> 
#> $var
#> [1] 22.76895
#> 
#> $sd
#> [1] 4.771682
#> 
#> $cv
#> [1] 0.5710837
#> 

# Compare with conditional variance at E[alpha] = 3
cond_var <- var_K_given_alpha(50, 3.0)
result$var > cond_var  # TRUE: marginal > conditional
#> [1] TRUE

# Verify mean bounds
1 <= result$mean && result$mean <= 50  # TRUE
#> [1] TRUE
```

# Compute Log Stirling Numbers (First Kind, Unsigned)

Computes the logarithm of unsigned Stirling numbers of the first kind
\\\|s(J,k)\|\\ for all \\J\\ from 0 to `J_max` and \\k\\ from 0 to
\\J\\. Uses log-space recursion to prevent numerical overflow.

## Usage

``` r
compute_log_stirling(J_max)
```

## Arguments

- J_max:

  Maximum J value (non-negative integer, at most 500).

## Value

A lower triangular matrix of dimension `(J_max+1) x (J_max+1)`. Entry
`[J+1, k+1]` contains \\\log\|s(J,k)\|\\ (using R's 1-based indexing).
Invalid entries (k \> J or k \< 1 when J \>= 1) contain `-Inf`.

## Details

The unsigned Stirling numbers of the first kind \\\|s(J,k)\|\\ count the
number of permutations of \\J\\ elements into exactly \\k\\ disjoint
cycles.

Uses the log-space recurrence relation: \$\$L\_{J,k} =
\text{logsumexp}(L\_{J-1,k-1}, \log(J-1) + L\_{J-1,k})\$\$

where \\L\_{J,k} = \log\|s(J,k)\|\\.

Boundary conditions:

- \\\|s(0,0)\| = 1\\ \\\rightarrow\\ `L[1,1] = 0`

- \\\|s(J,0)\| = 0\\ for \\J \geq 1\\ \\\rightarrow\\ `L[J+1,1] = -Inf`

- \\\|s(J,J)\| = 1\\ for all \\J\\ \\\rightarrow\\ `L[J+1,J+1] = 0`

Complexity: O(`J_max`^2) time and space. Results should be cached for
repeated use within a session.

## References

Antoniak, C. E. (1974). Mixtures of Dirichlet Processes with
Applications to Bayesian Nonparametric Problems. *The Annals of
Statistics*, 2(6), 1152-1174.

## See also

[`get_log_stirling`](https://joonho112.github.io/DPprior/reference/get_log_stirling.md)
for safe accessor with bounds checking

## Examples

``` r
# Compute Stirling numbers up to J=10
logS <- compute_log_stirling(10)

# Access |s(4,2)| = 11
exp(logS[5, 3])
#> [1] 11

# Access |s(5,3)| = 35
exp(logS[6, 4])
#> [1] 35
```

# PLN RNG

Random generation for the PLN model with latent mean equal to mu, latent
covariance matrix equal to Sigma and average depths (sum of counts in a
sample) equal to depths

## Usage

``` r
rPLN(
  n = 10,
  mu = rep(0, ncol(Sigma)),
  Sigma = diag(1, 5, 5),
  depths = rep(10000, n)
)
```

## Arguments

- n:

  the sample size

- mu:

  vectors of means of the latent variable

- Sigma:

  covariance matrix of the latent variable

- depths:

  Numeric vector of target depths. The first is recycled if there are
  not `n` values

## Value

a n \* p count matrix, with row-sums close to depths, with an attribute
"offsets" corresponding to the true generated offsets (in log-scale).

## Details

The default value for mu and Sigma assume equal abundances and no
correlation between the different species.

## Examples

``` r
## 10 samples of 5 species with equal abundances, no covariance and target depths of 10,000
rPLN()
#>       Y1   Y2   Y3   Y4    Y5
#> S1  1632 1771 1046 2086   515
#> S2  6213  548  551 1324   908
#> S3  8941  828  354 9925   213
#> S4  1563   46 2187 2669  3339
#> S5   537  515  643 1200 11571
#> S6   381 4359 1812 2993  6939
#> S7  3350 1384 3040 1247  1397
#> S8  2902  273 1231 6662  2663
#> S9   188 1251  515 2329  4458
#> S10 1405  161 4190  660   926
#> attr(,"offsets")
#>           [,1]     [,2]     [,3]     [,4]     [,5]
#>  [1,] 7.100902 7.100902 7.100902 7.100902 7.100902
#>  [2,] 7.100902 7.100902 7.100902 7.100902 7.100902
#>  [3,] 7.100902 7.100902 7.100902 7.100902 7.100902
#>  [4,] 7.100902 7.100902 7.100902 7.100902 7.100902
#>  [5,] 7.100902 7.100902 7.100902 7.100902 7.100902
#>  [6,] 7.100902 7.100902 7.100902 7.100902 7.100902
#>  [7,] 7.100902 7.100902 7.100902 7.100902 7.100902
#>  [8,] 7.100902 7.100902 7.100902 7.100902 7.100902
#>  [9,] 7.100902 7.100902 7.100902 7.100902 7.100902
#> [10,] 7.100902 7.100902 7.100902 7.100902 7.100902
## 2 samples of 10 highly correlated species with target depths 1,000 and 100,000
## very different abundances
mu <- rep(c(1, -1), each = 5)
Sigma <- matrix(0.8, 10, 10); diag(Sigma) <- 1
rPLN(n=2, mu = mu, Sigma = Sigma, depths = c(1e3, 1e5))
#>       Y1   Y2    Y3    Y4    Y5   Y6   Y7   Y8   Y9  Y10
#> S1   336  134   228   274   238   40   67   22   22   47
#> S2 10876 8772 10003 20058 26924 2073 4011 1163 2599 1681
#> attr(,"offsets")
#>          [,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]     [,8]
#> [1,] 3.671389 3.671389 3.671389 3.671389 3.671389 3.671389 3.671389 3.671389
#> [2,] 8.276560 8.276560 8.276560 8.276560 8.276560 8.276560 8.276560 8.276560
#>          [,9]    [,10]
#> [1,] 3.671389 3.671389
#> [2,] 8.276560 8.276560
```

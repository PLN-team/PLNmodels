# Extract model coefficients

Extracts model coefficients from objects returned by
[`ZIPLN()`](https://pln-team.github.io/PLNmodels/reference/ZIPLN.md) and
its variants

## Usage

``` r
# S3 method for class 'ZIPLNfit'
coef(object, type = c("count", "zero", "precision", "covariance"), ...)
```

## Arguments

- object:

  an R6 object with class
  [`ZIPLNfit`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit.md)

- type:

  type of parameter that should be extracted. Either "count" (default)
  for \\B\\, "zero" for \\B0\\, "precision" for \\\Omega\\, "covariance"
  for \\\Sigma\\

- ...:

  additional parameters for S3 compatibility. Not used

## Value

A matrix of coefficients extracted from the ZIPLNfit model.

## See also

[`sigma.ZIPLNfit()`](https://pln-team.github.io/PLNmodels/reference/sigma.ZIPLNfit.md)

## Examples

``` r
data(scRNA)
# data subsample: only 100 random cell and the 50 most varying transcript
subset <- sample.int(nrow(scRNA), 100)
myPLN  <- ZIPLN(counts[, 1:50] ~ 1 + offset(log(total_counts)), subset = subset, data = scRNA)
#> 
#>  Initialization...
#>  Adjusting a ZI-PLN model with full covariance model and single specific parameter(s) in Zero inflation component.
#>  DONE!
```

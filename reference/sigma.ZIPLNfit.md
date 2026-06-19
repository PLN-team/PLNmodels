# Extract variance-covariance of residuals 'Sigma'

Extract the variance-covariance matrix of the residuals, usually noted
\\\Sigma\\ in ZIPLN models.

## Usage

``` r
# S3 method for class 'ZIPLNfit'
sigma(object, ...)
```

## Arguments

- object:

  an R6 object with class
  [`ZIPLNfit`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit.md)

- ...:

  additional parameters for S3 compatibility. Not used

## Value

A semi definite positive matrix of size p, assuming there are p species
in the model.

## See also

[`coef.ZIPLNfit()`](https://pln-team.github.io/PLNmodels/reference/coef.ZIPLNfit.md)

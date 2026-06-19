# Extract log-likelihood of a fitted PLN model

Returns the variational lower bound of the log-likelihood as a
`"logLik"` object, compatible with
[`stats::AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`stats::BIC()`](https://rdrr.io/r/stats/AIC.html).

## Usage

``` r
# S3 method for class 'PLNfit'
logLik(object, ...)
```

## Arguments

- object:

  an R6 object with class
  [`PLNfit`](https://pln-team.github.io/PLNmodels/reference/PLNfit.md)

- ...:

  additional parameters for S3 compatibility. Not used

## Value

An object of class `"logLik"`. The numeric value is the variational
ELBO. Attributes `df` and `nobs` hold the number of parameters and
observations.

## Examples

``` r
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
model <- PLN(Abundance ~ 1, data = trichoptera)
#> 
#>  Initialization...
#>  Adjusting a full covariance PLN model with nlopt optimizer
#>  Post-treatments...
#>  DONE!
logLik(model)
#> 'log Lik.' -1129.624 (df=170)
```

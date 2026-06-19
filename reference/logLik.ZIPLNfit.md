# Extract log-likelihood of a fitted ZIPLN model

Returns the variational lower bound of the log-likelihood as a
`"logLik"` object, compatible with
[`stats::AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`stats::BIC()`](https://rdrr.io/r/stats/AIC.html).

## Usage

``` r
# S3 method for class 'ZIPLNfit'
logLik(object, ...)
```

## Arguments

- object:

  an R6 object with class
  [`ZIPLNfit`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit.md)

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
model <- ZIPLN(Abundance ~ 1, data = trichoptera)
#> 
#>  Initialization...
#>  Adjusting a ZI-PLN model with full covariance model and single specific parameter(s) in Zero inflation component.
#>  DONE!
logLik(model)
#> 'log Lik.' -1139.296 (df=171)
```

# Extract the regularization path of a PLNnetwork fit

Extract the regularization path of a PLNnetwork fit

## Usage

``` r
coefficient_path(Robject, precision = TRUE, corr = TRUE)
```

## Arguments

- Robject:

  an object with class
  [`Networkfamily`](https://pln-team.github.io/PLNmodels/reference/Networkfamily.md),
  i.e. an output from
  [`PLNnetwork()`](https://pln-team.github.io/PLNmodels/reference/PLNnetwork.md)

- precision:

  a logical, should the coefficients of the precision matrix Omega or
  the covariance matrix Sigma be sent back. Default is `TRUE`.

- corr:

  a logical, should the correlation (partial in case `precision = TRUE`)
  be sent back. Default is `TRUE`.

## Value

Sends back a tibble/data.frame.

## Examples

``` r
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
fits <- PLNnetwork(Abundance ~ 1, data = trichoptera)
#> 
#>  Initialization...
#>  Adjusting 30 PLN with sparse inverse covariance estimation
#>  Joint optimization alternating gradient descent and graphical-lasso
#>  sparsifying penalty = 4.479926  sparsifying penalty = 4.137977  sparsifying penalty = 3.822129  sparsifying penalty = 3.530389  sparsifying penalty = 3.260917  sparsifying penalty = 3.012014  sparsifying penalty = 2.78211   sparsifying penalty = 2.569754  sparsifying penalty = 2.373607  sparsifying penalty = 2.192431  sparsifying penalty = 2.025085  sparsifying penalty = 1.870512  sparsifying penalty = 1.727737  sparsifying penalty = 1.595861  sparsifying penalty = 1.47405   sparsifying penalty = 1.361537  sparsifying penalty = 1.257612  sparsifying penalty = 1.16162   sparsifying penalty = 1.072954  sparsifying penalty = 0.9910565     sparsifying penalty = 0.91541   sparsifying penalty = 0.8455375     sparsifying penalty = 0.7809984     sparsifying penalty = 0.7213854     sparsifying penalty = 0.6663227     sparsifying penalty = 0.6154629     sparsifying penalty = 0.5684851     sparsifying penalty = 0.5250931     sparsifying penalty = 0.4850132     sparsifying penalty = 0.4479926 
#>  Post-treatments
#>  DONE!
head(coefficient_path(fits))
#>   Node1 Node2 Coeff  Penalty    Edge
#> 1   Aga   Che     0 4.479926 Aga|Che
#> 2   Ath   Che     0 4.479926 Ath|Che
#> 3   Cea   Che     0 4.479926 Cea|Che
#> 4   Ced   Che     0 4.479926 Ced|Che
#> 5   All   Che     0 4.479926 All|Che
#> 6   Che   Hyc     0 4.479926 Che|Hyc
```

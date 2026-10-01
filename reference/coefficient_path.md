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
#>  sparsifying penalty = 2.847172  sparsifying penalty = 2.62985   sparsifying penalty = 2.429116  sparsifying penalty = 2.243704  sparsifying penalty = 2.072444  sparsifying penalty = 1.914256  sparsifying penalty = 1.768142  sparsifying penalty = 1.633182  sparsifying penalty = 1.508522  sparsifying penalty = 1.393378  sparsifying penalty = 1.287023  sparsifying penalty = 1.188785  sparsifying penalty = 1.098046  sparsifying penalty = 1.014233  sparsifying penalty = 0.9368178     sparsifying penalty = 0.8653113     sparsifying penalty = 0.7992629     sparsifying penalty = 0.7382558     sparsifying penalty = 0.6819054     sparsifying penalty = 0.6298561     sparsifying penalty = 0.5817797     sparsifying penalty = 0.537373  sparsifying penalty = 0.4963558     sparsifying penalty = 0.4584694     sparsifying penalty = 0.4234748     sparsifying penalty = 0.3911513     sparsifying penalty = 0.3612951     sparsifying penalty = 0.3337177     sparsifying penalty = 0.3082453     sparsifying penalty = 0.2847172 
#>  Post-treatments
#>  DONE!
head(coefficient_path(fits))
#>   Node1 Node2 Coeff  Penalty    Edge
#> 1   Aga   Che     0 2.847172 Aga|Che
#> 2   Ath   Che     0 2.847172 Ath|Che
#> 3   Cea   Che     0 2.847172 Cea|Che
#> 4   Ced   Che     0 2.847172 Ced|Che
#> 5   All   Che     0 2.847172 All|Che
#> 6   Che   Hyc     0 2.847172 Che|Hyc
```

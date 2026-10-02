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
#>  sparsifying penalty = 0.5246191     sparsifying penalty = 0.4845754     sparsifying penalty = 0.4475881     sparsifying penalty = 0.4134241     sparsifying penalty = 0.3818678     sparsifying penalty = 0.3527202     sparsifying penalty = 0.3257973     sparsifying penalty = 0.3009295     sparsifying penalty = 0.2779598     sparsifying penalty = 0.2567434     sparsifying penalty = 0.2371464     sparsifying penalty = 0.2190452     sparsifying penalty = 0.2023257     sparsifying penalty = 0.1868823     sparsifying penalty = 0.1726178     sparsifying penalty = 0.159442  sparsifying penalty = 0.1472719     sparsifying penalty = 0.1360308     sparsifying penalty = 0.1256477     sparsifying penalty = 0.1160571     sparsifying penalty = 0.1071986     sparsifying penalty = 0.09901618    sparsifying penalty = 0.09145836    sparsifying penalty = 0.08447742    sparsifying penalty = 0.07802933    sparsifying penalty = 0.07207342    sparsifying penalty = 0.06657212    sparsifying penalty = 0.06149072    sparsifying penalty = 0.05679719    sparsifying penalty = 0.05246191 
#>  Post-treatments
#>  DONE!
head(coefficient_path(fits))
#>   Node1 Node2 Coeff   Penalty    Edge
#> 1   Aga   Che     0 0.5246191 Aga|Che
#> 2   Ath   Che     0 0.5246191 Ath|Che
#> 3   Cea   Che     0 0.5246191 Cea|Che
#> 4   Ced   Che     0 0.5246191 Ced|Che
#> 5   All   Che     0 0.5246191 All|Che
#> 6   Che   Hyc     0 0.5246191 Che|Hyc
```

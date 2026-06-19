# Best model extraction from a collection of models

Best model extraction from a collection of models

## Usage

``` r
# S3 method for class 'PLNPCAfamily'
getBestModel(Robject, crit = c("ICL", "BIC"), ...)

getBestModel(Robject, crit, ...)

# S3 method for class 'PLNmixturefamily'
getBestModel(Robject, crit = c("ICL", "BIC"), ...)

# S3 method for class 'Networkfamily'
getBestModel(Robject, crit = c("BIC", "EBIC", "StARS"), ...)

# S3 method for class 'PLNnetworkfamily'
getBestModel(Robject, crit = c("BIC", "EBIC", "StARS"), ...)

# S3 method for class 'ZIPLNnetworkfamily'
getBestModel(Robject, crit = c("BIC", "EBIC", "StARS"), ...)
```

## Arguments

- Robject:

  an object with class PLNPCAfamilly ot PLNnetworkfamily

- crit:

  a character for the criterion used to performed the selection. Either
  "BIC", "ICL", "EBIC", "StARS", "R_squared". Default is `ICL` for
  `PLNPCA`, and `BIC` for `PLNnetwork`. If StARS (Stability Approach to
  Regularization Selection) is chosen and stability selection was not
  yet performed, the function will call the method
  [`stability_selection()`](https://pln-team.github.io/PLNmodels/reference/stability_selection.md)
  with default argument.

- ...:

  additional parameters for StARS criterion (only for `PLNnetwork`).
  `stability`, a scalar indicating the target stability (= 1 - 2 beta)
  at which the network is selected. Default is `0.9`.

## Value

Send back an object with class
[`PLNPCAfit`](https://pln-team.github.io/PLNmodels/reference/PLNPCAfit.md)
or
[`PLNnetworkfit`](https://pln-team.github.io/PLNmodels/reference/PLNnetworkfit.md)

## Methods (by class)

- `getBestModel(PLNPCAfamily)`: Model extraction for
  [`PLNPCAfamily`](https://pln-team.github.io/PLNmodels/reference/PLNPCAfamily.md)

- `getBestModel(PLNmixturefamily)`: Model extraction for
  [`PLNmixturefamily`](https://pln-team.github.io/PLNmodels/reference/PLNmixturefamily.md)

- `getBestModel(Networkfamily)`: Model extraction for
  [`PLNnetworkfamily`](https://pln-team.github.io/PLNmodels/reference/PLNnetworkfamily.md)
  or
  [`ZIPLNnetworkfamily`](https://pln-team.github.io/PLNmodels/reference/ZIPLNnetworkfamily.md)

- `getBestModel(PLNnetworkfamily)`: Model extraction for
  [`PLNnetworkfamily`](https://pln-team.github.io/PLNmodels/reference/PLNnetworkfamily.md)

- `getBestModel(ZIPLNnetworkfamily)`: Model extraction for
  [`ZIPLNnetworkfamily`](https://pln-team.github.io/PLNmodels/reference/ZIPLNnetworkfamily.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
myPCA <- PLNPCA(Abundance ~ 1 + offset(log(Offset)), data = trichoptera, ranks = 1:4)
myModel <- getBestModel(myPCA)
} # }
```

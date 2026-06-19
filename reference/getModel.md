# Model extraction from a collection of models

Model extraction from a collection of models

## Usage

``` r
# S3 method for class 'PLNPCAfamily'
getModel(Robject, var, index = NULL)

getModel(Robject, var, index)

# S3 method for class 'PLNmixturefamily'
getModel(Robject, var, index = NULL)

# S3 method for class 'Networkfamily'
getModel(Robject, var, index = NULL)

# S3 method for class 'PLNnetworkfamily'
getModel(Robject, var, index = NULL)

# S3 method for class 'ZIPLNnetworkfamily'
getModel(Robject, var, index = NULL)
```

## Arguments

- Robject:

  an R6 object with class
  [`PLNPCAfamily`](https://pln-team.github.io/PLNmodels/reference/PLNPCAfamily.md)
  or
  [`PLNnetworkfamily`](https://pln-team.github.io/PLNmodels/reference/PLNnetworkfamily.md)

- var:

  value of the parameter (`rank` for PLNPCA, `sparsity` for PLNnetwork)
  that identifies the model to be extracted from the collection. If no
  exact match is found, the model with closest parameter value is
  returned with a warning.

- index:

  Integer index of the model to be returned. Only the first value is
  taken into account.

## Value

Sends back an object with class
[`PLNPCAfit`](https://pln-team.github.io/PLNmodels/reference/PLNPCAfit.md)
or
[`PLNnetworkfit`](https://pln-team.github.io/PLNmodels/reference/PLNnetworkfit.md).

## Methods (by class)

- `getModel(PLNPCAfamily)`: Model extraction for
  [`PLNPCAfamily`](https://pln-team.github.io/PLNmodels/reference/PLNPCAfamily.md)

- `getModel(PLNmixturefamily)`: Model extraction for
  [`PLNmixturefamily`](https://pln-team.github.io/PLNmodels/reference/PLNmixturefamily.md)

- `getModel(Networkfamily)`: Model extraction for
  [`PLNnetworkfamily`](https://pln-team.github.io/PLNmodels/reference/PLNnetworkfamily.md)
  or
  [`ZIPLNnetworkfamily`](https://pln-team.github.io/PLNmodels/reference/ZIPLNnetworkfamily.md)

- `getModel(PLNnetworkfamily)`: Model extraction for
  [`PLNnetworkfamily`](https://pln-team.github.io/PLNmodels/reference/PLNnetworkfamily.md)

- `getModel(ZIPLNnetworkfamily)`: Model extraction for
  [`ZIPLNnetworkfamily`](https://pln-team.github.io/PLNmodels/reference/ZIPLNnetworkfamily.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
myPCA <- PLNPCA(Abundance ~ 1 + offset(log(Offset)), data = trichoptera, ranks = 1:5)
myModel <- getModel(myPCA, 2)
} # }
```

# An R6 Class to represent a collection of PLNPCAfit

The function
[`PLNPCA()`](https://pln-team.github.io/PLNmodels/reference/PLNPCA.md)
produces an instance of this class.

This class comes with a set of methods, some of them being useful for
the user: See the documentation for
[`getBestModel()`](https://pln-team.github.io/PLNmodels/reference/getBestModel.md),
[`getModel()`](https://pln-team.github.io/PLNmodels/reference/getModel.md)
and
[`plot()`](https://pln-team.github.io/PLNmodels/reference/plot.PLNPCAfamily.md).

## See also

The function
[`PLNPCA()`](https://pln-team.github.io/PLNmodels/reference/PLNPCA.md),
the class
[`PLNPCAfit()`](https://pln-team.github.io/PLNmodels/reference/PLNPCAfit.md)

## Super class

[`PLNfamily`](https://pln-team.github.io/PLNmodels/reference/PLNfamily.md)
-\> `PLNPCAfamily`

## Active bindings

- `ranks`:

  the dimensions of the successively fitted models

## Methods

### Public methods

- [`PLNPCAfamily$new()`](#method-PLNPCAfamily-initialize)

- [`PLNPCAfamily$optimize()`](#method-PLNPCAfamily-optimize)

- [`PLNPCAfamily$getModel()`](#method-PLNPCAfamily-getModel)

- [`PLNPCAfamily$getBestModel()`](#method-PLNPCAfamily-getBestModel)

- [`PLNPCAfamily$plot()`](#method-PLNPCAfamily-plot)

- [`PLNPCAfamily$show()`](#method-PLNPCAfamily-show)

- [`PLNPCAfamily$clone()`](#method-PLNPCAfamily-clone)

Inherited methods

- [`PLNfamily$postTreatment()`](https://pln-team.github.io/PLNmodels/reference/PLNfamily.html#method-postTreatment)
- [`PLNfamily$print()`](https://pln-team.github.io/PLNmodels/reference/PLNfamily.html#method-print)

------------------------------------------------------------------------

### `PLNPCAfamily$new()`

Initialize all models in the collection. A single SVD of the residual
matrix `M - X*B` is computed once and shared across all ranks. `M` and
`B` come from either a user-provided
[`PLNfit`](https://pln-team.github.io/PLNmodels/reference/PLNfit.md)
inception or a fast LM on log-transformed counts (default, controlled by
`init_method`).

#### Usage

    PLNPCAfamily$new(
      ranks,
      responses,
      covariates,
      offsets,
      weights,
      formula,
      control
    )

#### Arguments

- `ranks`:

  the dimensions of the successively fitted models

- `responses`:

  the matrix of responses common to every models

- `covariates`:

  the matrix of covariates common to every models

- `offsets`:

  the matrix of offsets common to every models

- `weights`:

  the vector of observation weights

- `formula`:

  model formula used for fitting, extracted from the formula in the
  upper-level call

- `control`:

  list controlling the optimization and the model

------------------------------------------------------------------------

### `PLNPCAfamily$optimize()`

Call to the C++ optimizer on all models of the collection

#### Usage

    PLNPCAfamily$optimize(config)

#### Arguments

- `config`:

  list controlling the optimization.

------------------------------------------------------------------------

### `PLNPCAfamily$getModel()`

Extract model from collection and add "PCA" class for compatibility with
[`factoextra::fviz()`](https://rdrr.io/pkg/factoextra/man/fviz.html)

#### Usage

    PLNPCAfamily$getModel(var, index = NULL)

#### Arguments

- `var`:

  value of the parameter (rank for PLNPCA, sparsity for PLNnetwork) that
  identifies the model to be extracted from the collection. If no exact
  match is found, the model with closest parameter value is returned
  with a warning.

- `index`:

  Integer index of the model to be returned. Only the first value is
  taken into account.

#### Returns

a
[`PLNPCAfit`](https://pln-team.github.io/PLNmodels/reference/PLNPCAfit.md)
object

------------------------------------------------------------------------

### `PLNPCAfamily$getBestModel()`

Extract best model in the collection

#### Usage

    PLNPCAfamily$getBestModel(crit = c("ICL", "BIC"))

#### Arguments

- `crit`:

  a character for the criterion used to performed the selection. Either
  "ICL", "BIC". Default is `ICL`

#### Returns

a
[`PLNPCAfit`](https://pln-team.github.io/PLNmodels/reference/PLNPCAfit.md)
object

------------------------------------------------------------------------

### `PLNPCAfamily$plot()`

Lineplot of selected criteria for all models in the collection

#### Usage

    PLNPCAfamily$plot(criteria = c("loglik", "BIC", "ICL"), reverse = FALSE)

#### Arguments

- `criteria`:

  A valid model selection criteria for the collection of models. Any of
  "loglik", "BIC" or "ICL" (all).

- `reverse`:

  A logical indicating whether to plot the value of the criteria in the
  "natural" direction (loglik - penalty) or in the "reverse" direction
  (-2 loglik + penalty). Default to FALSE, i.e use the natural
  direction, on the same scale as the log-likelihood.

#### Returns

A
[`ggplot2::ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object

------------------------------------------------------------------------

### `PLNPCAfamily$show()`

User friendly print method

#### Usage

    PLNPCAfamily$show()

------------------------------------------------------------------------

### `PLNPCAfamily$clone()`

The objects of this class are cloneable with this method.

#### Usage

    PLNPCAfamily$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
myPCAs <- PLNPCA(Abundance ~ 1 + offset(log(Offset)), data = trichoptera, ranks = 1:5)
#> 
#>  Initialization...
#> 
#>  Adjusting 5 PLN models for PCA analysis.
#>   Rank approximation = 1      Rank approximation = 2      Rank approximation = 3      Rank approximation = 4      Rank approximation = 5 
#>  Post-treatments
#>  DONE!
class(myPCAs)
#> [1] "PLNPCAfamily" "PLNfamily"    "R6"          
```

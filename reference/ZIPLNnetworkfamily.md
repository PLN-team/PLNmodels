# An R6 Class to represent a collection of ZIPLNnetwork

The function
[`ZIPLNnetwork()`](https://pln-team.github.io/PLNmodels/reference/ZIPLNnetwork.md)
produces an instance of this class.

This class comes with a set of methods, some of them being useful for
the user: See the documentation for
[`getBestModel()`](https://pln-team.github.io/PLNmodels/reference/getBestModel.md),
[`getModel()`](https://pln-team.github.io/PLNmodels/reference/getModel.md)
and
[plot()](https://pln-team.github.io/PLNmodels/reference/plot.Networkfamily.md)

## See also

The function
[`ZIPLNnetwork()`](https://pln-team.github.io/PLNmodels/reference/ZIPLNnetwork.md),
the class
[`ZIPLNfit_sparse`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit_sparse.md)

## Super classes

[`PLNfamily`](https://pln-team.github.io/PLNmodels/reference/PLNfamily.md)
-\>
[`Networkfamily`](https://pln-team.github.io/PLNmodels/reference/Networkfamily.md)
-\> `ZIPLNnetworkfamily`

## Public fields

- `covariates0`:

  the matrix of covariates included in the ZI component

## Methods

### Public methods

- [`ZIPLNnetworkfamily$new()`](#method-ZIPLNnetworkfamily-initialize)

- [`ZIPLNnetworkfamily$stability_selection()`](#method-ZIPLNnetworkfamily-stability_selection)

- [`ZIPLNnetworkfamily$clone()`](#method-ZIPLNnetworkfamily-clone)

Inherited methods

- [`PLNfamily$getModel()`](https://pln-team.github.io/PLNmodels/reference/PLNfamily.html#method-getModel)
- [`PLNfamily$postTreatment()`](https://pln-team.github.io/PLNmodels/reference/PLNfamily.html#method-postTreatment)
- [`PLNfamily$print()`](https://pln-team.github.io/PLNmodels/reference/PLNfamily.html#method-print)
- [`Networkfamily$coefficient_path()`](https://pln-team.github.io/PLNmodels/reference/Networkfamily.html#method-coefficient_path)
- [`Networkfamily$getBestModel()`](https://pln-team.github.io/PLNmodels/reference/Networkfamily.html#method-getBestModel)
- [`Networkfamily$optimize()`](https://pln-team.github.io/PLNmodels/reference/Networkfamily.html#method-optimize)
- [`Networkfamily$plot()`](https://pln-team.github.io/PLNmodels/reference/Networkfamily.html#method-plot)
- [`Networkfamily$plot_objective()`](https://pln-team.github.io/PLNmodels/reference/Networkfamily.html#method-plot_objective)
- [`Networkfamily$plot_stars()`](https://pln-team.github.io/PLNmodels/reference/Networkfamily.html#method-plot_stars)
- [`Networkfamily$show()`](https://pln-team.github.io/PLNmodels/reference/Networkfamily.html#method-show)

------------------------------------------------------------------------

### `ZIPLNnetworkfamily$new()`

Initialize all models in the collection

#### Usage

    ZIPLNnetworkfamily$new(penalties, data, control)

#### Arguments

- `penalties`:

  a vector of positive real number controlling the level of sparsity of
  the underlying network.

- `data`:

  a named list used internally to carry the data matrices

- `control`:

  a list for controlling the optimization.

#### Returns

Update current
[`PLNnetworkfit`](https://pln-team.github.io/PLNmodels/reference/PLNnetworkfit.md)
with smart starting values

------------------------------------------------------------------------

### `ZIPLNnetworkfamily$stability_selection()`

Compute the stability path by stability selection

#### Usage

    ZIPLNnetworkfamily$stability_selection(
      subsamples = NULL,
      control = ZIPLNnetwork_param()
    )

#### Arguments

- `subsamples`:

  a list of vectors describing the subsamples. The number of vectors (or
  list length) determines the number of subsamples used in the stability
  selection. Automatically set to 20 subsamples with size `10*sqrt(n)`
  if `n >= 144` and `0.8*n` otherwise following Liu et al. (2010)
  recommendations.

- `control`:

  a list controlling the main optimization process in each call to
  [`PLNnetwork()`](https://pln-team.github.io/PLNmodels/reference/PLNnetwork.md).
  See
  [`ZIPLNnetwork()`](https://pln-team.github.io/PLNmodels/reference/ZIPLNnetwork.md)
  and
  [`ZIPLN_param()`](https://pln-team.github.io/PLNmodels/reference/ZIPLN_param.md)
  for details.

------------------------------------------------------------------------

### `ZIPLNnetworkfamily$clone()`

The objects of this class are cloneable with this method.

#### Usage

    ZIPLNnetworkfamily$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

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
class(fits)
#> [1] "PLNnetworkfamily" "Networkfamily"    "PLNfamily"        "R6"              
```

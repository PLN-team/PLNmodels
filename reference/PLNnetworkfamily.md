# An R6 Class to represent a collection of [`PLNnetworkfit`](https://pln-team.github.io/PLNmodels/reference/PLNnetworkfit.md)s

The function
[`PLNnetwork()`](https://pln-team.github.io/PLNmodels/reference/PLNnetwork.md)
produces an instance of this class.

This class comes with a set of methods mostly used to compare network
fits (in terms of goodness of fit) or extract one from the family (based
on penalty parameter and/or goodness of it). See the documentation for
[`getBestModel()`](https://pln-team.github.io/PLNmodels/reference/getBestModel.md),
[`getModel()`](https://pln-team.github.io/PLNmodels/reference/getModel.md)
and
[plot()](https://pln-team.github.io/PLNmodels/reference/plot.Networkfamily.md)
for the user-facing ones.

## See also

The function
[`PLNnetwork()`](https://pln-team.github.io/PLNmodels/reference/PLNnetwork.md),
the class
[`PLNnetworkfit`](https://pln-team.github.io/PLNmodels/reference/PLNnetworkfit.md)

## Super classes

[`PLNfamily`](https://pln-team.github.io/PLNmodels/reference/PLNfamily.md)
-\>
[`Networkfamily`](https://pln-team.github.io/PLNmodels/reference/Networkfamily.md)
-\> `PLNnetworkfamily`

## Methods

### Public methods

- [`PLNnetworkfamily$new()`](#method-PLNnetworkfamily-initialize)

- [`PLNnetworkfamily$stability_selection()`](#method-PLNnetworkfamily-stability_selection)

- [`PLNnetworkfamily$clone()`](#method-PLNnetworkfamily-clone)

Inherited methods

- [`PLNfamily$getModel()`](https://pln-team.github.io/PLNmodels/reference/PLNfamily.html#method-getModel)
- [`PLNfamily$print()`](https://pln-team.github.io/PLNmodels/reference/PLNfamily.html#method-print)
- [`Networkfamily$coefficient_path()`](https://pln-team.github.io/PLNmodels/reference/Networkfamily.html#method-coefficient_path)
- [`Networkfamily$getBestModel()`](https://pln-team.github.io/PLNmodels/reference/Networkfamily.html#method-getBestModel)
- [`Networkfamily$optimize()`](https://pln-team.github.io/PLNmodels/reference/Networkfamily.html#method-optimize)
- [`Networkfamily$plot()`](https://pln-team.github.io/PLNmodels/reference/Networkfamily.html#method-plot)
- [`Networkfamily$plot_objective()`](https://pln-team.github.io/PLNmodels/reference/Networkfamily.html#method-plot_objective)
- [`Networkfamily$plot_stars()`](https://pln-team.github.io/PLNmodels/reference/Networkfamily.html#method-plot_stars)
- [`Networkfamily$postTreatment()`](https://pln-team.github.io/PLNmodels/reference/Networkfamily.html#method-postTreatment)
- [`Networkfamily$show()`](https://pln-team.github.io/PLNmodels/reference/Networkfamily.html#method-show)

------------------------------------------------------------------------

### `PLNnetworkfamily$new()`

Initialize all models in the collection

#### Usage

    PLNnetworkfamily$new(penalties, data, control)

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

### `PLNnetworkfamily$stability_selection()`

Compute the stability path by stability selection

#### Usage

    PLNnetworkfamily$stability_selection(
      subsamples = NULL,
      control = PLNnetwork_param()
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
  [`PLNnetwork()`](https://pln-team.github.io/PLNmodels/reference/PLNnetwork.md)
  and
  [`PLN_param()`](https://pln-team.github.io/PLNmodels/reference/PLN_param.md)
  for details.

------------------------------------------------------------------------

### `PLNnetworkfamily$clone()`

The objects of this class are cloneable with this method.

#### Usage

    PLNnetworkfamily$clone(deep = FALSE)

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
#>  sparsifying penalty = 2.847172  sparsifying penalty = 2.62985   sparsifying penalty = 2.429116  sparsifying penalty = 2.243704  sparsifying penalty = 2.072444  sparsifying penalty = 1.914256  sparsifying penalty = 1.768142  sparsifying penalty = 1.633182  sparsifying penalty = 1.508522  sparsifying penalty = 1.393378  sparsifying penalty = 1.287023  sparsifying penalty = 1.188785  sparsifying penalty = 1.098046  sparsifying penalty = 1.014233  sparsifying penalty = 0.9368178     sparsifying penalty = 0.8653113     sparsifying penalty = 0.7992629     sparsifying penalty = 0.7382558     sparsifying penalty = 0.6819054     sparsifying penalty = 0.6298561     sparsifying penalty = 0.5817797     sparsifying penalty = 0.537373  sparsifying penalty = 0.4963558     sparsifying penalty = 0.4584694     sparsifying penalty = 0.4234748     sparsifying penalty = 0.3911513     sparsifying penalty = 0.3612951     sparsifying penalty = 0.3337177     sparsifying penalty = 0.3082453     sparsifying penalty = 0.2847172 
#>  Post-treatments
#>  DONE!
class(fits)
#> [1] "PLNnetworkfamily" "Networkfamily"    "PLNfamily"        "R6"              
```

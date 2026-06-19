# An R6 Class to represent a PLNfit in a sparse inverse covariance framework

The function
[`PLNnetwork()`](https://pln-team.github.io/PLNmodels/reference/PLNnetwork.md)
produces a collection of models which are instances of object with class
`PLNnetworkfit`. This class comes with a set of methods, some of them
being useful for the user: See the documentation for
[`plot()`](https://pln-team.github.io/PLNmodels/reference/plot.PLNnetworkfit.md)
and methods inherited from
[`PLNfit`](https://pln-team.github.io/PLNmodels/reference/PLNfit.md).

## See also

The function
[`PLNnetwork()`](https://pln-team.github.io/PLNmodels/reference/PLNnetwork.md),
the class
[`PLNnetworkfamily`](https://pln-team.github.io/PLNmodels/reference/PLNnetworkfamily.md)

## Super classes

[`PLNfit`](https://pln-team.github.io/PLNmodels/reference/PLNfit.md) -\>
[`PLNfit_fixedcov`](https://pln-team.github.io/PLNmodels/reference/PLNfit_fixedcov.md)
-\> `PLNnetworkfit`

## Active bindings

- `vcov_model`:

  character: the model used for the residual covariance

- `penalty`:

  the global level of sparsity in the current model

- `penalty_weights`:

  a matrix of weights controlling the amount of penalty element-wise.

- `n_edges`:

  number of edges if the network (non null coefficient of the sparse
  precision matrix)

- `nb_param`:

  number of parameters in the current PLN model

- `pen_loglik`:

  variational lower bound of the l1-penalized loglikelihood

- `EBIC`:

  variational lower bound of the EBIC

- `density`:

  proportion of non-null edges in the network

- `criteria`:

  a vector with loglik, penalized loglik, BIC, EBIC, ICL, R_squared,
  number of parameters, number of edges and graph density

## Methods

### Public methods

- [`PLNnetworkfit$new()`](#method-PLNnetworkfit-initialize)

- [`PLNnetworkfit$optimize()`](#method-PLNnetworkfit-optimize)

- [`PLNnetworkfit$latent_network()`](#method-PLNnetworkfit-latent_network)

- [`PLNnetworkfit$plot_network()`](#method-PLNnetworkfit-plot_network)

- [`PLNnetworkfit$show()`](#method-PLNnetworkfit-show)

- [`PLNnetworkfit$clone()`](#method-PLNnetworkfit-clone)

Inherited methods

- [`PLNfit$optimize_vestep()`](https://pln-team.github.io/PLNmodels/reference/PLNfit.html#method-optimize_vestep)
- [`PLNfit$postTreatment()`](https://pln-team.github.io/PLNmodels/reference/PLNfit.html#method-postTreatment)
- [`PLNfit$predict()`](https://pln-team.github.io/PLNmodels/reference/PLNfit.html#method-predict)
- [`PLNfit$predict_cond()`](https://pln-team.github.io/PLNmodels/reference/PLNfit.html#method-predict_cond)
- [`PLNfit$print()`](https://pln-team.github.io/PLNmodels/reference/PLNfit.html#method-print)
- [`PLNfit$update()`](https://pln-team.github.io/PLNmodels/reference/PLNfit.html#method-update)

------------------------------------------------------------------------

### `PLNnetworkfit$new()`

Initialize a `PLNnetworkfit` object

#### Usage

    PLNnetworkfit$new(data, control)

#### Arguments

- `data`:

  a named list used internally to carry the data matrices

- `control`:

  a list for controlling the optimization.

------------------------------------------------------------------------

### `PLNnetworkfit$optimize()`

Call to the C++ optimizer and update of the relevant fields

#### Usage

    PLNnetworkfit$optimize(data, config)

#### Arguments

- `data`:

  a named list used internally to carry the data matrices

- `config`:

  a list for controlling the optimization

------------------------------------------------------------------------

### `PLNnetworkfit$latent_network()`

Extract interaction network in the latent space

#### Usage

    PLNnetworkfit$latent_network(type = c("partial_cor", "support", "precision"))

#### Arguments

- `type`:

  edge value in the network. Can be "support" (binary edges),
  "precision" (coefficient of the precision matrix) or "partial_cor"
  (partial correlation between species)

#### Returns

a square matrix of size `PLNnetworkfit$n`

------------------------------------------------------------------------

### `PLNnetworkfit$plot_network()`

plot the latent network.

#### Usage

    PLNnetworkfit$plot_network(
      type = c("partial_cor", "support"),
      output = c("igraph", "corrplot"),
      edge.color = c("#F8766D", "#00BFC4"),
      remove.isolated = FALSE,
      node.labels = NULL,
      layout = layout_in_circle,
      plot = TRUE
    )

#### Arguments

- `type`:

  edge value in the network. Either "precision" (coefficient of the
  precision matrix) or "partial_cor" (partial correlation between
  species).

- `output`:

  Output type. Either `igraph` (for the network) or `corrplot` (for the
  adjacency matrix)

- `edge.color`:

  Length 2 color vector. Color for positive/negative edges. Default is
  `c("#F8766D", "#00BFC4")`. Only relevant for igraph output.

- `remove.isolated`:

  if `TRUE`, isolated node are remove before plotting. Only relevant for
  igraph output.

- `node.labels`:

  vector of character. The labels of the nodes. The default will use the
  column names ot the response matrix.

- `layout`:

  an optional igraph layout. Only relevant for igraph output.

- `plot`:

  logical. Should the final network be displayed or only sent back to
  the user. Default is `TRUE`.

------------------------------------------------------------------------

### `PLNnetworkfit$show()`

User friendly print method

#### Usage

    PLNnetworkfit$show()

------------------------------------------------------------------------

### `PLNnetworkfit$clone()`

The objects of this class are cloneable with this method.

#### Usage

    PLNnetworkfit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
if (FALSE) { # \dontrun{
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
nets <- PLNnetwork(Abundance ~ 1, data = trichoptera)
myPLNnet <- getBestModel(nets)
class(myPLNnet)
print(myPLNnet)
} # }
```

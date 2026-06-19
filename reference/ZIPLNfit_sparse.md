# An R6 Class to represent a ZIPLNfit in a standard, general framework, with sparse inverse residual covariance

An R6 Class to represent a ZIPLNfit in a standard, general framework,
with sparse inverse residual covariance

## Super class

[`ZIPLNfit`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit.md)
-\> `ZIPLNfit_sparse`

## Active bindings

- `penalty`:

  the global level of sparsity in the current model

- `penalty_weights`:

  a matrix of weights controlling the amount of penalty element-wise.

- `n_edges`:

  number of edges if the network (non null coefficient of the sparse
  precision matrix)

- `nb_param_pln`:

  number of parameters in the PLN part of the current model

- `vcov_model`:

  character: the model used for the residual covariance

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

- [`ZIPLNfit_sparse$new()`](#method-ZIPLNfit_sparse-initialize)

- [`ZIPLNfit_sparse$latent_network()`](#method-ZIPLNfit_sparse-latent_network)

- [`ZIPLNfit_sparse$plot_network()`](#method-ZIPLNfit_sparse-plot_network)

- [`ZIPLNfit_sparse$clone()`](#method-ZIPLNfit_sparse-clone)

Inherited methods

- [`ZIPLNfit$optimize()`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit.html#method-optimize)
- [`ZIPLNfit$optimize_vestep()`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit.html#method-optimize_vestep)
- [`ZIPLNfit$predict()`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit.html#method-predict)
- [`ZIPLNfit$print()`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit.html#method-print)
- [`ZIPLNfit$show()`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit.html#method-show)
- [`ZIPLNfit$update()`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit.html#method-update)

------------------------------------------------------------------------

### `ZIPLNfit_sparse$new()`

Initialize a
[`ZIPLNfit_fixed`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit_fixed.md)
model

#### Usage

    ZIPLNfit_sparse$new(data, control)

#### Arguments

- `data`:

  a named list used internally to carry the data matrices

- `control`:

  a list for controlling the optimization. See details.

------------------------------------------------------------------------

### `ZIPLNfit_sparse$latent_network()`

Extract interaction network in the latent space

#### Usage

    ZIPLNfit_sparse$latent_network(type = c("partial_cor", "support", "precision"))

#### Arguments

- `type`:

  edge value in the network. Can be "support" (binary edges),
  "precision" (coefficient of the precision matrix) or "partial_cor"
  (partial correlation between species)

#### Returns

a square matrix of size `ZIPLNfit_sparse$n`

------------------------------------------------------------------------

### `ZIPLNfit_sparse$plot_network()`

plot the latent network.

#### Usage

    ZIPLNfit_sparse$plot_network(
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

### `ZIPLNfit_sparse$clone()`

The objects of this class are cloneable with this method.

#### Usage

    ZIPLNfit_sparse$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
if (FALSE) { # \dontrun{
# See other examples in function ZIPLN
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
myPLN <- ZIPLN(Abundance ~ 1, data = trichoptera, control=  ZIPLN_param(penalty = 1))
class(myPLN)
print(myPLN)
plot(myPLN)
} # }
```

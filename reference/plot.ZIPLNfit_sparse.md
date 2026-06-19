# Extract and plot the network (partial correlation, support or inverse covariance) from a [`ZIPLNfit_sparse`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit_sparse.md) object

Extract and plot the network (partial correlation, support or inverse
covariance) from a
[`ZIPLNfit_sparse`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit_sparse.md)
object

## Usage

``` r
# S3 method for class 'ZIPLNfit_sparse'
plot(
  x,
  type = c("partial_cor", "support"),
  output = c("igraph", "corrplot"),
  edge.color = c("#F8766D", "#00BFC4"),
  remove.isolated = FALSE,
  node.labels = NULL,
  layout = layout_in_circle,
  plot = TRUE,
  ...
)
```

## Arguments

- x:

  an R6 object with class
  [`ZIPLNfit_sparse`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit_sparse.md)

- type:

  character. Value of the weight of the edges in the network, either
  "partial_cor" (partial correlation) or "support" (binary). Default is
  `"partial_cor"`.

- output:

  the type of output used: either 'igraph' or 'corrplot'. Default is
  `'igraph'`.

- edge.color:

  Length 2 color vector. Color for positive/negative edges. Default is
  `c("#F8766D", "#00BFC4")`. Only relevant for igraph output.

- remove.isolated:

  if `TRUE`, isolated node are remove before plotting. Only relevant for
  igraph output.

- node.labels:

  vector of character. The labels of the nodes. The default will use the
  column names ot the response matrix.

- layout:

  an optional igraph layout. Only relevant for igraph output.

- plot:

  logical. Should the final network be displayed or only sent back to
  the user. Default is `TRUE`.

- ...:

  Not used (S3 compatibility).

## Value

Send back an invisible object (igraph or Matrix, depending on the output
chosen) and optionally displays a graph (via igraph or corrplot for
large ones)

## Examples

``` r
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
fit <- ZIPLN(Abundance ~ 1, data = trichoptera, control = ZIPLN_param(penalty = 0.1))
#> 
#>  Initialization...
#>  Adjusting a ZI-PLN model with sparse covariance model and single specific parameter(s) in Zero inflation component.
#>  DONE!
if (FALSE) { # \dontrun{
plot(fit)
} # }
```

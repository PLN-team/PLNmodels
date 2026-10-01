# Extract and plot the network (partial correlation, support or inverse covariance) from a [`PLNnetworkfit`](https://pln-team.github.io/PLNmodels/reference/PLNnetworkfit.md) object

Extract and plot the network (partial correlation, support or inverse
covariance) from a
[`PLNnetworkfit`](https://pln-team.github.io/PLNmodels/reference/PLNnetworkfit.md)
object

## Usage

``` r
# S3 method for class 'PLNnetworkfit'
plot(
  x,
  type = c("partial_cor", "support"),
  output = c("igraph", "corrplot"),
  edge.color = c("#F8766D", "#00BFC4"),
  remove.isolated = FALSE,
  node.labels = NULL,
  layout = layout_in_circle,
  edge.alpha = 0.2,
  plot = TRUE,
  ...
)
```

## Arguments

- x:

  an R6 object with class
  [`PLNnetworkfit`](https://pln-team.github.io/PLNmodels/reference/PLNnetworkfit.md)

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

- edge.alpha:

  opacity of the weakest edge, the strongest one being fully opaque, so
  that the strength of an edge can be read off a dense network. Default
  is `0.2`. Set it to `1` for uniformly opaque edges. Only relevant for
  igraph output with `type = "partial_cor"`.

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
fits <- PLNnetwork(Abundance ~ 1, data = trichoptera)
#> 
#>  Initialization...
#>  Adjusting 30 PLN with sparse inverse covariance estimation
#>  Joint optimization alternating gradient descent and graphical-lasso
#>  sparsifying penalty = 2.847172  sparsifying penalty = 2.62985   sparsifying penalty = 2.429116  sparsifying penalty = 2.243704  sparsifying penalty = 2.072444  sparsifying penalty = 1.914256  sparsifying penalty = 1.768142  sparsifying penalty = 1.633182  sparsifying penalty = 1.508522  sparsifying penalty = 1.393378  sparsifying penalty = 1.287023  sparsifying penalty = 1.188785  sparsifying penalty = 1.098046  sparsifying penalty = 1.014233  sparsifying penalty = 0.9368178     sparsifying penalty = 0.8653113     sparsifying penalty = 0.7992629     sparsifying penalty = 0.7382558     sparsifying penalty = 0.6819054     sparsifying penalty = 0.6298561     sparsifying penalty = 0.5817797     sparsifying penalty = 0.537373  sparsifying penalty = 0.4963558     sparsifying penalty = 0.4584694     sparsifying penalty = 0.4234748     sparsifying penalty = 0.3911513     sparsifying penalty = 0.3612951     sparsifying penalty = 0.3337177     sparsifying penalty = 0.3082453     sparsifying penalty = 0.2847172 
#>  Post-treatments
#>  DONE!
myNet <- getBestModel(fits)
if (FALSE) { # \dontrun{
plot(myNet)
} # }
```

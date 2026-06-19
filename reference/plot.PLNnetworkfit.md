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
#>  sparsifying penalty = 4.479926  sparsifying penalty = 4.137977  sparsifying penalty = 3.822129  sparsifying penalty = 3.530389  sparsifying penalty = 3.260917  sparsifying penalty = 3.012014  sparsifying penalty = 2.78211   sparsifying penalty = 2.569754  sparsifying penalty = 2.373607  sparsifying penalty = 2.192431  sparsifying penalty = 2.025085  sparsifying penalty = 1.870512  sparsifying penalty = 1.727737  sparsifying penalty = 1.595861  sparsifying penalty = 1.47405   sparsifying penalty = 1.361537  sparsifying penalty = 1.257612  sparsifying penalty = 1.16162   sparsifying penalty = 1.072954  sparsifying penalty = 0.9910565     sparsifying penalty = 0.91541   sparsifying penalty = 0.8455375     sparsifying penalty = 0.7809984     sparsifying penalty = 0.7213854     sparsifying penalty = 0.6663227     sparsifying penalty = 0.6154629     sparsifying penalty = 0.5684851     sparsifying penalty = 0.5250931     sparsifying penalty = 0.4850132     sparsifying penalty = 0.4479926 
#>  Post-treatments
#>  DONE!
myNet <- getBestModel(fits)
if (FALSE) { # \dontrun{
plot(myNet)
} # }
```

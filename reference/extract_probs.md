# Extract edge selection frequency in bootstrap subsamples

Extracts edge selection frequency in networks reconstructed from
bootstrap subsamples during the stars stability selection procedure, as
either a matrix or a named vector. In the latter case, edge names follow
igraph naming convention.

## Usage

``` r
extract_probs(
  Robject,
  penalty = NULL,
  index = NULL,
  crit = c("StARS", "BIC", "EBIC"),
  format = c("matrix", "vector"),
  tol = 1e-05
)
```

## Arguments

- Robject:

  an object with class
  [`PLNnetworkfamily`](https://pln-team.github.io/PLNmodels/reference/PLNnetworkfamily.md),
  i.e. an output from
  [`PLNnetwork()`](https://pln-team.github.io/PLNmodels/reference/PLNnetwork.md)

- penalty:

  penalty used for the bootstrap subsamples

- index:

  Integer index of the model to be returned. Only the first value is
  taken into account.

- crit:

  a character for the criterion used to performed the selection. Either
  "BIC", "ICL", "EBIC", "StARS", "R_squared". Default is `ICL` for
  `PLNPCA`, and `BIC` for `PLNnetwork`. If StARS (Stability Approach to
  Regularization Selection) is chosen and stability selection was not
  yet performed, the function will call the method
  [`stability_selection()`](https://pln-team.github.io/PLNmodels/reference/stability_selection.md)
  with default argument.

- format:

  output format. Either a matrix (default) or a named vector.

- tol:

  tolerance for rounding error when comparing penalties.

## Value

Either a matrix or named vector of edge-wise probabilities. In the
latter case, edge names follow igraph convention.

## Examples

``` r
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
nets <- PLNnetwork(Abundance ~ 1 + offset(log(Offset)), data = trichoptera)
#> 
#>  Initialization...
#>  Adjusting 30 PLN with sparse inverse covariance estimation
#>  Joint optimization alternating gradient descent and graphical-lasso
#>  sparsifying penalty = 2.123258  sparsifying penalty = 1.961192  sparsifying penalty = 1.811496  sparsifying penalty = 1.673226  sparsifying penalty = 1.54551   sparsifying penalty = 1.427542  sparsifying penalty = 1.318579  sparsifying penalty = 1.217933  sparsifying penalty = 1.124969  sparsifying penalty = 1.039101  sparsifying penalty = 0.9597877     sparsifying penalty = 0.886528  sparsifying penalty = 0.81886   sparsifying penalty = 0.7563572     sparsifying penalty = 0.6986251     sparsifying penalty = 0.6452996     sparsifying penalty = 0.5960444     sparsifying penalty = 0.5505489     sparsifying penalty = 0.508526  sparsifying penalty = 0.4697106     sparsifying penalty = 0.433858  sparsifying penalty = 0.400742  sparsifying penalty = 0.3701537     sparsifying penalty = 0.3419002     sparsifying penalty = 0.3158032     sparsifying penalty = 0.2916982     sparsifying penalty = 0.2694332     sparsifying penalty = 0.2488676     sparsifying penalty = 0.2298717     sparsifying penalty = 0.2123258 
#>  Post-treatments
#>  DONE!
if (FALSE) { # \dontrun{
stability_selection(nets)
probs <- extract_probs(nets, crit = "StARS", format = "vector")
probs
} # }

if (FALSE) { # \dontrun{
## Add edge attributes to graph using igraph
net_stars <- getBestModel(nets, "StARS")
g <- plot(net_stars, type = "partial_cor", plot=F)
library(igraph)
E(g)$prob <- probs[as_ids(E(g))]
g
} # }
```

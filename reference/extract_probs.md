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
#>  sparsifying penalty = 0.8883632     sparsifying penalty = 0.8205552     sparsifying penalty = 0.7579229     sparsifying penalty = 0.7000713     sparsifying penalty = 0.6466354     sparsifying penalty = 0.5972783     sparsifying penalty = 0.5516886     sparsifying penalty = 0.5095787     sparsifying penalty = 0.470683  sparsifying penalty = 0.4347561     sparsifying penalty = 0.4015716     sparsifying penalty = 0.37092   sparsifying penalty = 0.342608  sparsifying penalty = 0.316457  sparsifying penalty = 0.2923021     sparsifying penalty = 0.2699909     sparsifying penalty = 0.2493827     sparsifying penalty = 0.2303476     sparsifying penalty = 0.2127653     sparsifying penalty = 0.1965251     sparsifying penalty = 0.1815246     sparsifying penalty = 0.1676689     sparsifying penalty = 0.1548709     sparsifying penalty = 0.1430497     sparsifying penalty = 0.1321309     sparsifying penalty = 0.1220454     sparsifying penalty = 0.1127298     sparsifying penalty = 0.1041253     sparsifying penalty = 0.09617746    sparsifying penalty = 0.08883632 
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

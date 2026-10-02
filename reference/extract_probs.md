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
#>  sparsifying penalty = 0.3985143     sparsifying penalty = 0.368096  sparsifying penalty = 0.3399996     sparsifying penalty = 0.3140477     sparsifying penalty = 0.2900767     sparsifying penalty = 0.2679354     sparsifying penalty = 0.2474841     sparsifying penalty = 0.2285939     sparsifying penalty = 0.2111455     sparsifying penalty = 0.1950289     sparsifying penalty = 0.1801425     sparsifying penalty = 0.1663924     sparsifying penalty = 0.1536918     sparsifying penalty = 0.1419607     sparsifying penalty = 0.1311249     sparsifying penalty = 0.1211163     sparsifying penalty = 0.1118716     sparsifying penalty = 0.1033325     sparsifying penalty = 0.09544523    sparsifying penalty = 0.08815998    sparsifying penalty = 0.0814308     sparsifying penalty = 0.07521526    sparsifying penalty = 0.06947414    sparsifying penalty = 0.06417124    sparsifying penalty = 0.0592731     sparsifying penalty = 0.05474884    sparsifying penalty = 0.05056991    sparsifying penalty = 0.04670995    sparsifying penalty = 0.04314462    sparsifying penalty = 0.03985143 
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

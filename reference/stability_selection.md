# Compute the stability path by stability selection

This function computes the StARS stability criteria over a path of
penalties. If a path has already been computed, the functions stops with
a message unless `force = TRUE` has been specified.

## Usage

``` r
stability_selection(
  Robject,
  subsamples = NULL,
  control = PLNnetwork_param(),
  force = FALSE
)
```

## Arguments

- Robject:

  an object with class
  [`PLNnetworkfamily`](https://pln-team.github.io/PLNmodels/reference/PLNnetworkfamily.md)
  or
  [`ZIPLNnetworkfamily`](https://pln-team.github.io/PLNmodels/reference/ZIPLNnetworkfamily.md),
  i.e. an output from
  [`PLNnetwork()`](https://pln-team.github.io/PLNmodels/reference/PLNnetwork.md)
  or
  [`ZIPLNnetwork()`](https://pln-team.github.io/PLNmodels/reference/ZIPLNnetwork.md)

- subsamples:

  a list of vectors describing the subsamples. The number of vectors (or
  list length) determines th number of subsamples used in the stability
  selection. Automatically set to 20 subsamples with size `10*sqrt(n)`
  if `n >= 144` and `0.8*n` otherwise following Liu et al. (2010)
  recommendations.

- control:

  a list controlling the main optimization process in each call to
  [`PLNnetwork()`](https://pln-team.github.io/PLNmodels/reference/PLNnetwork.md)
  or
  [`ZIPLNnetwork()`](https://pln-team.github.io/PLNmodels/reference/ZIPLNnetwork.md).
  See
  [`PLN_param()`](https://pln-team.github.io/PLNmodels/reference/PLN_param.md)
  or
  [`ZIPLN_param()`](https://pln-team.github.io/PLNmodels/reference/ZIPLN_param.md)
  for details.

- force:

  force computation of the stability path, even if a previous one has
  been detected.

## Value

the list of subsamples. The estimated probabilities of selection of the
edges are stored in the fields `stability_path` of the initial Robject
with class
[`Networkfamily`](https://pln-team.github.io/PLNmodels/reference/Networkfamily.md)

## Examples

``` r
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
fits <- PLNnetwork(Abundance ~ 1, data = trichoptera)
#> 
#>  Initialization...
#>  Adjusting 30 PLN with sparse inverse covariance estimation
#>  Joint optimization alternating gradient descent and graphical-lasso
#>  sparsifying penalty = 0.5246191     sparsifying penalty = 0.4845754     sparsifying penalty = 0.4475881     sparsifying penalty = 0.4134241     sparsifying penalty = 0.3818678     sparsifying penalty = 0.3527202     sparsifying penalty = 0.3257973     sparsifying penalty = 0.3009295     sparsifying penalty = 0.2779598     sparsifying penalty = 0.2567434     sparsifying penalty = 0.2371464     sparsifying penalty = 0.2190452     sparsifying penalty = 0.2023257     sparsifying penalty = 0.1868823     sparsifying penalty = 0.1726178     sparsifying penalty = 0.159442  sparsifying penalty = 0.1472719     sparsifying penalty = 0.1360308     sparsifying penalty = 0.1256477     sparsifying penalty = 0.1160571     sparsifying penalty = 0.1071986     sparsifying penalty = 0.09901618    sparsifying penalty = 0.09145836    sparsifying penalty = 0.08447742    sparsifying penalty = 0.07802933    sparsifying penalty = 0.07207342    sparsifying penalty = 0.06657212    sparsifying penalty = 0.06149072    sparsifying penalty = 0.05679719    sparsifying penalty = 0.05246191 
#>  Post-treatments
#>  DONE!
if (FALSE) { # \dontrun{
n <- nrow(trichoptera)
subs <- replicate(10, sample.int(n, size = n/2), simplify = FALSE)
stability_selection(nets, subsamples = subs)
} # }
```

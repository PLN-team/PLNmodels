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
#>  sparsifying penalty = 4.479926  sparsifying penalty = 4.137977  sparsifying penalty = 3.822129  sparsifying penalty = 3.530389  sparsifying penalty = 3.260917  sparsifying penalty = 3.012014  sparsifying penalty = 2.78211   sparsifying penalty = 2.569754  sparsifying penalty = 2.373607  sparsifying penalty = 2.192431  sparsifying penalty = 2.025085  sparsifying penalty = 1.870512  sparsifying penalty = 1.727737  sparsifying penalty = 1.595861  sparsifying penalty = 1.47405   sparsifying penalty = 1.361537  sparsifying penalty = 1.257612  sparsifying penalty = 1.16162   sparsifying penalty = 1.072954  sparsifying penalty = 0.9910565     sparsifying penalty = 0.91541   sparsifying penalty = 0.8455375     sparsifying penalty = 0.7809984     sparsifying penalty = 0.7213854     sparsifying penalty = 0.6663227     sparsifying penalty = 0.6154629     sparsifying penalty = 0.5684851     sparsifying penalty = 0.5250931     sparsifying penalty = 0.4850132     sparsifying penalty = 0.4479926 
#>  Post-treatments
#>  DONE!
if (FALSE) { # \dontrun{
n <- nrow(trichoptera)
subs <- replicate(10, sample.int(n, size = n/2), simplify = FALSE)
stability_selection(nets, subsamples = subs)
} # }
```

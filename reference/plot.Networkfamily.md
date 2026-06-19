# Display various outputs (goodness-of-fit criteria, robustness, diagnostic) associated with a collection of network fits (either [`PLNnetworkfamily`](https://pln-team.github.io/PLNmodels/reference/PLNnetworkfamily.md) or [`ZIPLNnetworkfamily`](https://pln-team.github.io/PLNmodels/reference/ZIPLNnetworkfamily.md))

Display various outputs (goodness-of-fit criteria, robustness,
diagnostic) associated with a collection of network fits (either
[`PLNnetworkfamily`](https://pln-team.github.io/PLNmodels/reference/PLNnetworkfamily.md)
or
[`ZIPLNnetworkfamily`](https://pln-team.github.io/PLNmodels/reference/ZIPLNnetworkfamily.md))

## Usage

``` r
# S3 method for class 'Networkfamily'
plot(
  x,
  type = c("criteria", "stability", "diagnostic"),
  criteria = c("loglik", "pen_loglik", "BIC", "EBIC"),
  reverse = FALSE,
  log.x = TRUE,
  stability = 0.9,
  ...
)

# S3 method for class 'PLNnetworkfamily'
plot(
  x,
  type = c("criteria", "stability", "diagnostic"),
  criteria = c("loglik", "pen_loglik", "BIC", "EBIC"),
  reverse = FALSE,
  log.x = TRUE,
  stability = 0.9,
  ...
)

# S3 method for class 'ZIPLNnetworkfamily'
plot(
  x,
  type = c("criteria", "stability", "diagnostic"),
  criteria = c("loglik", "pen_loglik", "BIC", "EBIC"),
  reverse = FALSE,
  log.x = TRUE,
  stability = 0.9,
  ...
)
```

## Arguments

- x:

  an R6 object with class
  [`PLNnetworkfamily`](https://pln-team.github.io/PLNmodels/reference/PLNnetworkfamily.md)
  or
  [`ZIPLNnetworkfamily`](https://pln-team.github.io/PLNmodels/reference/ZIPLNnetworkfamily.md)

- type:

  a character, either "criteria", "stability" or "diagnostic" for the
  type of plot.

- criteria:

  Vector of criteria to plot, to be selected among "loglik"
  (log-likelihood), "BIC", "ICL", "R_squared", "EBIC" and "pen_loglik"
  (penalized log-likelihood). Default is c("loglik", "pen_loglik",
  "BIC", "EBIC"). Only used when `type = "criteria"`.

- reverse:

  A logical indicating whether to plot the value of the criteria in the
  "natural" direction (loglik - 0.5 penalty) or in the "reverse"
  direction (-2 loglik + penalty). Default to FALSE, i.e use the natural
  direction, on the same scale as the log-likelihood.

- log.x:

  logical: should the x-axis be represented in log-scale? Default is
  `TRUE`.

- stability:

  scalar: the targeted level of stability in stability plot. Default is
  .9.

- ...:

  additional parameters for S3 compatibility. Not used

## Value

Produces either a diagnostic plot (with `type = 'diagnostic'`), a
stability plot (with `type = 'stability'`) or the evolution of the
criteria of the different models considered (with `type = 'criteria'`,
the default).

## Details

The BIC and ICL criteria have the form 'loglik - 1/2 \* penalty' so that
they are on the same scale as the model log-likelihood. You can change
this direction and use the alternate form '-2\*loglik + penalty', as
some authors do, by setting `reverse = TRUE`.

## Functions

- `plot(PLNnetworkfamily)`: Display various outputs associated with a
  collection of network fits

- `plot(ZIPLNnetworkfamily)`: Display various outputs associated with a
  collection of network fits

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
plot(fits)
} # }
```

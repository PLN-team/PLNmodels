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
#>  sparsifying penalty = 2.847172  sparsifying penalty = 2.62985   sparsifying penalty = 2.429116  sparsifying penalty = 2.243704  sparsifying penalty = 2.072444  sparsifying penalty = 1.914256  sparsifying penalty = 1.768142  sparsifying penalty = 1.633182  sparsifying penalty = 1.508522  sparsifying penalty = 1.393378  sparsifying penalty = 1.287023  sparsifying penalty = 1.188785  sparsifying penalty = 1.098046  sparsifying penalty = 1.014233  sparsifying penalty = 0.9368178     sparsifying penalty = 0.8653113     sparsifying penalty = 0.7992629     sparsifying penalty = 0.7382558     sparsifying penalty = 0.6819054     sparsifying penalty = 0.6298561     sparsifying penalty = 0.5817797     sparsifying penalty = 0.537373  sparsifying penalty = 0.4963558     sparsifying penalty = 0.4584694     sparsifying penalty = 0.4234748     sparsifying penalty = 0.3911513     sparsifying penalty = 0.3612951     sparsifying penalty = 0.3337177     sparsifying penalty = 0.3082453     sparsifying penalty = 0.2847172 
#>  Post-treatments
#>  DONE!
if (FALSE) { # \dontrun{
plot(fits)
} # }
```

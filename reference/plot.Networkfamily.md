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
#>  sparsifying penalty = 0.5246191     sparsifying penalty = 0.4845754     sparsifying penalty = 0.4475881     sparsifying penalty = 0.4134241     sparsifying penalty = 0.3818678     sparsifying penalty = 0.3527202     sparsifying penalty = 0.3257973     sparsifying penalty = 0.3009295     sparsifying penalty = 0.2779598     sparsifying penalty = 0.2567434     sparsifying penalty = 0.2371464     sparsifying penalty = 0.2190452     sparsifying penalty = 0.2023257     sparsifying penalty = 0.1868823     sparsifying penalty = 0.1726178     sparsifying penalty = 0.159442  sparsifying penalty = 0.1472719     sparsifying penalty = 0.1360308     sparsifying penalty = 0.1256477     sparsifying penalty = 0.1160571     sparsifying penalty = 0.1071986     sparsifying penalty = 0.09901618    sparsifying penalty = 0.09145836    sparsifying penalty = 0.08447742    sparsifying penalty = 0.07802933    sparsifying penalty = 0.07207342    sparsifying penalty = 0.06657212    sparsifying penalty = 0.06149072    sparsifying penalty = 0.05679719    sparsifying penalty = 0.05246191 
#>  Post-treatments
#>  DONE!
if (FALSE) { # \dontrun{
plot(fits)
} # }
```

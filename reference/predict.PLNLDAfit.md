# Predict group of new samples

Predict group of new samples

## Usage

``` r
# S3 method for class 'PLNLDAfit'
predict(
  object,
  newdata,
  type = c("posterior", "response", "scores"),
  scale = c("log", "prob"),
  prior = NULL,
  control = PLN_param(backend = "nlopt"),
  ...
)
```

## Arguments

- object:

  an R6 object with class
  [`PLNLDAfit`](https://pln-team.github.io/PLNmodels/reference/PLNLDAfit.md)

- newdata:

  A data frame in which to look for variables, offsets and counts with
  which to predict.

- type:

  The type of prediction required. The default are posterior
  probabilities for each group (in either unnormalized log-scale or
  natural probabilities, see "scale" for details), "response" is the
  group with maximal posterior probability and "scores" is the average
  score along each separation axis in the latent space, with weights
  equal to the posterior probabilities.

- scale:

  The scale used for the posterior probability. Either log-scale ("log",
  default) or natural probabilities summing up to 1 ("prob").

- prior:

  User-specified prior group probabilities in the new data. If NULL
  (default), prior probabilities are computed from the learning set.

- control:

  a list for controlling the optimization. See
  [`PLN()`](https://pln-team.github.io/PLNmodels/reference/PLN.md) for
  details.

- ...:

  additional parameters for S3 compatibility. Not used

## Value

A matrix of posterior probabilities for each group (if type =
"posterior"), a matrix of (average) scores in the latent space (if type
= "scores") or a vector of predicted groups (if type = "response").

## Examples

``` r
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
myLDA <- PLNLDA(Abundance ~ 0 + offset(log(Offset)),
                grouping = Group,
                data = trichoptera)
#> 
#>  Performing discriminant Analysis...
#>  DONE!
if (FALSE) { # \dontrun{
post_probs <- predict(myLDA, newdata = trichoptera, type = "posterior", scale = "prob")
head(round(post_probs, digits = 3))
predicted_group <- predict(myLDA, newdata = trichoptera, type = "response")
table(predicted_group, trichoptera$Group, dnn = c("predicted", "true"))
} # }
```

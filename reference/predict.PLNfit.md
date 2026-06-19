# Predict counts of a new sample

Predict counts of a new sample

## Usage

``` r
# S3 method for class 'PLNfit'
predict(
  object,
  newdata,
  responses = NULL,
  level = 1,
  type = c("link", "response"),
  ...
)
```

## Arguments

- object:

  an R6 object with class
  [`PLNfit`](https://pln-team.github.io/PLNmodels/reference/PLNfit.md)

- newdata:

  A data frame in which to look for variables and offsets with which to
  predict

- responses:

  Optional data frame containing the count of the observed variables
  (matching the names of the provided as data in the PLN function),
  assuming the interest in in testing the model.

- level:

  Optional integer value the level to be used in obtaining the
  predictions. Level zero corresponds to the population predictions
  (default if `responses` is not provided) while level one (default)
  corresponds to predictions after evaluating the variational parameters
  for the new data.

- type:

  The type of prediction required. The default is on the scale of the
  linear predictors (i.e. log average count)

- ...:

  additional parameters for S3 compatibility. Not used

## Value

A matrix of predicted log-counts (if `type = "link"`) or predicted
counts (if `type = "response"`).

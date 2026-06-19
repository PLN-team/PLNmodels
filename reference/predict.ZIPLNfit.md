# Predict counts of a new sample

Predict counts of a new sample

## Usage

``` r
# S3 method for class 'ZIPLNfit'
predict(
  object,
  newdata,
  responses = NULL,
  level = 1,
  type = c("link", "response", "deflated"),
  ...
)
```

## Arguments

- object:

  an R6 object with class
  [`ZIPLNfit`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit.md)

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

  Scale used for the prediction. Either `"link"` (default, predicted
  positions in the latent space), `"response"` (predicted average
  counts, accounting for zero-inflation) or `"deflated"` (predicted
  average counts, not accounting for zero-inflation and using only the
  PLN part of the model).

- ...:

  additional parameters for S3 compatibility. Not used

## Value

A matrix of predicted log-counts (if `type = "link"`) or predicted
counts, accounting for zero-inflation (if `type = "response"`) or not
(if `type = "deflated"`).

## Details

Note that `level = 1` can only be used if responses are provided, as the
variational parameters can't be estimated otherwise. In the absence of
responses, `level` is ignored and the fitted values are returned

Note also that when `type = "response"` corresponds to predicting values
with \\(1 - \pi)A\\, where \\A\\ is the average count in the PLN part of
the model and \\\pi\\ the probability of zero-inflation, whereas
`type = "deflated"` corresponds to \\A\\.

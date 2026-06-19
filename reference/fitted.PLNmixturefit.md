# Extracts model fitted values from objects returned by [`PLNmixture()`](https://pln-team.github.io/PLNmodels/reference/PLNmixture.md) and its variants

Extracts model fitted values from objects returned by
[`PLNmixture()`](https://pln-team.github.io/PLNmodels/reference/PLNmixture.md)
and its variants

## Usage

``` r
# S3 method for class 'PLNmixturefit'
fitted(object, ...)
```

## Arguments

- object:

  an R6 object with class
  [`PLNmixturefit`](https://pln-team.github.io/PLNmodels/reference/PLNmixturefit.md)

- ...:

  additional parameters for S3 compatibility. Not used

## Value

A matrix of Fitted values extracted from the object object.

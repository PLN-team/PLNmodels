# Extracts model fitted values from objects returned by [`ZIPLN()`](https://pln-team.github.io/PLNmodels/reference/ZIPLN.md) and its variants

Extracts model fitted values from objects returned by
[`ZIPLN()`](https://pln-team.github.io/PLNmodels/reference/ZIPLN.md) and
its variants

## Usage

``` r
# S3 method for class 'ZIPLNfit'
fitted(object, ...)
```

## Arguments

- object:

  an R6 object with class
  [`ZIPLNfit`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit.md)

- ...:

  additional parameters for S3 compatibility. Not used

## Value

A matrix of Fitted values extracted from the object object.

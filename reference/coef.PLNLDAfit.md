# Extracts model coefficients from objects returned by [`PLNLDA()`](https://pln-team.github.io/PLNmodels/reference/PLNLDA.md)

The method for objects returned by
[`PLNLDA()`](https://pln-team.github.io/PLNmodels/reference/PLNLDA.md)
only returns coefficients associated to the \$\$\Theta\$\$ part of the
model (see the PLNLDA vignette for mathematical details).

## Usage

``` r
# S3 method for class 'PLNLDAfit'
coef(object, ...)
```

## Arguments

- object:

  an R6 object with class PLNLDAfit

- ...:

  additional parameters for S3 compatibility. Not used

## Value

Either NULL or a matrix of coefficients extracted from the PLNLDAfit
model.

## Examples

``` r
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
myPLNLDA <- PLNLDA(Abundance ~ Wind, grouping = Group, data = trichoptera)
#> 
#>  Performing discriminant Analysis...
#>  DONE!
coef(myPLNLDA)
#>            Wind
#> Che -0.08824877
#> Hyc  1.02627396
#> Hym  0.06692173
#> Hys -0.35104838
#> Psy  0.33811516
#> Aga  0.24191742
#> Glo  0.13439612
#> Ath  0.05742486
#> Cea  0.60508093
#> Ced  0.13311326
#> Set  0.38657209
#> All  0.10828058
#> Han  0.16578080
#> Hfo  0.37992782
#> Hsp  0.31509351
#> Hve  0.15861340
#> Sta  0.26039977
```

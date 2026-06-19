# Mixture visualization of a [`PLNmixturefit`](https://pln-team.github.io/PLNmodels/reference/PLNmixturefit.md) object

Represent the result of the clustering either by coloring the individual
in a two-dimension PCA factor map, or by representing the expected
matrix of count reorder according to the clustering.

## Usage

``` r
# S3 method for class 'PLNmixturefit'
plot(x, type = c("pca", "matrix"), main = NULL, plot = TRUE, ...)
```

## Arguments

- x:

  an R6 object with class
  [`PLNmixturefit`](https://pln-team.github.io/PLNmodels/reference/PLNmixturefit.md)

- type:

  character for the type of plot, either "pca", for or "matrix". Default
  is `"pca"`.

- main:

  character. A title for the plot. If NULL (the default), an hopefully
  appropriate title will be used.

- plot:

  logical. Should the plot be displayed or sent back as
  [`ggplot2::ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
  object

- ...:

  Not used (S3 compatibility).

## Value

a
[`ggplot2::ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
graphic

## Examples

``` r
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
myPLN <- PLNmixture(Abundance ~ 1 + offset(log(Offset)),
           data = trichoptera, control = PLNmixture_param(smoothing = "none"))  %>% getBestModel()
#> 
#>  Initialization...
#> 
#>  Adjusting 5 PLN mixture models.
#>  number of cluster = 1   number of cluster = 2   number of cluster = 3   number of cluster = 4   number of cluster = 5 
#>  Post-treatments
#>  DONE!
if (FALSE) { # \dontrun{
plot(myPLN, "pca")
plot(myPLN, "matrix")
} # }
```

# LDA visualization (individual and/or variable factor map(s)) for a [`PLNPCAfit`](https://pln-team.github.io/PLNmodels/reference/PLNPCAfit.md) object

LDA visualization (individual and/or variable factor map(s)) for a
[`PLNPCAfit`](https://pln-team.github.io/PLNmodels/reference/PLNPCAfit.md)
object

## Usage

``` r
# S3 method for class 'PLNLDAfit'
plot(
  x,
  map = c("both", "individual", "variable"),
  nb_axes = min(3, x$rank),
  axes = seq.int(min(2, x$rank)),
  var_cols = "var_colors",
  plot = TRUE,
  main = NULL,
  ...
)
```

## Arguments

- x:

  an R6 object with class PLNPCAfit

- map:

  the type of output for the PCA visualization: either "individual",
  "variable" or "both". Default is "both".

- nb_axes:

  scalar: the number of axes to be considered when map = "both". The
  default is min(3,rank).

- axes:

  numeric, the axes to use for the plot when map = "individual" or
  "variable". Default it c(1,min(rank))

- var_cols:

  a character or factor to define the color associated with the
  variables. By default, all variables receive the default color of the
  current palette.

- plot:

  logical. Should the plot be displayed or sent back as
  [`ggplot2::ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
  object

- main:

  character. A title for the single plot (individual or variable factor
  map). If NULL (the default), an hopefully appropriate title will be
  used.

- ...:

  Not used (S3 compatibility).

## Value

displays an individual and/or variable factor maps for the corresponding
axes, and/or sends back a
[`ggplot2::ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
or gtable object

## Examples

``` r
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
myPLNLDA <- PLNLDA(Abundance ~ 1, grouping = Group, data = trichoptera)
#> 
#>  Performing discriminant Analysis...
#>  DONE!
if (FALSE) { # \dontrun{
plot(myPLNLDA, map = "individual", nb_axes = 2)
} # }
```

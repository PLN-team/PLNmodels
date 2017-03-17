# PLNmodels: Poisson lognormal models

## Description

The Poisson lognormal model and variants can be used for a variety of multivariate problems when count data are at play (including PCA for count data and network inference). This package implements an efficient algorithm to fit such models accompanied with a set of functions for vizualisation and diagnostic.    
## Installation

```
devtools::install_github("jchiquet/PLNmodels", build_vignettes=TRUE)
```

## Use and example

See the package [vignette](https://github.com/jchiquet/PLNmodels/blob/master/vignettes/trichoptera.Rmd) running PLNPCA on the *ade4* Trichoptera data set.

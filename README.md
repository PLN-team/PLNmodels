# PLNmodels: Poisson lognormal models

## Description

The Poisson lognormal model and variants can be used for a variety of multivariate problems when count data are at play (including PCA for count data and network inference). This package implements an efficient algorithm to fit such models accompanied with a set of functions for vizualisation and diagnostic.    
## Installation

```
devtools::install_github("jchiquet/PLNmodels")
```

## Use and example

See the package [vignette](https://github.com/jchiquet/PLNmodels/blob/master/vignettes/trichoptera.Rmd) running PLNPCA and PLNnetwork on the *ade4* Trichoptera data set. To build the vignettes on intallation, you need the *ade4* package installed and to run the following code:

```
devtools::install_github("jchiquet/PLNmodels", build_vignettes=TRUE)
```

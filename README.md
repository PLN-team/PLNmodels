# PLNmodels: Poisson lognormal models <img src="inst/sticker/PLNmodels.png" align="right" width="155" height="180"/>

## Description

The Poisson lognormal model and variants can be used for a variety of multivariate problems when count data are at play (including PCA for count data and network inference). This package implements an efficient algorithm to fit such models accompanied with a set of functions for vizualisation and diagnostic.    

## Installation

For now, only works on Linux. Installation from source on Linux requires libnlopt 2.4-2. On Debian or Ubuntu use libnlopt-dev:

```bash
sudo apt-get install libnlopt-dev
```

On Fedora we need NLopt-devel:


```bash
sudo yum install NLopt-devel
```

Then you can install from github

```
devtools::install_github("jchiquet/PLNmodels")
```

## Use and example

See the package [vignette](https://github.com/jchiquet/PLNmodels/blob/master/vignettes/trichoptera.Rmd) running PLNPCA and PLNnetwork on the *ade4* Trichoptera data set. To build the vignettes on intallation, you need the *ade4* package installed and to run the following code:

```
devtools::install_github("jchiquet/PLNmodels", build_vignettes=TRUE)
```

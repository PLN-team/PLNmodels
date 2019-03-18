
# PLNmodels: Poisson lognormal models <img src="man/figures/logo.png" align="right" width="155" height="180"/>

[![Travis-CI build
status](https://travis-ci.org/jchiquet/PLNmodels.svg?branch=master)](https://travis-ci.org/jchiquet/PLNmodels)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/jchiquet/PLNmodels?branch=master&svg=true)](https://ci.appveyor.com/project/jchiquet/PLNmodels)
[![Codacy
Badge](https://api.codacy.com/project/badge/Grade/c031ad73ccdb4c88ba11dfd74fab1255)](https://www.codacy.com/app/jchiquet/PLNmodels?utm_source=github.com&utm_medium=referral&utm_content=jchiquet/PLNmodels&utm_campaign=Badge_Grade)
[![Coverage
status](https://codecov.io/gh/jchiquet/PLNmodels/branch/master/graph/badge.svg)](https://codecov.io/github/jchiquet/PLNmodels?branch=master)

> The Poisson lognormal model and variants can be used for a variety of
> multivariate problems when count data are at play (including PCA, LDA
> and network inference for count data). This package implements
> efficient algorithms to fit such models accompanied with a set of
> functions for vizualisation and diagnostic.

## Installation

### R Package installation

#### CRAN dependencies

**PLNmodels** needs the following CRAN R packages, so check that their
are installed on your computer.

``` r
library(R6)
library(glassoFast)
library(Matrix)
library(Rcpp)
library(RcppArmadillo)
library(nloptr)
library(igraph)
library(grid)
library(gridExtra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(corrplot)
library(magrittr)
# use install.packages() if needed
```

#### Bioconductor dependencies

**PLNmodels** also needs two BioConductor packages

``` r
library(phyloseq)
library(biomformat)
## use BiocManager::install() if needed
```

#### Installing PLNmodels

  - For the development version, use

<!-- end list -->

``` r
devtools::install_github("jchiquet/PLNmodels")
```

  - For the last tagged release, use

<!-- end list -->

``` r
devtools::install_github("jchiquet/PLNmodels@v0.8.2")
```

  - Windows users may want to use [the binary version of the
    package](https://github.com/jchiquet/PLNmodels/releases/download/v0.8.2/PLNmodels_0.8.2.zip)

## Usage and main fitting functions

The package comes with a ecological data to present the functionality

``` r
library(PLNmodels)
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
```

The main fitting functions work with the usual `R formula` notations,
with mutivariate responses on the left hand side. You probably want to
start by one of them. Check the corresponding vignette and documentation
page. There is a dedicated vignettes for each model in the package (See
<http://jchiquet.github.io/PLNmodels/articles/>).

### Unpenalized Poisson lognormal model (aka PLN)

``` r
myPLN <- PLN(Abundance ~ 1, data = trichoptera)
```

### Rank Contraint Poisson lognormal for Poisson Principal Component Analysis (aka PLNPCA)

``` r
myPCA <- PLNPCA(Abundance ~ 1, data = trichoptera, ranks = 1:8)
```

### Poisson lognormal discriminant analysis (aka PLNLDA)

``` r
myLDA <- PLNLDA(Abundance ~ 1, grouping = Group, data = trichoptera)
```

### Sparse Poisson lognormal model for sparse covariance inference for counts (aka PLNnetwork)

``` r
myPLNnetwork <- PLNnetwork(Abundance ~ 1, data = trichoptera)
```

## References

Please cite our work using the following references:

  - J. Chiquet, M. Mariadassou and S. Robin: Variational inference for
    probabilistic Poisson PCA, the Annals of Applied Statistics, 12:
    2674–2698, 2018. [link](http://dx.doi.org/10.1214/18-AOAS1177)

  - J. Chiquet, M. Mariadassou and S. Robin: Variational inference for
    sparse network reconstruction from count data, arXiv preprint, 2018.
    [link](https://arxiv.org/abs/1806.03120)

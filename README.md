
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

### System Requirements

Installation requires a system version of
[nlopt 2.4-2](https://nlopt.readthedocs.io/)

On **Debian** or **Ubuntu** use `libnlopt-dev`:

``` bash
sudo apt-get install libnlopt-dev
```

On **Debian testing** use `libnlopt-cxx-dev`:

``` bash
sudo apt-get install libnlopt-cxx-dev
```

On **Fedora** or similar use `NLopt-devel`:

``` bash
sudo yum install NLopt-devel
```

With **Mac OS X**, install `nlopt` via [homebrew](https://brew.sh/)

``` bash
brew install nlopt
```

On **Windows**, the package builds and installs correctly by [including
static libraries](https://github.com/rwinlib/nlopt) on compilation.

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

  - For the last tagged release, use

<!-- end list -->

``` r
devtools::install_github("jchiquet/PLNmodels@v0.7.0.2")
```

  - Windows users may want to use [the binary version of the
    package](https://github.com/jchiquet/PLNmodels/releases/download/v0.7.0.2/PLNmodels_0.7.0.2.zip)

  - For the development version, use

<!-- end list -->

``` r
devtools::install_github("jchiquet/PLNmodels")
```

## Usage and main fitting functions

The package comes with a ecological data to present the functionality

``` r
library(PLNmodels)
data(trichoptera)
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
myLDA <- PLNLDA(Abundance ~ 1, grouping = trichoptera$Group, data = trichoptera)
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

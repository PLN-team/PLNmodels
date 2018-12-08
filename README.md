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
> multivariate problems when count data are at play (including PCA for
> count data and network inference). This package implements an
> efficient algorithm to fit such models accompanied with a set of
> functions for vizualisation and diagnostic.

## Installation

### System Requirements

Installation requires a system version of
[nlopt 2.4-2](https://nlopt.readthedocs.io/)

On **Debian** or **Ubuntu** use `libnlopt-dev`:

``` bash
sudo apt-get install libnlopt-dev
```

On **Fedora** or similar use `NLopt-devel`:

``` bash
sudo yum install NLopt-devel
```

With **Mac OS X**, install `nlopt` via [homebrew](https://brew.sh/)

``` bash
brew install nlopt
```

On **Windows**, the package now builds and installs correctly by
[including static libraries](https://github.com/rwinlib/nlopt) on
compilation. For the binary version of the package, [check this
link](https://ci.appveyor.com/project/jchiquet/plnmodels/build/artifacts)

### R Package installation

``` r
## w/o vignettes
devtools::install_github("jchiquet/PLNmodels")
devtools::install_github("jchiquet/PLNmodels", build_vignettes = TRUE)
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
`vignette(package = "PLNmodels")` or
<http://jchiquet.github.io/PLNmodels/articles/>).

### Unpenalized Poisson lognormal model (aka PLN)

``` r
myPLN <- PLN(Abundance ~ 1, data = trichoptera)
```

### Rank Contraint Poisson lognormal for Poisson Principal Component Analysis (ala PLNPCA)

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
    probabilistic Poisson PCA, the Annals of Applied Statistics, to
    appear. [link](https://arxiv.org/abs/1703.06633)

  - J. Chiquet, M. Mariadassou and S. Robin: Variational inference for
    sparse network reconstruction from count data, arXiv preprint, 2018.
    [link](https://arxiv.org/abs/1806.03120)

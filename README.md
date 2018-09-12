
PLNmodels: Poisson lognormal models <img src="inst/sticker/PLNmodels.png" align="right" width="155" height="180"/>
==================================================================================================================

[![Travis-CI build status](https://travis-ci.org/jchiquet/PLNmodels.svg?branch=master)](https://travis-ci.org/jchiquet/PLNmodels)

> The Poisson lognormal model and variants can be used for a variety of multivariate problems when count data are at play (including PCA for count data and network inference). This package implements an efficient algorithm to fit such models accompanied with a set of functions for vizualisation and diagnostic. Learn more in the vignettes: `vignette(package = "PLNmodels")`.

Installation
------------

###  System Requirements

Installation requires a system version of [nlopt 2.4-2](https://nlopt.readthedocs.io/)

-   On Debian or Ubuntu use `libnlopt-dev`:

``` bash
sudo apt-get install libnlopt-dev
```

-   On Fedora or similar use `NLopt-devel`:

``` bash
sudo yum install NLopt-devel
```

-   With Mac OS X, install `nlopt` via [homebrew](https://brew.sh/)

``` bash
brew install nlopt
```

-   Windows is not supported yet: feel free to give it a try and send us some feedbacks.

###  R Package installation

``` r
## w/o vignettes
devtools::install_github("jchiquet/PLNmodels")
devtools::install_github("jchiquet/PLNmodels", build_vignettes = TRUE)
```

References
----------

Please cite our work using the following references:

-   J. Chiquet, M. Mariadassou and S. Robin: Variational inference for probabilistic Poisson PCA, the Annals of Applied Statistics, to appear. [link](https://arxiv.org/abs/1703.06633)

-   J. Chiquet, M. Mariadassou and S. Robin: Variational inference for sparse network reconstruction from count data, arXiv preprint, 2018. [link](https://arxiv.org/abs/1806.03120)


# PLNmodels: Poisson lognormal models <img src="man/figures/logo.png" align="right" width="155" height="180"/>

<!-- badges: start -->

[![R build
status](https://github.com/pln-team/PLNmodels/workflows/R-CMD-check/badge.svg)](https://github.com/pln-team/PLNmodels/actions)
[![Coverage
status](https://codecov.io/gh/pln-team/PLNmodels/branch/master/graph/badge.svg)](https://codecov.io/github/pln-team/PLNmodels?branch=master)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/PLNmodels)](https://cran.r-project.org/package=PLNmodels)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![](https://img.shields.io/github/last-commit/pln-team/PLNmodels.svg)](https://github.com/pln-team/PLNmodels/commits/master)
<!-- badges: end -->

> The Poisson lognormal model and variants can be used for a variety of
> multivariate problems when count data are at play (including PCA, LDA
> and network inference for count data). This package implements
> efficient algorithms to fit such models accompanied with a set of
> functions for visualization and diagnostic. See [this deck of
> slides](https://pln-team.github.io/slideshow/) for a
> comprehensive introduction.

## Installation

**PLNmodels** is available on
[CRAN](https://cran.r-project.org/package=PLNmodels). The development
version is available on [Github](https://github.com/pln-team/PLNmodels).

### R Package installation

<!-- #### CRAN dependencies -->
<!-- **PLNmodels** needs the following CRAN R packages, so check that they are are installed on your computer. -->
<!-- ```{r CRAN dependencies, eval = FALSE} -->
<!-- required_CRAN <- c("R6", "glassoFast", "Matrix", "Rcpp", "RcppArmadillo", -->
<!--                    "nloptr", "igraph", "grid", "gridExtra", "dplyr", -->
<!--                    "tidyr", "ggplot2", "corrplot", "magrittr", "devtools") -->
<!-- not_installed_CRAN <- setdiff(required_CRAN, rownames(installed.packages())) -->
<!-- if (length(not_installed_CRAN) > 0) install.packages(not_installed_CRAN) -->
<!-- ``` -->
<!-- #### Bioconductor dependencies -->
<!-- **PLNmodels** also needs two BioConductor packages -->
<!-- ```{r Bioconductor dependencies, eval = FALSE} -->
<!-- required_BioC <- c("phyloseq", "biomformat") -->
<!-- not_installed_BioC <- setdiff(required_BioC, rownames(installed.packages())) -->
<!-- if (length(not_installed_BioC) > 0) BiocManager::install(not_installed_BioC) -->
<!-- ``` -->

#### Installing PLNmodels

-   For the last stable version, use the CRAN version

``` r
install.packages("PLNmodels")
```

-   For the development version, use the github install

``` r
remotes::install_github("pln-team/PLNmodels")
```

-   For a specific tagged release, use

``` r
remotes::install_github("pln-team/PLNmodels@tag_number")
```

## Usage and main fitting functions

The package comes with an ecological data set to present the
functionality

``` r
library(PLNmodels)
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
```

The main fitting functions work with the usual `R formula` notations,
with mutivariate responses on the left hand side. You probably want to
start by one of them. Check the corresponding vignette and documentation
page. There is a dedicated vignettes for each model in the package (See
<https://pln-team.github.io/PLNmodels/articles/>).

### Unpenalized Poisson lognormal model (aka PLN)

``` r
myPLN <- PLN(Abundance ~ 1, data = trichoptera)
```

### Rank Constrained Poisson lognormal for Poisson Principal Component Analysis (aka PLNPCA)

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

### Mixture of Poisson lognormal models for model-based clustering of counts (aka PLNmixture)

``` r
myPLNmixture <- PLNmixture(Abundance ~ 1, data = trichoptera)
```

## References

Please cite our work using the following references:

-   J. Chiquet, M. Mariadassou and S. Robin: The Poisson-lognormal model
    as a versatile framework for the joint analysis of species
    abundances, Frontiers in Ecology and Evolution, 2021.
    [link](https://www.frontiersin.org/articles/10.3389/fevo.2021.588292/full)

-   J. Chiquet, M. Mariadassou and S. Robin: Variational inference for
    sparse network reconstruction from count data, Proceedings of the
    36th International Conference on Machine Learning (ICML), 2019.
    [link](http://proceedings.mlr.press/v97/chiquet19a.html)

-   J. Chiquet, M. Mariadassou and S. Robin: Variational inference for
    probabilistic Poisson PCA, the Annals of Applied Statistics, 12:
    2674–2698, 2018. [link](http://dx.doi.org/10.1214/18%2DAOAS1177)

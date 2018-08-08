# PLNmodels: Poisson lognormal models <img src="inst/sticker/PLNmodels.png" align="right" width="155" height="180"/>

## Description

> The Poisson lognormal model and variants can be used for a variety of multivariate problems when count data are at play (including PCA for count data and network inference). This package implements an efficient algorithm to fit such models accompanied with a set of functions for vizualisation and diagnostic.    

## Installation

[![Build Status](https://travis-ci.org/jchiquet/PLNmodels.svg?branch=master)](https://travis-ci.org/jchiquet/PLNmodels)

### Linux

Installation from source on Linux requires libnlopt 2.4-2. On Debian or Ubuntu use libnlopt-dev:

```bash
sudo apt-get install libnlopt-dev
```

On Fedora we need NLopt-devel:


```bash
sudo yum install NLopt-devel
```

Then you can install from github

```R
devtools::install_github("jchiquet/PLNmodels")
```

### MAC OS

Installation requires nlopt, which can be installed via [homebrew](https://brew.sh/)

```bash
brew install nlopt
```

Finally install the package via

```R
devtools::install_github("jchiquet/PLNmodels")
```

#### Troubleshooting

If you experience problems due to the lack of OpenMP support in Clang like 
```bash
clang: error: unsupported option '-fopenmp'
make: *** [RcppExports.o] Error 1
ERROR: compilation failed for package ‘PLNmodels’
```
[have a look a this page](http://thecoatlessprofessor.com/programming/openmp-in-r-on-os-x/) or this [issue](https://github.com/jchiquet/PLNmodels/issues/12)

A simple fix (for R 3.4.x) consists in installing the [dev tool chains](https://github.com/coatless/r-macos-rtools/releases) for Mac OS X: https://github.com/coatless/r-macos-rtools/releases/download/v1.1.0/macos-rtools-1.1.0.pkg

For R 3.5.x, you can also install the CRAN Mac OS toolchain from https://cran.r-project.org/bin/macosx/tools/, especially the [Clang 6.0.0 compiler](https://cran.r-project.org/bin/macosx/tools/clang-6.0.0.pkg) and the [GNU Fortran 6.1 compiler](https://cran.r-project.org/bin/macosx/tools/gfortran-6.1.pkg) and modify your `Makevars` to make it look like

```bash
F77 = /usr/local/bin/gfortran
FC = $F77
CXX =  /usr/local/clang6/bin/clang++  -Wall
LDFLAGS=-L/usr/local/clang6/lib
CC=  /usr/local/clang6/bin/clang
SHLIB_CXXLD=ccache /usr/local/clang6/bin/clang++
CXX11 =  /usr/local/clang6/bin/clang++
CXX98 =  /usr/local/clang6/bin/clang++
CXX14 =  /usr/local/clang6/bin/clang++
```

### Windows

Not supported yet...

## Use and example

See the package [vignette](https://github.com/jchiquet/PLNmodels/blob/master/vignettes/trichoptera.Rmd) running PLNPCA and PLNnetwork on the *ade4* Trichoptera data set. To build the vignettes on intallation, you need the *ade4* package installed and to run the following code:

```
devtools::install_github("jchiquet/PLNmodels", build_vignettes=TRUE)
```

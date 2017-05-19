##' PLNmodels
##'
##' Implements probabilistic PCA for count data via a Poisson log-normal model. Can use covariates
##' and offset on top of the main observation matrix in a GLM setup. The algorithm relies on variational
##' inference which is performed by maximizing a bi-convex function with the L-BFGS-B method.
##'
##' Try vignette("trichoptera", package="PLNmodels") for an (not yet comprehensive) example.
##'
##' The main function is \code{\link{PLNPCA}} which produces two kinds of Reference Class objects,
##' \code{\link[=PLNfamily-class]{PLNfamily}} and \code{\link[=PLNfit.PCA-class]{PLNfit.PCA}}. See their
##' documentation and the associated methods for manipulation.
##'
##' @author Julien Chiquet \email{julien.chiquet@@inra.fr}
##' @author Mahendra Mariadassou \email{mahendra.mariadassou@@inra.fr}
##' @author Stéphane Robin \email{stephane.robin@@inra.fr}
##'
##' @import methods parallel ggplot2 reshape2 grid gridExtra R6
##' @docType package
##' @name PLNmodels
NULL

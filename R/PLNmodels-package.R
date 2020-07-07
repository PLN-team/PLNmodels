##' PLNmodels
##'
##'The Poisson lognormal model and variants can be used for a variety of multivariate problems when count data are at play (including
##' PCA or LDA for count data, network inference). This package implements efficient
##' variational algorithms to fit such models accompanied with a set of functions for visualization and diagnostic.
##'
##' @section Multivariate Poisson lognormal model (aka PLN):
##'
##' See the main function [PLN()] and the associated methods for manipulation.
##'
##' Also try vignette("PLN_trichoptera", package="PLNmodels") for an overview.
##'
##' @section Rank Constrained Poisson lognormal for Poisson Principal Component Analysis (aka PLNPCA):
##'
##' See the main function [PLNPCA()] and the associated methods for manipulation.
##'
##' The Poisson PCA and the associated variational inference is fully explained in Chiquet el al (2018), see reference below.
##'
##' Also try vignette("PLNPCA_trichoptera", package="PLNmodels") for an overview.
##'
##' @section Sparse Poisson lognormal model for sparse covariance inference for counts (aka PLNnetwork):
##'
##' See the main function [PLNnetwork()] and the associated methods for manipulation.
##'
##' Also try vignette("PLNnetwork_trichoptera", package="PLNmodels") for an overview.
##'
##' @section Poisson lognormal discriminant analysis (aka PLNLDA):
##'
##' See the main function [PLNLDA()] and the associated methods for manipulation.
##'
##' Also try vignette("PLNLDA_trichoptera", package="PLNmodels") for an overview.
##'
##' @author Julien Chiquet \email{julien.chiquet@@inrae.fr}
##' @author Mahendra Mariadassou \email{mahendra.mariadassou@@inrae.fr}
##' @author Stéphane Robin \email{stephane.robin@@inrae.fr}
##'
##' @import methods R6 parallel Matrix nloptr
##' @importFrom Rcpp sourceCpp
##' @import dplyr
##' @importFrom tidyr gather
##' @import ggplot2
##' @importFrom igraph V E plot.igraph layout_in_circle graph_from_adjacency_matrix degree delete.vertices
##' @importFrom corrplot corrplot
##' @useDynLib PLNmodels
##' @docType package
##' @name PLNmodels
NULL

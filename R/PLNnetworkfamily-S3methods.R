## =========================================================================================
##
## PUBLIC S3 METHODS FOR PLNnetworkfamily + OTHER FUNCTIONS
##
## =========================================================================================


## Auxiliary functions to check the given class of an objet
isPLNnetworkfamily <- function(Robject) {inherits(Robject, "PLNnetworkfamily")}

#' Display various ouputs (goodness-of-fit criteria, robustness, diagnostic) associated with a collection of PLNnetwork fits (a PLNnetworkfamily)
#'
#' @name plot.PLNnetworkfamily
#'
#' @param x an R6 object with class PLNfamily
#' @param type a character, either "criteria", "stability" or "diagnostic" for the type of plot.
#' @param criteria vector of characters. The criteria to plot in c("loglik", "BIC", "ICL", "R_squared", "EBIC", "pen_loglik").
#' Default is  c("loglik", "pen_loglik", "BIC", "EBIC"). Only relevant when type equals "criteria".
#' @param log.x logical: should the x-axis be repsented in log-scale? Default is \code{TRUE}.
#' @param stability scalar: the targeted level of stability in stability plot. Default is .9.
#' @param annotate logical: should the value of approximated R squared be added to the plot of criteria? Default is \code{TRUE}.
#' @param ... additional parameters for S3 compatibility. Not used
#'
#' @return Produces a plot  representing the evolution of the criteria of the different models considered,
#' highlighting the best model in terms of BIC and EBIC.
#'
#' @export
plot.PLNnetworkfamily <-
  function(x,
           type     = c("criteria", "stability", "diagnostic"),
           criteria = c("loglik", "pen_loglik", "BIC", "EBIC"),
           log.x    = TRUE,
           stability = 0.9,
           annotate = TRUE, ...) {
  stopifnot(isPLNnetworkfamily(x))
  type <- match.arg(type)
  if (type == "criteria")
    p <- x$plot(criteria, annotate)
  if (type == "stability")
    p <- x$plot_stars(stability, log.x)
  if (type == "diagnostic")
    p <- x$plot_objective()

  p
}

#' @describeIn getModel Model extraction for PLNnetworkfamily
#' @export
getModel.PLNnetworkfamily <- function(Robject, var, index = NULL) {
  stopifnot(isPLNnetworkfamily(Robject))
  Robject$getModel(var, index = NULL)
}

#' @describeIn getBestModel Model extraction for PLNnetworkfamily
#' @export
getBestModel.PLNnetworkfamily <- function(Robject, crit = c("BIC", "loglik", "R_squared", "EBIC", "StARS"), ...) {
  stopifnot(isPLNnetworkfamily(Robject))
  stability <- list(...)[["stability"]]
  if (is.null(stability)) stability <- 0.9
  Robject$getBestModel(match.arg(crit), stability)
}


#' Extract the regularization path of a PLNnetwork fit
#'
#' @name coefficient_path
#' @param Robject an object with class PLNnetworkfamily, i.e. an output from \code{\link{PLNnetwork}}
#' @param precision a logical, should the coefficients of the precision matrix Omega or the covariance matrice Sigma be sent back. Default is \code{TRUE}.
#' @param corr a logical, should the correlation (partial in case  precision = TRUE) be sent back. Default is \code{TRUE}.
#'
#' @return  Send back a tibble/data.frame.
#'
#' @export
coefficient_path <- function(Robject, precision = TRUE, corr = TRUE) {
  stopifnot(isPLNnetworkfamily(Robject))
  Robject$coefficient_path(precision, corr)
}

#' Compute the stability path by stability selection
#'
#' @name stability_selection
#' @param Robject an object with class PLNnetworkfamily, i.e. an output from \code{\link{PLNnetwork}}
#' @param subsamples a list of vectors describing the subsamples. The number of vectors (or list length) determines th number of subsamples used in the stability selection. Automatically set to 20 subsamples with size \code{10*sqrt(n)} if \code{n >= 144} and \code{0.8*n} otherwise following Liu et al. (2010) recommandations.
#' @param control a list controling the main optimization process in each call to PLNnetwork. See \code{\link[=PLNnetwork]{PLNnetwork}} for details.
#' @param mc.cores the number of cores to used. Default is 1.
#'
#' @return the list of subsamples. The estimated probabilities of selection of the edges are stored in the fields stability_path of the initial Robject with class PLNnetworkfamily
stability_selection <- function(Robject, subsamples = NULL, control = list(), mc.cores = 1) {
  stopifnot(isPLNnetworkfamily(Robject))
  Robject$stability_selection(subsamples, control, mc.cores)
}

## =========================================================================================
##
## PUBLIC S3 METHODS FOR PLNnetworkfamily
## =========================================================================================


## Auxiliary functions to check the given class of an objet
isPLNnetworkfamily <- function(Robject) {inherits(Robject, "PLNnetworkfamily")}

#' Display the criteria associated with a collection of PLNnetwork fits (a PLNnetworkfamily)
#'
#' @name plot.PLNnetworkfamily
#'
#' @param x an R6 object with class PLNfamily
#' @param criteria vector of characters. The criteria to plot in c("loglik", "BIC", "ICL", "R_squared", "EBIC", "pen_loglik").
#' Default is  c("loglik", "pen_loglik", "BIC", "EBIC").
#' @param log.x logical: should the x-axis be repsented in log-scale? Default is \code{TRUE}.
#' @param annotate logical: should the value of approximated R squared be added to the plot? Default is \code{TRUE}.
#' @param ... additional parameters for S3 compatibility. Not used
#'
#' @return Produces a plot  representing the evolution of the criteria of the different models considered,
#' highlighting the best model in terms of BIC and EBIC.
#'
#' @export
plot.PLNnetworkfamily <- function(x, criteria = c("loglik", "pen_loglik", "BIC", "EBIC"), log.x = TRUE, annotate = TRUE,...) {
  stopifnot(isPLNnetworkfamily(x))
  x$plot(criteria, annotate)
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



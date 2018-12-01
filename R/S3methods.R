## =========================================================================================
##
## PUBLIC S3 METHODS FOR THE USERS
##
## =========================================================================================

## Auxiliary functions to check the given class of an objet
isPLNfit           <- function(Robject) {inherits(Robject, "PLNfit"          )}
isPLNLDAfit        <- function(Robject) {inherits(Robject, "PLNLDAfit"       )}
isPLNfamily        <- function(Robject) {inherits(Robject, 'PLNfamily'       )}
isPLNnetworkfamily <- function(Robject) {inherits(Robject, "PLNnetworkfamily")}
isPLNPCAfamily     <- function(Robject) {inherits(Robject, "PLNPCAfamily"    )}

## =========================================================================================
##
## METHODS FOR PLNfit
##
## =========================================================================================

#' Predict counts of a new sample
#'
#' @name predict.PLNfit
#'
#' @param object an R6 object with class PLNfit
#' @param newdata A data frame in which to look for variables and offsets with which to predict
#' @param type The type of prediction required. The default is on the scale of the linear predictors (i.e. log average count)
#' @param ... additional parameters for S3 compatibility. Not used
#' @return A matrix of predicted log-counts (if type = "link") or predicted counts (if type = "response").
#' @export
predict.PLNfit <- function(object, newdata, type = c("link", "response"), ...) {
  stopifnot(isPLNfit(object))
  object$predict(newdata, type, parent.frame())
}

#' Extracts model coefficients from objects returned by \code{\link[=PLN]{PLN}} and its variants
#'
#' @name coef.PLNfit
#'
#' @param object an R6 object with class PLNfit
#' @param ... additional parameters for S3 compatibility. Not used
#' @return A matrix of coefficients extracted from the PLNfit model.
#'
#' @export
coef.PLNfit <- function(object, ...) {
  stopifnot(isPLNfit(object))
  object$model_par$Theta
}

#' Extracts model covariance from objects returned by \code{\link[=PLN]{PLN}} and its variants
#'
#' @name vcov.PLNfit
#'
#' @param object an R6 object with class PLNfit
#' @param ... additional parameters for S3 compatibility. Not used
#' @return A matrix of variance/covariance extracted from the PLNfit model.
#'
#' @export
vcov.PLNfit <- function(object, ...) {
  stopifnot(isPLNfit(object))
  object$model_par$Sigma
}

## =========================================================================================
##
## METHODS FOR PLNLDAfit
##
## =========================================================================================

#' Predict group of new samples
#'
#' @name predict.PLNLDAfit
#'
#' @param object an R6 object with class PLNLDAfit
#' @param newdata A data frame in which to look for variables with which to predict.
#' @param newOffsets A matrix in which to look for offsets with which to predict.
#' @param newCounts A matrix in which to look for counts with to predict
#' @param type The type of prediction required. The default are posterior probabilities for each group (in log-scale),
#'             the alternative "response" is the group with maximal posterior probability.
#' @param control a list for controlling the optimization. See \code{\link[=PLN]{PLN}} for details.
#' @param ... additional parameters for S3 compatibility. Not used
#' @return A matrix of predicted scores for each group (if type = "score") or a vector of predicted
#'         groups (if type = "response").
#' @export
predict.PLNLDAfit <- function(object, newdata, newOffsets, newCounts,
                              type = c("posterior", "response"), control = list(), ...) {
  stopifnot(isPLNLDAfit(object))
  object$predict(newdata, newOffsets, newCounts, type, control)
}

## =========================================================================================
##
## METHODS FOR PLNfamily
##
## =========================================================================================

## S3 methods declared generic for children of PLNfamily

#' Best model extraction from a collection of models
#'
#' @param Robject an object with class PLNPCAfamilly ot PLNnetworkfamily
#' @param crit a character for the criterion used to performed the selection. Either
#' "BIC", "ICL", "EBIC", "StARS", "R_squared". Default is \code{ICL} for PLNPCA, and \code{BIC} for PLNnetwork.
#'  If StARS (Stability Approach to Regularization Selection) is chosen and stability selection
#'  was not yet performed, the function will call the method stability_selection with default argument.
#' @param ... additional parameters for StARS criterion (only for PLNnetwork). stability, a scalar indicating the target stability (= 1 - 2 beta) at which the network is selected. Default is \code{0.9}.
#' @return  Send back a object with class \code{\link[=PLNPCAfit]{PLNPCAfit}} or \code{\link[=PLNnetworkfit]{PLNnetworkfit}}
#'
#' @export
getBestModel <- function(Robject, crit, ...) {UseMethod("getBestModel", Robject)}

#' Model extraction from a collection of models
#'
#' @param Robject an R6 object with class PLNPCAfamily or PLNnetworkfamily
#' @param var value of the parameter (rank for PLNPCA, sparisty for PLNnetwork) that identifies the model to be extracted from the collection. If no exact match is found, the model with closest parameter value is returned with a warning.
#' @param index Integer index of the model to be returned. Only the first value is taken into account.
#'
#' @return Sends back a object with class \code{\link[=PLNPCAfit]{PLNPCAfit}} or \code{\link[=PLNPCAfit]{PLNnetworkfit}}.
#'
#' @export
getModel <- function(Robject, var , index) {UseMethod("getModel"    , Robject)}

## =========================================================================================
##
## METHODS FOR PLNPCAfamily
##
## =========================================================================================

#' Display the criteria associated with a collection of PLNPCA fits (a PLNPCAfamily)
#'
#' @name plot.PLNPCAfamily
#'
#' @param x an R6 object with class PLNfamily
#' @param criteria vector of characters. The criteria to plot in c("loglik", "BIC", "ICL", "R_squared").
#' Default is  c("loglik", "BIC", "ICL").
#' @param annotate logical: should the value of approximated R squared be added to the plot?
#' @param ... additional parameters for S3 compatibility. Not used
#'
#' @return Produces a plot  representing the evolution of the criteria of the different models considered,
#' highlighting the best model in terms of BIC and ICL.
#'
#' @export
plot.PLNPCAfamily <- function(x, criteria = c("loglik", "BIC", "ICL"), annotate = TRUE, ...) {
  stopifnot(isPLNfamily(x))
  x$plot(criteria, annotate)
}

#' @describeIn getModel Model extraction for PLNPCAfamily
#' @export
getModel.PLNPCAfamily <- function(Robject, var, index = NULL) {
  stopifnot(isPLNPCAfamily(Robject))
  Robject$getModel(var, index = NULL)
}

#' @describeIn getBestModel Model extraction for PLNPCAfamily
#' @export
getBestModel.PLNPCAfamily <- function(Robject, crit = c("ICL", "BIC", "R_squared"), ...) {
  stopifnot(isPLNPCAfamily(Robject))
  Robject$getBestModel(match.arg(crit))
}

## =========================================================================================
##
## METHODS FOR PLNnetworkfamily
## =========================================================================================

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



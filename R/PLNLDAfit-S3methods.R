## =========================================================================================
##
## PUBLIC S3 METHODS FOR PLNLDAfit
##
## =========================================================================================

## Auxiliary functions to check the given class of an objet
isPLNLDAfit <- function(Robject) {inherits(Robject, "PLNLDAfit"       )}

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

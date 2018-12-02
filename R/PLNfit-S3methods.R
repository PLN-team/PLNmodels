## =========================================================================================
##
## PUBLIC S3 METHODS FOR PLNfit
##
## =========================================================================================

## Auxiliary functions to check the given class of an objet
isPLNfit <- function(Robject) {inherits(Robject, "PLNfit"          )}

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

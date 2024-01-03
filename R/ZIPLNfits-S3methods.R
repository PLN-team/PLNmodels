## =========================================================================================
##
## PUBLIC S3 METHODS FOR ZIPLNfit
##
## =========================================================================================

## Auxiliary functions to check the given class of an objet
isZIPLNfit <- function(Robject) {inherits(Robject, "ZIPLNfit")}

# #' Predict counts of a new sample
# #'
# #' @name predict.ZIPLNfit
# #'
# #' @param object an R6 object with class [`ZIPLNfit`]
# #' @param newdata A data frame in which to look for variables and offsets with which to predict
# #' @param type The type of prediction required. The default is on the scale of the linear predictors (i.e. log average count)
# #' @param ... additional parameters for S3 compatibility. Not used
# #' @return A matrix of predicted log-counts (if `type = "link"`) or predicted counts (if `type = "response"`).
# #' @importFrom stats predict
# #' @export
# predict.ZIPLNfit <- function(object, newdata, type = c("link", "response"), ...) {
#   stopifnot(isZIPLNfit(object))
#   object$predict(newdata, type, parent.frame())
# }

#' Extract model coefficients
#'
#' @description Extracts model coefficients from objects returned by [ZIPLN()] and its variants
#'
#' @name coef.ZIPLNfit
#'
#' @param object an R6 object with class [`ZIPLNfit`]
#' @param type type of parameter that should be extracted. Either "mainPLN" (default) for \deqn{\Theta},
#' "mainZI" (default) for \deqn{\Theta0} or "precision" for \deqn{\Omega}
#' @param ... additional parameters for S3 compatibility. Not used
#' @return A matrix of coefficients extracted from the ZIPLNfit model.
#'
#' @importFrom stats coef
#' @seealso [sigma.ZIPLNfit()]
#'
#' @examples
#' data(scRNA)
#' # data subsample: only 100 random cell and the 50 most varying transcript
#' subset <- sample.int(nrow(scRNA), 100)
#' myPLN  <- ZIPLN(counts[, 1:50] ~ 1 + offset(log(total_counts)), subset = subset, data = scRNA)
#'
#' @export
coef.ZIPLNfit <- function(object, type = c("mainPLN", "mainZI", "precision"), ...) {
  stopifnot(isZIPLNfit(object))
  switch(match.arg(type),
         mainPLN   = object$model_par$Theta,
         mainZI    = object$model_par$Theta0,
         precision = object$model_par$Omega)
}

#' Extracts model fitted values from objects returned by [ZIPLN()] and its variants
#'
#' @name fitted.ZIPLNfit
#'
#' @inheritParams coef.ZIPLNfit
#' @return A matrix of Fitted values extracted from the object object.
#'
#' @importFrom stats fitted
#' @export
fitted.ZIPLNfit <- function(object, ...) {
  stopifnot(isZIPLNfit(object))
  object$fitted
}

#' Extract variance-covariance of residuals 'Sigma'
#'
#' @name sigma.ZIPLNfit
#' @description Extract the variance-covariance matrix of the residuals, usually noted \deqn{\Sigma} in ZIPLN models.
#'
#' @inheritParams coef.ZIPLNfit
#'
#' @return A semi definite positive matrix of size p, assuming there are p species in the model.
#'
#' @export
#'
#' @importFrom stats sigma
sigma.ZIPLNfit <- function(object, ...) {
  stopifnot(isZIPLNfit(object))
  solve(object$model_par$Omega)
}


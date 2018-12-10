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

#' Extracts (one-data) Fisher information matrix of Theta from objects returned by \code{\link[=PLN]{PLN}} and its variants.
#' Fisher matrix is computed using one of two approximation scheme: wald (default, conservative, gives large confidence interval)
#' or louis (anticonservative).
#'
#' @name fisher.PLNfit
#'
#' @param object an R6 object with class PLNfit
#' @param type Either `wald` (default) or `louis`. Approxomation scheme used to compute the
#' Fisher information matrix
#' @param ... additional parameters for S3 compatibility. Not used
#' @return A block-diagonal matrix with p (number of species) blocks of size d (number of covariates), assuming
#' Theta is a matrix of size d * p.
#'
#' @export
fisher.PLNfit <- function(object, type = c("wald", "louis"), ...) {
  stopifnot(isPLNfit(object))
  type <- match.arg(type)
  if (type != object$fisher$type) {
    stop(paste("Fisher information was not computed using the", type, "approximation. Try another approximation scheme."))
  }
  object$fisher$mat
}

#' Extracts univariate standard errors for the estimated coefficient of Theta. Standard errors are computed from the (approximate)
#' Fisher information matrix. See \code{\link[=fisher.PLNfit]{fisher}} for more details on the approximations.
#'
#' @name standard_error.PLNfit
#'
#' @param object an R6 object with class PLNfit
#' @param type Either `Wald` (default) or `Louis`. Approxomation scheme used to compute the
#' Fisher information matrix
#' @param ... additional parameters for S3 compatibility. Not used
#'
#' @return A p * d positive matrix (same size as Theta) with standard errors.
#' @export
#'
standard_error.PLNfit <- function(object, type = c("wald", "louis"), ...) {
  stopifnot(isPLNfit(object))
  type <- match.arg(type)
  if (type != object$fisher$type) {
    stop(paste("Standard errors were not computed using the", type, "approximation. Try another approximation scheme."))
  }
  object$std_error
}

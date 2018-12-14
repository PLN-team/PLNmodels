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

#' Extracts model fitted values from objects returned by \code{\link[=PLN]{PLN}} and its variants
#'
#' @name fitted.PLNfit
#'
#' @param object an R6 object with class PLNfit
#' @param ... additional parameters for S3 compatibility. Not used
#' @return A matrix of Fitted values extracted from the object object.
#'
#' @export
fitted.PLNfit <- function(object, ...) {
  stopifnot(isPLNfit(object))
  object$fitted
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

# Fisher information matrix for Theta
#
# @description Compute or extract component-wise standard errors of Theta in multivariate (generalized) linear models of the form \deqn{g(E[Y|X]) = X\Theta}.
# Useful to create confidence intervals and (multivariate) confidence regions under a Gaussian approximation of \eqn{\Theta}. Note that the Fisher information matrix is the one-data version (not scaled by the number of observations).
#
# @param object An `R` object. Currently there are methods for \code{\link{PLNfit}} (and its variants) objects.
# @param ... Further arguments passed to or from other methods.
#
# @return The fisher information matrix of \eqn{\Theta} in generalized linear models.
# @export
#

#' Fisher information matrix for Theta
#'
#' @description Extracts Fisher information matrix of \eqn{\Theta} from objects returned by \code{\link[=PLN]{PLN}} and its variants. Fisher matrix is computed using one of two approximation scheme: wald (default, conservative, gives large confidence interval) or louis (anticonservative). Note that the Fisher information matrix is the one-data version (not scaled by the number of observations).
#'
#' @param object an R6 object with class PLNfit
#' @param type Either `wald` (default) or `louis`. Approxomation scheme used to compute the
#' Fisher information matrix
#' @return A block-diagonal matrix with p (number of species) blocks of size d (number of covariates), assuming
#' \eqn{\Theta} is a matrix of size d * p.
#'
#' @seealso \code{\link[=standard_error.PLNfit]{standard_error}} for standard errors
#'
#' @export
fisher <- function(object, type) {
  UseMethod("fisher", object)
}

#' @describeIn fisher Fisher information matrix for PLNfit
#' @export
fisher.PLNfit <- function(object, type = c("wald", "louis")) {
  stopifnot(isPLNfit(object))
  type <- match.arg(type)
  if (type != object$fisher$type) {
    stop(paste("Fisher information was not computed using the", type, "approximation. Try another approximation scheme."))
  }
  object$fisher$mat
}

# Component-wise standard errors of Theta
#
# @description Compute or extract component-wise standard errors of Theta in multivariate (generalized) linear models of the form \deqn{g(E[Y|X]) = X\Theta}. Useful to compute Z-scores and p-values under a gaussian/student approximation of \eqn{\Theta}
#
# @param object An `R` object. Currently there are methods for \code{\link{PLNfit}} (and its variants) objects.
# @param ... Further arguments passed to or from other methods.
#
# @return The standard errors associated with coefficients of \eqn{\Theta}
# @export
#

#' Component-wise standard errors of Theta
#'
#' @description Extracts univariate standard errors for the estimated coefficient of Theta. Standard errors are computed from the (approximate) Fisher information matrix. See \code{\link[=fisher.PLNfit]{fisher}} for more details on the approximations.
#'
#' @param object an R6 object with class PLNfit
#' @param type Either `Wald` (default) or `Louis`. Approximation scheme used to compute the Fisher information matrix
#'
#' @seealso \code{\link[=fisher.PLNfit]{fisher}} for the complete Fisher information matrix
#'
#' @return A p * d positive matrix (same size as \eqn{\Theta}) with standard errors for the coefficients of \eqn{\Theta}
#'
#' @export
standard_error <- function(object, type) {
  UseMethod("standard_error", object)
}

#' @describeIn standard_error Component-wise standard errors of Theta in PLNfit
#' @export
standard_error.PLNfit <- function(object, type = c("wald", "louis")) {
  stopifnot(isPLNfit(object))
  type <- match.arg(type)
  if (type != object$fisher$type) {
    stop(paste("Standard errors were not computed using the", type, "approximation. Try another approximation scheme."))
  }
  object$std_err
}

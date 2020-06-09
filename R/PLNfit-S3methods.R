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
#' @param object an R6 object with class [`PLNfit`]
#' @param newdata A data frame in which to look for variables and offsets with which to predict
#' @param type The type of prediction required. The default is on the scale of the linear predictors (i.e. log average count)
#' @param ... additional parameters for S3 compatibility. Not used
#' @return A matrix of predicted log-counts (if `type = "link"`) or predicted counts (if `type = "response"`).
#' @export
predict.PLNfit <- function(object, newdata, type = c("link", "response"), ...) {
  stopifnot(isPLNfit(object))
  object$predict(newdata, type, parent.frame())
}

#' Extract model coefficients
#'
#' @description Extracts model coefficients from objects returned by [PLN()] and its variants
#'
#' @name coef.PLNfit
#'
#' @param object an R6 object with class [`PLNfit`]
#' @param type type of parameter that should be extracted. Either "main" (default) for \deqn{\Theta} or "covariance" for \deqn{\Sigma}
#' @param ... additional parameters for S3 compatibility. Not used
#' @return A matrix of coefficients extracted from the PLNfit model.
#'
#' @seealso [sigma.PLNfit()], [vcov.PLNfit()], [standard_error.PLNfit()]
#'
#' @export
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPLN <- PLN(Abundance ~ 1 + offset(log(Offset)), data = trichoptera)
#' coef(myPLN) ## Theta
#' coef(myPLN, type = "covariance") ## Sigma
coef.PLNfit <- function(object, type = c("main", "covariance"), ...) {
  stopifnot(isPLNfit(object))
  switch(match.arg(type),
         main       = object$model_par$Theta,
         covariance = object$model_par$Sigma)
}

#' Extracts model fitted values from objects returned by [PLN()] and its variants
#'
#' @name fitted.PLNfit
#'
#' @inheritParams coef.PLNfit
#' @return A matrix of Fitted values extracted from the object object.
#'
#' @export
fitted.PLNfit <- function(object, ...) {
  stopifnot(isPLNfit(object))
  object$fitted
}

#' Calculate Variance-Covariance Matrix for a fitted [PLN()] model object
#'
#' @name vcov.PLNfit
#'
#' @description Returns the variance-covariance matrix of the main parameters of a fitted [PLN()] model object. The main parameters of the model correspond to \deqn{\Theta}, as returned by [coef.PLNfit()]. The function can also be used to return the variance-covariance matrix of the residuals. The latter matrix can also be accessed via [sigma.PLNfit()]
#'
#' @inheritParams coef.PLNfit
#' @return A matrix of variance/covariance extracted from the PLNfit model. If type="main" and \eqn{\Theta} is a matrix of size d * p, the result is a block-diagonal matrix with p (number of species) blocks of size d (number of covariates). if type="main", it is a symmetric matrix of size p.
#' .
#'
#' @seealso [sigma.PLNfit()], [coef.PLNfit()], [standard_error.PLNfit()]
#'
#' @export
#'
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPLN <- PLN(Abundance ~ 1 + offset(log(Offset)), data = trichoptera)
#' vcov(myPLN) ## variance-covariance of Theta
#' vcov(myPLN, type = "covariance") ## Sigma
vcov.PLNfit <- function(object, type = c("main", "covariance"), ...) {
  stopifnot(isPLNfit(object))
  switch(match.arg(type),
         main       = object$fisher$mat,
         covariance = object$model_par$Sigma)
}


#' Extract variance-covariance of residuals 'Sigma'
#'
#' @name sigma.PLNfit
#' @description Extract the variance-covariance matrix of the residuals, usually noted \deqn{\Sigma} in PLN models. This captures the correlation between the species in the latent space.
#'
#' @inheritParams coef.PLNfit
#'
#' @return A semi definite positive matrix of size p, assuming there are p species in the model.
#'
#' @export
#'
#' @seealso [coef.PLNfit()], [standard_error.PLNfit()] and [vcov.PLNfit()] for other ways to access \deqn{\Sigma}.
#'
#' @importFrom stats sigma
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPLN <- PLN(Abundance ~ 1 + offset(log(Offset)), data = trichoptera)
#' sigma(myPLN) ## Sigma
sigma.PLNfit <- function(object, ...) {
  stopifnot(isPLNfit(object))
  object$model_par$Sigma
}

#' Component-wise standard errors of Theta
#'
#' @description Extracts univariate standard errors for the estimated coefficient of Theta. Standard errors are computed from the (approximate) Fisher information matrix. See [fisher.PLNfit()] for more details on the approximations.
#'
#' @param object an R6 object with class PLNfit
#' @param type Either `Wald` (default) or `Louis`. Approximation scheme used to compute the Fisher information matrix
#'
#' @seealso [vcov.PLNfit()] for the complete Fisher information matrix
#'
#' @return A p * d positive matrix (same size as \eqn{\Theta}) with standard errors for the coefficients of \eqn{\Theta}
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPLN <- PLN(Abundance ~ 1 + offset(log(Offset)), data = trichoptera)
#' standard_error(myPLN, "wald")
#' @export
standard_error <- function(object, type) {
  UseMethod("standard_error", object)
}

#' @describeIn standard_error Component-wise standard errors of Theta in [`PLNfit`]
#' @export
standard_error.PLNfit <- function(object, type = c("wald", "louis")) {
  stopifnot(isPLNfit(object))
  type <- match.arg(type)
  if (type != object$fisher$type) {
    stop(paste("Standard errors were not computed using the", type, "approximation. Try another approximation scheme."))
  }
  object$std_err
}

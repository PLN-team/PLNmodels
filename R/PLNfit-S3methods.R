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

#' Predict counts conditionally
#'
#' @name predict_cond
#' @description Predict counts of a new sample conditionally on a (set of) observed variables
#' @param object an R6 object with class [`PLNfit`]
#' @param cond_responses a data frame containing the counts of the observed variables (matching the names provided as data in the PLN function)
#' @param newdata A data frame in which to look for variables and offsets with which to predict
#' @param var_par Boolean. Should new estimations of the variational parameters of mean and variance be sent back, as attributes of the matrix of predictions. Default to \code{FALSE}.
#' @param type The type of prediction required. The default is on the scale of the linear predictors (i.e. log average count)
#' @return A list containing:
#' \item{pred}{A matrix of predicted log-counts (if `type = "link"`) or predicted counts (if `type = "response"`)}
#' \item{M}{A matrix containing E(Z_uncond | Y_c) for each given site.}
#' \item{S}{A matrix containing Var(Z_uncond | Y_c) for each given site (sites are the third dimension of the array)}
#' @export
#' @examples
#' data(trichoptera)
#' trichoptera_prep <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPLN <- PLN(Abundance ~ Temperature + Wind, trichoptera_prep)
#' #Condition on the set of the first two species in the dataset (Hym, Hys) at the ten first sites
#' Yc <- trichoptera$Abundance[1:10, c(1, 2), drop=FALSE]
#' newX <- cbind(1, trichoptera$Covariate[1:10, c("Temperature", "Wind")])
#' pred <- predict_cond(myPLN, newX, Yc, type = "response")
predict_cond <- function(object, newdata, cond_responses, type = c("link", "response"), var_par = FALSE) {
  UseMethod("predict_cond", object)
}

#' @describeIn predict_cond Predict counts of a new sample conditionally on a (set of) observed variables for a [`PLNfit`]
#' @export
predict_cond.PLNfit = function(object, newdata, cond_responses, type = c("link", "response"), var_par = FALSE){
  stopifnot(isPLNfit(object))
  object$predict_cond(newdata, cond_responses, type, var_par, parent.frame())
}

#' Extract model coefficients
#'
#' @description Extracts model coefficients from objects returned by [PLN()] and its variants
#'
#' @name coef.PLNfit
#'
#' @param object an R6 object with class [`PLNfit`]
#' @param type type of parameter that should be extracted. Either "main" (default) for \deqn{B} or "covariance" for \deqn{\Sigma}
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
#' coef(myPLN) ## B
#' coef(myPLN, type = "covariance") ## Sigma
coef.PLNfit <- function(object, type = c("main", "covariance"), ...) {
  stopifnot(isPLNfit(object))
  switch(match.arg(type),
         main       = object$model_par$B,
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
#' @description Returns the variance-covariance matrix of the main parameters of a fitted [PLN()] model object. The main parameters of the model correspond to \deqn{B}, as returned by [coef.PLNfit()]. The function can also be used to return the variance-covariance matrix of the residuals. The latter matrix can also be accessed via [sigma.PLNfit()]
#'
#' @inheritParams coef.PLNfit
#' @return A matrix of variance/covariance extracted from the PLNfit model. If type="main" and \eqn{B} is a matrix of size d * p, the result is a block-diagonal matrix with p (number of species) blocks of size d (number of covariates). if type="main", it is a symmetric matrix of size p.
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
#' vcov(myPLN, type = "covariance") ## Sigma
vcov.PLNfit <- function(object, type = c("main", "covariance"), ...) {
  stopifnot(isPLNfit(object))
  type <- match.arg(type)
  if (type == "main" & is.null(attr(object$model_par$B, "vcov_variational")))
    stop("Variational estimation not available: rerun by setting `variational_var = TRUE` in the control list.")
  switch(type,
         main       = attr(object$model_par$B, "vcov_variational"),
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

#' Component-wise standard errors of B
#'
#' @description Extracts univariate standard errors for the estimated coefficient of B. Standard errors are computed from the (approximate) Fisher information matrix.
#'
#' @param object an R6 object with class PLNfit
#' @param type string describing the type of variance approximation: "variational", "jackknife", "sandwich" (only for fixed covariance). Default is "variational".
#' @param parameter string describing the target parameter: either B (regression coefficients) or Omega (inverse residual covariance)
#'
#' @seealso [vcov.PLNfit()] for the complete variance covariance estimation of the coefficient
#'
#' @return A p * d positive matrix (same size as \eqn{B}) with standard errors for the coefficients of \eqn{B}
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPLN <- PLN(Abundance ~ 1 + offset(log(Offset)), data = trichoptera,
#'               control = PLN_param(config_post = list(variational_var = TRUE)))
#' standard_error(myPLN)
#' @export
standard_error <- function(object, type = c("variational", "jackknife", "sandwich"), parameter = c("B", "Omega")) {
  UseMethod("standard_error", object)
}

#' @describeIn standard_error Component-wise standard errors of B in [`PLNfit`]
#' @export
standard_error.PLNfit <- function(object, type = c("variational", "jackknife", "bootstrap", "sandwich"), parameter = c("B", "Omega")) {
  type <- match.arg(type)
  par  <- match.arg(parameter)
  if (type == "variational" & is.null(attr(object$model_par$B, "variance_variational")))
    stop("Variational estimation not available: rerun by setting `variational_var = TRUE` in the control list.")
  if (type == "jackknife" & is.null(attr(object$model_par$B, "variance_jackknife")))
    stop("Jackknife estimation not available: rerun by setting `jackknife = TRUE` in the control list.")
  if (type == "bootstrap" & is.null(attr(object$model_par$B, "variance_bootstrap")))
    stop("Bootstrap estimation not available: rerun by setting `bootstrap > 0` in the control list.")
  if (type == "sandwich")
    stop("Sandwich estimator is only available for fixed covariance / precision matrix.")
  attr(object$model_par[[par]], paste0("variance_", type)) %>% sqrt()
}

#' @describeIn standard_error Component-wise standard errors of B in [`PLNfit_fixedcov`]
#' @export
standard_error.PLNfit_fixedcov <- function(object, type = c("variational", "jackknife", "bootstrap", "sandwich"), parameter = c("B", "Omega")) {
  type <- match.arg(type)
  par  <- match.arg(parameter)
  if (par == "Omega")
    stop("Omega is not estimated for fixed covariance model")
  if (type == "variational" & is.null(attr(object$model_par$B, "variance_variational")))
    stop("Variational estimation not available: rerun by setting `variational_var = TRUE` in the control list.")
  if (type == "jackknife" & is.null(attr(object$model_par$B, "variance_jackknife")))
    stop("Jackknife estimation not available: rerun by setting `jackknife = TRUE` in the control list.")
  if (type == "bootstrap" & is.null(attr(object$model_par$B, "variance_bootstrap")))
    stop("Bootstrap estimation not available: rerun by setting `bootstrap > 0` in the control list.")
  attr(object$model_par[[par]], paste0("variance_", type)) %>% sqrt()
}

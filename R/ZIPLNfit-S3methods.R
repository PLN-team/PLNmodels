## =========================================================================================
##
## PUBLIC S3 METHODS FOR ZIPLNfit
##
## =========================================================================================

## Auxiliary functions to check the given class of an objet
isZIPLNfit <- function(Robject) {inherits(Robject, "ZIPLNfit")}

#' Predict counts of a new sample
#'
#' @name predict.ZIPLNfit
#' @inheritParams predict.PLNfit
#' @param type Scale used for the prediction. Either `"link"` (default, predicted positions in the latent space), `"response"` (predicted average counts, accounting for zero-inflation) or `"deflated"` (predicted average counts, not accounting for zero-inflation and using only the PLN part of the model).
#'
#' @details
#' Note that `level = 1` can only be used if responses are provided,
#' as the variational parameters can't be estimated otherwise. In the absence of responses, `level` is ignored and the fitted values are returned
#'
#' Note also that when `type = "response"` corresponds to predicting
#' values with \eqn{(1 - \pi)A}, where \eqn{A} is the average count in
#' the PLN part of the model and \eqn{\pi} the probability of zero-inflation,
#' whereas `type = "deflated"` corresponds to \eqn{A}.
#'
#' @param object an R6 object with class [`ZIPLNfit`]
#' @return A matrix of predicted log-counts (if `type = "link"`) or predicted counts,
#' accounting for zero-inflation (if `type = "response"`) or not (if `type = "deflated"`).
#' @export
predict.ZIPLNfit <- function(object, newdata, responses = NULL, level = 1, type = c("link", "response", "deflated"), ...) {
  stopifnot(isZIPLNfit(object))
  object$predict(newdata = newdata, type = type, envir = parent.frame(), level = level, responses = responses)
}



#' Extract model coefficients
#'
#' @description Extracts model coefficients from objects returned by [ZIPLN()] and its variants
#'
#' @name coef.ZIPLNfit
#'
#' @param object an R6 object with class [`ZIPLNfit`]
#' @param type type of parameter that should be extracted. Either "count" (default) for \eqn{B},
#' "zero" for \eqn{B0}, "precision" for \eqn{\Omega}, "covariance" for \eqn{\Sigma}
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
coef.ZIPLNfit <- function(object, type = c("count", "zero", "precision", "covariance"), ...) {
  stopifnot(isZIPLNfit(object))
  switch(match.arg(type),
         count    = object$model_par$B,
         zero     = object$model_par$B0,
         precision  = object$model_par$Omega,
         covariance = object$model_par$Sigma)
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
#' @description Extract the variance-covariance matrix of the residuals, usually noted \eqn{\Sigma} in ZIPLN models.
#'
#' @inheritParams coef.ZIPLNfit
#'
#' @return A semi definite positive matrix of size p, assuming there are p species in the model.
#'
#' @export
#' @seealso [coef.ZIPLNfit()]
#'
#' @importFrom stats sigma
sigma.ZIPLNfit <- function(object, ...) {
  stopifnot(isZIPLNfit(object))
  object$model_par$Sigma
}

#' Extract log-likelihood of a fitted ZIPLN model
#'
#' @name logLik.ZIPLNfit
#' @description Returns the variational lower bound of the log-likelihood as a `"logLik"` object,
#' compatible with [stats::AIC()] and [stats::BIC()].
#'
#' @param object an R6 object with class [`ZIPLNfit`]
#' @param ... additional parameters for S3 compatibility. Not used
#' @return An object of class `"logLik"`. The numeric value is the variational ELBO.
#' Attributes `df` and `nobs` hold the number of parameters and observations.
#'
#' @importFrom stats logLik
#' @export
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' model <- ZIPLN(Abundance ~ 1, data = trichoptera)
#' logLik(model)
logLik.ZIPLNfit <- function(object, ...) {
  stopifnot(isZIPLNfit(object))
  structure(object$loglik, class = "logLik", df = object$nb_param, nobs = object$n)
}

#' Akaike Information Criterion for a fitted ZIPLN model
#'
#' @name AIC.ZIPLNfit
#' @description Computes the variational AIC as `loglik - nb_param` (larger is better).
#' This follows the maximization convention used throughout PLNmodels.
#'
#' @param object an R6 object with class [`ZIPLNfit`]
#' @param k not used, present for S3 compatibility.
#' @param ... additional parameters for S3 compatibility. Not used
#' @return A scalar: the variational AIC (larger is better).
#'
#' @importFrom stats AIC
#' @export
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' model <- ZIPLN(Abundance ~ 1, data = trichoptera)
#' AIC(model)
AIC.ZIPLNfit <- function(object, ..., k = 2) {
  stopifnot(isZIPLNfit(object))
  object$AIC
}

#' Bayesian Information Criterion for a fitted ZIPLN model
#'
#' @name BIC.ZIPLNfit
#' @description Computes the variational BIC as `loglik - 0.5 * log(n) * nb_param` (larger is better).
#' This follows the maximization convention used throughout PLNmodels.
#'
#' @param object an R6 object with class [`ZIPLNfit`]
#' @param ... additional parameters for S3 compatibility. Not used
#' @return A scalar: the variational BIC (larger is better).
#'
#' @importFrom stats BIC
#' @export
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' model <- ZIPLN(Abundance ~ 1, data = trichoptera)
#' BIC(model)
BIC.ZIPLNfit <- function(object, ...) {
  stopifnot(isZIPLNfit(object))
  object$BIC
}

#' @rdname ICL
#' @description `ICL.ZIPLNfit`: ICL for a fitted [`ZIPLNfit`].
#' @param object an R6 object with class [`ZIPLNfit`]
#' @export
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' model <- ZIPLN(Abundance ~ 1, data = trichoptera)
#' ICL(model)
ICL.ZIPLNfit <- function(object, ...) {
  stopifnot(isZIPLNfit(object))
  object$ICL
}

## =========================================================================================
##
## PUBLIC S3 METHODS FOR ZIPLNfit_sparse
##
## =========================================================================================

## Auxiliary functions to check the given class of an objet
isZIPLNfit_sparse <- function(Robject) {inherits(Robject, "ZIPLNfit_sparse")}

#' Extract and plot the network (partial correlation, support or inverse covariance) from a [`ZIPLNfit_sparse`] object
#'
#' @name plot.ZIPLNfit_sparse
#' @inheritParams plot.PLNnetworkfit
#' @param x an R6 object with class [`ZIPLNfit_sparse`]
#'
#' @inherit plot.PLNnetworkfit return
#'
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' fit <- ZIPLN(Abundance ~ 1, data = trichoptera, control = ZIPLN_param(penalty = 0.1))
#' \dontrun{
#' plot(fit)
#' }
#' @export
plot.ZIPLNfit_sparse <-
  function(x,
           type            = c("partial_cor", "support"),
           output          = c("igraph", "corrplot"),
           edge.color      = c("#F8766D", "#00BFC4"),
           remove.isolated = FALSE,
           node.labels     = NULL,
           layout          = layout_in_circle,
           plot            = TRUE, ...) {
    stopifnot(isZIPLNfit_sparse(x))
    invisible(x$plot_network(type, output, edge.color, remove.isolated, node.labels, layout, plot))
  }

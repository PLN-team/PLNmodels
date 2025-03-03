## =========================================================================================
##
## PUBLIC S3 METHODS FOR PLNmixturefit
##
## =========================================================================================

## Auxiliary functions to check the given class of an objet
isPLNmixturefit <- function(Robject) {inherits(Robject, "PLNmixturefit")}

#' Mixture visualization of a [`PLNmixturefit`] object
#'
#' Represent the result of the clustering either by coloring the individual in a two-dimension PCA factor map,
#' or by representing the expected matrix  of count reorder according to the clustering.
#'
#' @name plot.PLNmixturefit
#'
#' @param x an R6 object with class [`PLNmixturefit`]
#' @param type character for the type of plot, either "pca", for or "matrix". Default is `"pca"`.
#' @param main character. A title for the  plot. If NULL (the default), an hopefully appropriate title will be used.
#' @param plot logical. Should the plot be displayed or sent back as [`ggplot2::ggplot`] object
#' @param ... Not used (S3 compatibility).
#'
#' @return a [`ggplot2::ggplot`] graphic
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPLN <- PLNmixture(Abundance ~ 1 + offset(log(Offset)),
#'            data = trichoptera, control = PLNmixture_param(smoothing = "none"))  %>% getBestModel()
#' \dontrun{
#' plot(myPLN, "pca")
#' plot(myPLN, "matrix")
#' }
#' @export
plot.PLNmixturefit <-
  function(x,
           type           = c("pca", "matrix"),
           main           = NULL,
           plot           = TRUE, ...) {

    if (is.null(main))
      main <- switch(match.arg(type),
           "pca"    = "Clustering labels in Individual Factor Map",
           "matrix" = "Expected counts reorder by clustering")
    p <- switch(match.arg(type),
           "pca"    = x$plot_clustering_pca(main = main, plot = FALSE),
           "matrix" = x$plot_clustering_data(main = main, plot = FALSE))
    if (plot) print(p)
    invisible(p)
  }


#' Prediction for a [`PLNmixturefit`] object
#'
#' Predict either posterior probabilities for each group or latent positions based on new samples
#'
#' @param object an R6 object with class [`PLNmixturefit`]
#' @param newdata A data frame in which to look for variables, offsets and counts with which to predict.
#' @param type The type of prediction required. The default `posterior` are posterior probabilities for each group ,
#'  `response` is the group with maximal posterior probability and `latent` is the averaged latent in the latent space,
#'  with weights equal to the posterior probabilities.
#' @param prior User-specified prior group probabilities in the new data. The default uses a uniform prior.
#' @param control a list-like structure for controlling the fit. See [PLNmixture_param()] for details.
#' @param ... additional parameters for S3 compatibility. Not used
#' @return A matrix of posterior probabilities for each group (if type = "posterior"), a matrix of (average) position in the
#' latent space (if type = "position") or a vector of predicted groups (if type = "response").
#' @export
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPLN <- PLNmixture(Abundance ~ 1 + offset(log(Offset)),
#'            data = trichoptera, control = PLNmixture_param(smoothing = "none"))  %>% getBestModel()
#' predict(myPLN, trichoptera, "posterior")
#' predict(myPLN, trichoptera, "position")
#' predict(myPLN, trichoptera, "response")
predict.PLNmixturefit <-
  function(object, newdata,
           type = c("posterior", "response", "position"),
           prior = matrix(rep(1/object$k, object$k), nrow(newdata), object$k, byrow = TRUE),
           control = PLNmixture_param(), ...) {

  stopifnot(isPLNmixturefit(object))
  object$predict(newdata, type, prior, control, parent.frame())

}

#' Extract model coefficients
#'
#' @description Extracts model coefficients from objects returned by [PLN()] and its variants
#'
#' @name coef.PLNmixturefit
#'
#' @param object an R6 object with class [`PLNmixturefit`]
#' @param type type of parameter that should be extracted. Either "main" (default) for \deqn{\Theta},
#' "means" for \deqn{\mu}, "mixture" for \deqn{\pi} or "covariance" for \deqn{\Sigma}
#' @param ... additional parameters for S3 compatibility. Not used
#' @return A matrix of coefficients extracted from the PLNfit model.
#'
#' @seealso [sigma.PLNmixturefit()]
#'
#' @export
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPLN <- PLNmixture(Abundance ~ 1 + offset(log(Offset)),
#'            data = trichoptera, control = PLNmixture_param(smoothing = "none"))  %>% getBestModel()
#' coef(myPLN) ## Theta - empty here
#' coef(myPLN, type = "mixture") ## pi
#' coef(myPLN, type = "means") ## mu
#' coef(myPLN, type = "covariance") ## Sigma
coef.PLNmixturefit <- function(object, type = c("main", "means", "covariance", "mixture"), ...) {
  stopifnot(isPLNmixturefit(object))
  switch(match.arg(type),
         main       = object$model_par$Theta,
         means      = object$group_means,
         mixture    = object$mixtureParam,
         covariance = object$model_par$Sigma)
}

#' Extracts model fitted values from objects returned by [PLNmixture()] and its variants
#'
#' @name fitted.PLNmixturefit
#'
#' @inheritParams coef.PLNmixturefit
#' @return A matrix of Fitted values extracted from the object object.
#'
#' @export
fitted.PLNmixturefit <- function(object, ...) {
  stopifnot(isPLNmixturefit(object))
  object$fitted
}

#' Extract variance-covariance of residuals 'Sigma'
#'
#' @name sigma.PLNmixturefit
#' @description Extract the variance-covariance matrix of the residuals, usually noted \deqn{\Sigma} in PLN models. This captures the correlation between the species in the latent space.
#' or PLNmixture, it is a weighted mean of the variance-covariance matrices of each component.
#'
#' @inheritParams coef.PLNmixturefit
#'
#' @return A semi definite positive matrix of size p, assuming there are p species in the model.
#'
#' @export
#'
#' @seealso [coef.PLNmixturefit()] for other ways to access \deqn{\Sigma}.
#'
#' @importFrom stats sigma
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPLN <- PLNmixture(Abundance ~ 1 + offset(log(Offset)),
#'            data = trichoptera, control = PLNmixture_param(smoothing = "none"))  %>% getBestModel()
#' sigma(myPLN) ## Sigma
sigma.PLNmixturefit <- function(object, ...) {
  stopifnot(isPLNmixturefit(object))
  object$model_par$Sigma
}

#' @describeIn standard_error Component-wise standard errors of B in [`PLNmixturefit`] (not implemented yet)
#' @export
standard_error.PLNmixturefit <- function(object, type = c("variational", "jackknife", "sandwich"), parameter = c("B", "Omega")) {
  par  <- match.arg(parameter)
  if (par == "Omega")
    stop("Omega is not estimated as such in PLNmixture models")
  if (par == "B") {
    warning("Standard error of B is not implemented yet for PLNmixture models")
    theta_sd <- coef.PLNmixturefit(object)
    theta_sd[ , ] <- NA
    theta_sd
  }
}

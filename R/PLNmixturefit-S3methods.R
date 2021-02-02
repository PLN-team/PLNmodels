## =========================================================================================
##
## PUBLIC S3 METHODS FOR PLNmixturefit
##
## =========================================================================================

## Auxiliary functions to check the given class of an objet
isPLNmixturefit <- function(Robject) {inherits(Robject, "PLNmixturefit")}

#' Mixture vizualisation of a [`PLNmixturefit`] object
#'
#' Represent the result of the clustering either by coloring the individual in a two-dimension PCA factor map,
#' or by representing the expected matrix  of count reorder according to the clustering.
#'
#' @name plot.PLNmixturefit
#'
#' @param x an R6 object with class [`PLNmixturefit`]
#' @param type character for the type of plot, either "pca", for or "matrix". Default is `"pca"`.
#' @param main character. A title for the  plot. If NULL (the default), an hopefully appropriate title will be used.
#' @param plot logical. Should the plot be displayed or sent back as [`ggplot`] object
#' @param ... Not used (S3 compatibility).
#'
#' @return a [`ggplot`] graphic
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' fits <- PLNmixture(Abundance ~ 1, clusters = 1:4, data = trichoptera)
#' myMixture <- getBestModel(fits)
#' \dontrun{
#' plot(myMixture, "pca")
#' plot(myMixture, "matrix")
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
#' @param newdata A data frame in which to look for variables, offsets and counts with which to predict.
#' @param type The type of prediction required. The default `posterior` are posterior probabilities for each group ,
#'  `response` is the group with maximal posterior probability and `latent` is the averaged latent in the latent space,
#'  with weights equal to the posterior probabilities.
#' @param prior User-specified prior group probabilities in the new data. The default uses a uniform prior.
#' @param control a list for controlling the optimization. See [PLN()] for details.
#' @param ... additional parameters for S3 compatibility. Not used
#' @return A matrix of posterior probabilities for each group (if type = "posterior"), a matrix of (average) position in the
#' latent space (if type = "position") or a vector of predicted groups (if type = "response").
#' @export
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' fits <- PLNmixture(Abundance ~ 1, clusters = 1:4, data = trichoptera)
#' myMixture <- getBestModel(fits)
#' predict(myMixture, "posterior")
#' predict(myMixture, "latent")
#' predict(myMixture, "response")
predict.PLNmixturefit <-
  function(object, newdata,
           type = c("posterior", "response", "position"),
           prior = matrix(rep(1/object$k, object$k), nrow(newdata), object$k, byrow = TRUE),
           control = list(), ...) {

  stopifnot(isPLNmixturefit(object))
  object$predict(newdata, type, prior, control, parent.frame())

}

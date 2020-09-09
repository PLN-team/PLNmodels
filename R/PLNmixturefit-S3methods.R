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
#' fits <- PLNmixture(Abundance ~ 1, data = trichoptera)
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

## =========================================================================================
##
## PUBLIC S3 METHODS FOR PLNLDAfit
##
## =========================================================================================

## Auxiliary functions to check the given class of an objet
isPLNLDAfit <- function(Robject) {inherits(Robject, "PLNLDAfit"       )}

#' LDA visualization (individual and/or variable factor map(s)) for a [`PLNPCAfit`] object
#'
#' @name plot.PLNLDAfit
#'
#' @param x an R6 object with class PLNPCAfit
#' @param map the type of output for the PCA visualization: either "individual", "variable" or "both". Default is "both".
#' @param nb_axes scalar: the number of axes to be considered when map = "both". The default is min(3,rank).
#' @param axes numeric, the axes to use for the plot when map = "individual" or "variable". Default it c(1,min(rank))
#' @param var_cols a character or factor to define the color associated with the variables. By default, all variables receive the default color of the current palette.
#' @param plot logical. Should the plot be displayed or sent back as [`ggplot2`] object
#' @param main character. A title for the single plot (individual or variable factor map). If NULL (the default), an hopefully appropriate title will be used.
#' @param ... Not used (S3 compatibility).
#'
#' @return displays an individual and/or variable factor maps for the corresponding axes, and/or sends back a [`ggplot2`] or gtable object
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPLNLDA <- PLNLDA(Abundance ~ 1, grouping = Group, data = trichoptera)
#' \dontrun{
#' plot(myPLNLDA, map = "individual", nb_axes = 2)
#' }
#' @export
plot.PLNLDAfit <-
  function(x,
           map      = c("both", "individual", "variable"),
           nb_axes  = min(3, x$rank),
           axes     = seq.int(min(2,x$rank)),
           var_cols = "var_colors",
           plot     = TRUE,
           main = NULL,
           ...) {

    stopifnot(isPLNLDAfit(x))

    map <- match.arg(map)

    if (is.null(main)) {
      if (map == "individual")
        main <- "Individual factor map"
      if (map == "variable")
        main <- "Variable factor map"
    }

    if (map == "individual")
      p <- x$plot_individual_map(main = main, plot = plot)
    if (map == "variable")
      p <- x$plot_correlation_map(cols = var_cols, main = main, plot = plot)
    if (map == "both")
      p <- x$plot_LDA(nb_axes = nb_axes, var_cols = var_cols, plot = plot)

    invisible(p)
}

#' Predict group of new samples
#'
#' @name predict.PLNLDAfit
#'
#' @param object an R6 object with class [`PLNLDAfit`]
#' @param newdata A data frame in which to look for variables, offsets and counts  with which to predict.
#' @param type The type of prediction required. The default are posterior probabilities for each group (in either unnormalized log-scale or natural probabilities, see "scale" for details), "response" is the group with maximal posterior probability and "scores" is the average score along each separation axis in the latent space, with weights equal to the posterior probabilities.
#' @param scale The scale used for the posterior probability. Either log-scale ("log", default) or natural probabilities summing up to 1 ("prob").
#' @param prior User-specified prior group probabilities in the new data. If NULL (default), prior probabilities are computed from the learning set.
#' @param control a list for controlling the optimization. See [PLN()] for details.
#' @param ... additional parameters for S3 compatibility. Not used
#' @return A matrix of posterior probabilities for each group (if type = "posterior"), a matrix of (average) scores in the latent space (if type = "scores") or a vector of predicted groups (if type = "response").
#' @export
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myLDA <- PLNLDA(Abundance ~ 0 + offset(log(Offset)),
#'                 grouping = Group,
#'                 data = trichoptera)
#' \dontrun{
#' post_probs <- predict(myLDA, newdata = trichoptera, type = "posterior", scale = "prob")
#' head(round(post_probs, digits = 3))
#' predicted_group <- predict(myLDA, newdata = trichoptera, type = "response")
#' table(predicted_group, trichoptera$Group, dnn = c("predicted", "true"))
#' }
predict.PLNLDAfit <- function(object, newdata,
                              type = c("posterior", "response", "scores"),
                              scale = c("log", "prob"),
                              prior = NULL,
                              control = list(), ...) {
  stopifnot(isPLNLDAfit(object))
  object$predict(newdata, type, scale, prior, control, parent.frame())
}

#' Extracts model coefficients from objects returned by [PLNLDA()]
#'
#' @description The method for objects returned by [PLNLDA()] only returns
#'              coefficients associated to the \deqn{\Theta} part of the model (see the PLNLDA vignette
#'              for mathematical details).
#'
#' @name coef.PLNLDAfit
#'
#' @param object an R6 object with class PLNfit
#' @param ... additional parameters for S3 compatibility. Not used
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPLNLDA <- PLNLDA(Abundance ~ Wind, grouping = Group, data = trichoptera)
#' coef(myPLNLDA)
#' @return Either NULL or a matrix of coefficients extracted from the PLNLDAfit model.
#'
#' @export
coef.PLNLDAfit <- function(object, ...) {
  stopifnot(isPLNLDAfit(object))
  n.covariates <- object$d - object$rank - 1
  if (n.covariates == 0) return(NULL)
  object$model_par$Theta[ , 1:n.covariates, drop = FALSE]
}


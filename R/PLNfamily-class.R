#' An R6 Class to represent a collection of PLNfit
#'
#' @description super class for [`PLNPCAfamily`] and [`PLNnetworkfamily`]
#'
#'
#'
#' @details This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for \code{\link[=getBestModel.PLNnetworkfamily]{getBestModel()}},
#' \code{\link[=getModel.PLNnetworkfamily]{getModel()}} and  \code{\link[=plot.PLNnetworkfamily]{plot()}}.
#'
#' @md
#' @include PLNfamily-class.R
#' @importFrom R6 R6Class
PLNfamily <-
  R6Class(classname = "PLNfamily",
    public = list(
      #' @field responses the matrix of responses common to every models
      responses  = NULL, # the Y matrix
      #' @field covariates the matrix of covariates common to every models
      covariates = NULL, # the X matrix
      #' @field offsets the matrix of offsets common to every models
      offsets    = NULL, # the O matrix
      #' @field weights the vector of observation weights
      weights    = NULL, # the vector of obervation weights
      #' @field inception a [PLNfit] object, obtained when no sparsifying penalty is applied.
      inception  = NULL, # the basic model in the collection (no regularization, nor sparsity, nor rank)
      #' @field models a list of [PLNfit] object, one per penalty.
      models     = NULL  # the collection of fitted models
    ),
    private = list(
      params     = NULL, # vector of parameters that indexes the models (either sparsity, rank, etc.)
      n          = NULL, # number of samples
      p          = NULL, # number of responses
      d          = NULL  # number of covariates
    ),
    active = list(
      #' @field criteria a data frame with the values of some criteria (variational lower bound J, BIC, ICL and R2) for the collection of models / fits
      criteria = function() {
        res <- do.call(rbind, lapply(self$models, function(model) {model$criteria}))
        data.frame(param = private$params, res)
      },
      #' @field convergence sends back a data frame with some convergence diagnostics associated with the optimization process (method, optimal value, etc)
      convergence = function() {
        res <- do.call(rbind, lapply(self$models, function(model) {
          c(nb_param = model$nb_param, sapply(model$optim_par, function(x) x[length(x)]))
        }))
        data.frame(param = private$params, res, stringsAsFactors = FALSE)
      }
    )
)

#' @description
#' Create a new [`PLNfamily`] object.
#' @param responses the matrix of responses common to every models
#' @param covariates the matrix of covariates common to every models
#' @param offsets the matrix of offsets common to every models
#' @param weights the vector of observation weights
#' @param control a list for controlling the optimization. See details.
#'
#' @inherit PLN details
#' @md
#'
#' @return A new [`PLNfamily`] object
PLNfamily$set("public", "initialize",
  function(responses, covariates, offsets, weights, control) {

    ## set data matrice and dimension
    self$responses  <- responses
    self$covariates <- covariates
    self$offsets    <- offsets
    self$weights    <- weights
    private$n <- nrow(responses)
    private$p <- ncol(responses)
    private$d <- ncol(covariates)

    ## set names of the data matrices
    if (is.null(rownames(responses)))  rownames(self$responses)  <- 1:private$n
    if (is.null(colnames(responses)))  colnames(self$responses)  <- 1:private$p
    if (is.null(rownames(covariates))) rownames(self$covariates) <- 1:private$n
    if (is.null(colnames(covariates)) & private$d > 0) colnames(self$covariates) <- 1:private$d

})

## a method to compute and set fields after optimization
#' @description
#' Update and set some fields (`R^2`, `vcov`, etc) after optimization
#' @md
#' @return An updated [`PLNfamily`] object
PLNfamily$set("public", "postTreatment",
function() {
  nullModel <- nullModelPoisson(self$responses, self$covariates, self$offsets, self$weights)
  for (model in self$models)
    model$postTreatment(
      self$responses,
      self$covariates,
      self$offsets,
      self$weights,
      nullModel = nullModel
    )
})

## a method to compute and set fields after optimization
#' @description
#' Extract a model from a collection of models
#' @inheritParams getModel
#' @md
#' @return A [`PLNfit`] object (potentially a [`PLNPCAfit`] or [`PLNnetworkfit`])
#' @seealso [getModel()]
PLNfamily$set("public", "getModel",
function(var, index = NULL) {
  ## Extraction by index
  if (!is.null(index) && index <= length(self$models)) {
    return(self$models[[index[1]]]$clone())
  }
  ## Extraction by parameter value
  id <- match(var, private$params)
  if (!is.na(id)) { ## Exact match
    return(self$models[[id]]$clone())
  } else { ## No exact match
    id <- which.min(abs(var - private$params)) ## closest model (in terms of parameter value)
    warning(paste("No such a model in the collection. Acceptable parameter values can be found via",
               "$ranks() (for PCA)",
               "$penalties() (for network)",
               paste("Returning model with closest value. Requested:", var, ", returned:", private$params[id]),
               sep = "\n"))
    return(self$models[[id]]$clone())
  }
})


## a method to plot a collection of models
#' @description
#' Lineplot of selected criteria for all models in the collection
#' @param criteria A valid model selection criteria for the collection of models. Includes loglik, BIC (all), ICL (PLNPCA) and pen_loglik, EBIC (PLNnetwork)
#' @param annotate Logical. Should R2 be added to the plot (defaults to `TRUE`)
#' @md
#' @return A [`ggplot2`] object
PLNfamily$set("public", "plot",
function(criteria, annotate = TRUE) {
  stopifnot(!anyNA(self$criteria[criteria]))
  dplot <- self$criteria %>%
    dplyr::select(c("param", criteria)) %>%
    tidyr::gather(key = "criterion", value = "value", -param) %>%
    dplyr::group_by(criterion)
  p <- ggplot(dplot, aes(x = param, y = value, group = criterion, colour = criterion)) +
       geom_line() + geom_point() + ggtitle("Model selection criteria") + theme_bw()

  if (annotate)
    p <- p + annotate("text", x = self$criteria$param, y = min(dplot$value), hjust=-.1, angle = 90,
           label = paste("R2 =", round(self$criteria$R_squared, 2)), size = 3, alpha = 0.7)
  p
})

## Print method
#' @description
#' Base print method for collections of models.
#' @md
PLNfamily$set("public", "show",
function() {
  cat("--------------------------------------------------------\n")
  cat("COLLECTION OF", length(self$models), "POISSON LOGNORMAL MODELS\n")
  cat("--------------------------------------------------------\n")
})

## Print method
#' @description
#' Base print method for collections of models.
#' @md
PLNfamily$set("public", "print", function() self$show())


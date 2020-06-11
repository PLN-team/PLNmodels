#' An R6 Class to represent a collection of PLNfit
#'
#' @description super class for [`PLNPCAfamily`] and [`PLNnetworkfamily`]. The R6 class benefits from S3 methods such as [getBestModel()], [getModel()] and [`plot()`][plot.PLNPCAfamily()].
#'
#' @param responses the matrix of responses common to every models
#' @param covariates the matrix of covariates common to every models
#' @param offsets the matrix of offsets common to every models
#' @param weights the vector of observation weights
#' @param control a list for controlling the optimization. See details.
#' @param var value of the parameter (`rank` for PLNPCA, `sparsity` for PLNnetwork) that identifies the model to be extracted from the collection. If no exact match is found, the model with closest parameter value is returned with a warning.
#' @param index Integer index of the model to be returned. Only the first value is taken into account.
#'
#' @inherit PLN details
#'
#' @seealso
#'
#' @rdname PLNfamily
#' @include PLNfamily-class.R
#' @importFrom R6 R6Class
PLNfamily <-
  R6Class(
    classname = "PLNfamily",
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## PUBLIC MEMBERS ------
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    public = list(
      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Public fields ---------------------
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
      models     = NULL,  # the collection of fitted models

      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Creator functions -----------------
      #' @description Create a new [`PLNfamily`] object.
      #' @return A new [`PLNfamily`] object
      initialize = function(responses, covariates, offsets, weights, control) {

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

      },

      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Post treatment --------------------
      #' @description Update fields after optimization
      postTreatment = function() {
        nullModel <- nullModelPoisson(self$responses, self$covariates, self$offsets, self$weights)
        for (model in self$models)
          model$postTreatment(
            self$responses,
            self$covariates,
            self$offsets,
            self$weights,
            nullModel = nullModel
          )
      },

      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Extractors ------------------------
      #' @description
      #' Extract a model from a collection of models
      # @inheritParams getModel
      #' @return A [`PLNfit`] object
      #' @seealso [getModel()]
      getModel = function(var, index = NULL) {
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
      },

      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Graphical methods -----------------
      #' @description
      #' Lineplot of selected criteria for all models in the collection
      #' @param criteria A valid model selection criteria for the collection of models. Includes loglik, BIC (all), ICL (PLNPCA) and pen_loglik, EBIC (PLNnetwork)
      #' @param annotate Logical. Should R2 be added to the plot (defaults to `TRUE`)
      #' @return A [`ggplot2`] object
      plot = function(criteria, annotate = TRUE) {
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
      },

      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Print methods ---------------------
      #' @description User friendly print method
      show = function() {
        cat("--------------------------------------------------------\n")
        cat("COLLECTION OF", length(self$models), "POISSON LOGNORMAL MODELS\n")
        cat("--------------------------------------------------------\n")
      },

      #' @description User friendly print method
      print = function() { self$show() }

      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Other public members --------------


    ),

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## PRIVATE MEMBERS ----
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    private = list(
      params     = NULL, # vector of parameters that indexes the models (either sparsity, rank, etc.)
      n          = NULL, # number of samples
      p          = NULL, # number of responses
      d          = NULL  # number of covariates
    ),

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## ACTIVE BINDINGS ----
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ##  END OF CLASS ----
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  )



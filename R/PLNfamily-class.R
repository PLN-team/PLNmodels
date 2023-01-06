#' An R6 Class to represent a collection of PLNfit
#'
#' @description super class for [`PLNPCAfamily`] and [`PLNnetworkfamily`].
#'
#' @param responses the matrix of responses common to every models
#' @param covariates the matrix of covariates common to every models
#' @param offsets the matrix of offsets common to every models
#' @param weights the vector of observation weights
#' @param var value of the parameter (`rank` for PLNPCA, `sparsity` for PLNnetwork) that identifies the model to be extracted from the collection. If no exact match is found, the model with closest parameter value is returned with a warning.
#' @param control list controlling the optimization and the model
#' @param index Integer index of the model to be returned. Only the first value is taken into account.
#'
#' @inherit PLN details
#'
#' @rdname PLNfamily
#' @include PLNfamily-class.R
#' @importFrom R6 R6Class
#' @importFrom purrr map reduce map_chr
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
      weights    = NULL, # the vector of observation weights
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
      #' @param config a list for controlling the post-treatment.
      postTreatment = function(config) {
        nullModel <- nullModelPoisson(self$responses, self$covariates, self$offsets, self$weights)
        for (model in self$models)
          model$postTreatment(
            self$responses,
            self$covariates,
            self$offsets,
            self$weights,
            config,
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
                        "$ranks (for PCA)",
                        "$clusters (for mixture models)",
                        "$penalties (for network)",
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
      #' @param reverse A logical indicating whether to plot the value of the criteria in the "natural" direction
      #' (loglik - penalty) or in the "reverse" direction (-2 loglik + penalty). Default to FALSE, i.e use the natural direction, on
      #' the same scale as the log-likelihood.
      #' @return A [`ggplot2`] object
      plot = function(criteria, reverse) {
        stopifnot(!anyNA(self$criteria[criteria]))
        dplot <- self$criteria %>%
          dplyr::select(dplyr::all_of(c("param", criteria))) %>%
          tidyr::gather(key = "criterion", value = "value", -param) %>%
          {if (reverse)
            dplyr::mutate(
              .data = .,
              value = -2 * value,
              criterion = dplyr::if_else(grepl('loglik', criterion), paste0("-2", criterion), criterion))
            else . } %>%
          dplyr::group_by(criterion)
        p <- ggplot(dplot, aes(x = param, y = value, group = criterion, colour = criterion)) +
          geom_line() + geom_point() +
          ggtitle(label    = "Model selection criteria",
                  subtitle = if (reverse) "Lower is better" else "Higher is better") +
          theme_bw()
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
      params = NULL, # vector of parameters that indexes the models (either sparsity, rank, number of cluster, etc.)
      n      = NULL, # number of samples
      p      = NULL, # number of responses
      d      = NULL  # number of covariates
    ),

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## ACTIVE BINDINGS ----
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    active = list(
      #' @field criteria a data frame with the values of some criteria (approximated log-likelihood, BIC, ICL, etc.) for the collection of models / fits
      #' BIC and ICL are defined so that they are on the same scale as the model log-likelihood, i.e. with the form, loglik - 0.5 penalty
      criteria = function() {
        res <- purrr::map(self$models, 'criteria') %>%
          purrr::reduce(rbind)
        data.frame(param = private$params, res)
      },
      #' @field convergence sends back a data frame with some convergence diagnostics associated with the optimization process (method, optimal value, etc)
      convergence = function() {
        res <- purrr::map(self$models, function(model) {
          c(nb_param = model$nb_param, purrr::map_chr(model$optim_par, tail, 1))
        }) %>% purrr::reduce(rbind)
        data.frame(param = private$params, res, stringsAsFactors = FALSE)
      }
    )

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ##  END OF CLASS ----
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  )


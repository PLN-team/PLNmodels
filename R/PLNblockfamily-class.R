#' An R6 Class to represent a collection of PLNblockfit
#'
#' @description The function [PLNblock()] produces an instance of this class.
#'
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for [getBestModel()],
#' [getModel()] and [plot()][plot.PLNnetworkfamily()]
#'
## Parameters shared by many methods
#' @param nb_blocks a vector of positive real number controlling the level of sparsity of the underlying network.
#' @param sparsity tuning parameter for controlling the sparsity level of Omega/Sigma
#' @param responses the matrix of responses common to every models
#' @param covariates the matrix of covariates common to every models
#' @param offsets the matrix of offsets common to every models
#' @param weights the vector of observation weights
#' @param formula model formula used for fitting, extracted from the formula in the upper-level call
#' @param control a list for controlling the optimization.
#' @param var value of the parameter (`rank` for PLNPCA, `sparsity` for PLNnetwork, `nb_bloc` for PLNblock) that identifies the model to be extracted from the collection. If no exact match is found, the model with closest parameter value is returned with a warning.
#' @param index Integer index of the model to be returned. Only the first value is taken into account
#'
#' @include PLNfamily-class.R
#' @importFrom R6 R6Class
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' fits <- PLNblock(Abundance ~ 1, nb_blocks = 1:5, data = trichoptera)
#' class(fits)
#' @seealso The function [PLNblock()], the class [`PLNblockfit`]
PLNblockfamily <- R6Class(
  classname = "PLNblockfamily",
  inherit = PLNfamily,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ------
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Creation functions ----------------
    #' @description Initialize all models in the collection
    #' @return Update current [`PLNblockfit`] with smart starting values
    initialize = function(nb_blocks, sparsity, responses, covariates, offsets, weights, formula, control) {

      ## Initialize fields shared by the super class
      super$initialize(responses, covariates, offsets, weights, control)
      private$params  <- nb_blocks

      ## Default clustering is obtained by CAH on the variational parameters of the means of a fully parametrized PLN
      control_init <- control
      control_init$config_optim <- config_default_nlopt
      control_init$backend <- "nlopt"

      if (is.numeric(control$init_cl)) {
        blocks <- control$init_cl
      }else{
        myPLN <- PLNfit$new(responses, covariates, offsets, rep(1, nrow(responses)), formula, control_init)
        myPLN$optimize(responses, covariates, offsets, weights, control_init$config_optim)
        if(control$init_cl=="kmeans"){
          Means <- t(myPLN$var_par$M)
          blocks <- lapply(seq(1:nb_blocks), function(k) kmeans(Means, centers=k)$clusters)
        }else{
          D <- 1 - abs(cov2cor(myPLN$model_par$Sigma))
          ## D <- diag(diag(myPLN$model_par$Sigma)) - abs(cov(myPLN$model_par$Sigma))
          blocks <- hclust(as.dist(D), method = "ward.D2") %>% cutree(nb_blocks) %>% as.data.frame() %>% as.list()
          # blocks <- lapply(nb_blocks, function(k) kmeans(D, centers = k, nstart = 30)$cl)
        }
      }

      ## instantiate as many models as cluterings
      self$models <-
      blocks %>%
        map(as_indicator) %>%
        map(.check_boundaries) %>%
        map(function(blocks_) {
          if (sparsity > 0) {
            model <- PLNblockfit_sparse$new(blocks_, sparsity, responses, covariates, offsets, weights, formula, control)
          } else {
            model <- PLNblockfit$new(blocks_, responses, covariates, offsets, weights, formula, control)
          }
        })
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Optimization ----------------------
    #' @description Call to the C++ optimizer on all models of the collection
    #' @param config a list for controlling the optimization.
    optimize = function(config) {
      self$models <- future.apply::future_lapply(self$models, function(model) {
        if (config$trace >= 1) {
          cat("\tnumber of blocks =", model$nb_block, "\r")
          flush.console()
        }
        model$optimize(self$responses, self$covariates, self$offsets, self$weights, config)
        model
      }, future.seed = TRUE, future.scheduling = structure(TRUE, ordering = "random"))
    },

    #' @description Extract best model in the collection
    #' @param crit a character for the criterion used to performed the selection. Either
    #' "BIC", "ICL" or "loglik". Default is `ICL`
    #' @return a [`PLNblockfit`] object
    getBestModel = function(crit = c("BIC", "ICL", "loglik")){
      crit <- match.arg(crit)
      stopifnot(!anyNA(self$criteria[[crit]]))
      id <- 1
      if (length(self$criteria[[crit]]) > 1) {
        id <- which.max(self$criteria[[crit]])
      }
      model <- self$models[[id]]$clone()
      model
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Graphical methods -----------------
    #' @description Display various outputs (goodness-of-fit criteria, robustness, diagnostic) associated with a collection of PLNblocks fits (a [`PLNblockfamily`])
    #' @param criteria vector of characters. The criteria to plot in `c("loglik", "BIC", "ICL")`. Defaults to all of them.
    #' @param reverse A logical indicating whether to plot the value of the criteria in the "natural" direction
    #' (loglik - 0.5 penalty) or in the "reverse" direction (-2 loglik + penalty). Default to FALSE, i.e use the
    #' natural direction, on the same scale as the log-likelihood..
    #' @return a [`ggplot`] graph
    plot = function(criteria = c("loglik", "ICL", "BIC"), reverse = FALSE) {
      vlines <- sapply(intersect(criteria, c("BIC", "ICL")) , function(crit) self$getBestModel(crit)$nb_block)
      p <- super$plot(criteria, reverse) + xlab("# blocks") + geom_vline(xintercept = vlines, linetype = "dashed", alpha = 0.25)
      p
    },

    #' @description Plot objective value of the optimization problem along the penalty path
    #' @return a [`ggplot`] graph
    plot_objective = function() {
      objective <- unlist(lapply(self$models, function(model) model$optim_par$objective))
      changes <- cumsum(unlist(lapply(self$models, function(model) model$optim_par$outer_iterations)))
      dplot <- data.frame(iteration = 1:length(objective), objective = objective)
      p <- ggplot(dplot, aes(x = iteration, y = objective)) + geom_line() +
        geom_vline(xintercept = changes, linetype="dashed", alpha = 0.25) +
        ggtitle("Objective along the alternate algorithm") + xlab("iteration (+ changes of model)") +
        annotate("text", x = changes, y = min(dplot$objective), angle = 90,
                 label = paste("nb blocks=",format(self$criteria$param, digits = 1)), hjust = -.1, size = 3, alpha = 0.7) + theme_bw()
      p
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Print methods ---------------------
    #' @description User friendly print method
    show = function() {
      super$show()
      cat(" Task: Blockwise Covariance \n")
      cat("========================================================\n")
      cat(" -", length(self$nb_blocks) , "blocks considered: from", min(self$nb_blocks), "to", max(self$nb_blocks), "\n")
      cat(" - Best model (greater BIC): nb blocks =", format(self$getBestModel("BIC")$nb_block, digits = 3), "\n")
      cat(" - Best model (greater ICL): nb blocks =", format(self$getBestModel("ICL")$nb_block, digits = 3), "\n")
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## ACTIVE BINDINGS ------
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field nb_blocks vector of number of blocks for regrouping the variables/responses (dimension of the residual covariance)
    nb_blocks = function() private$params
  )

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## END OF CLASS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
)


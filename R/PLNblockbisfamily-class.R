# library(Rmixmod)
#' An R6 Class to represent a collection of PLNblockbisfit
#'
#' @description The function [PLNblockbis()] produces an instance of this class.
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
#' @param var value of the parameter (`rank` for PLNPCA, `sparsity` for PLNnetwork, `nb_bloc` for PLNblockbis) that identifies the model to be extracted from the collection. If no exact match is found, the model with closest parameter value is returned with a warning.
#' @param index Integer index of the model to be returned. Only the first value is taken into account
#'
#' @include PLNfamily-class.R
#' @importFrom R6 R6Class
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' fits <- PLNblockbis(Abundance ~ 1, nb_blocks = 1:17, data = trichoptera)
#' class(fits)
#' @seealso The function [PLNblockbis()], the class [`PLNblockbisfit`]
PLNblockbisfamily <- R6Class(
  classname = "PLNblockbisfamily",
  inherit = PLNfamily,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ------
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Creation functions ----------------
    #' @description Initialize all models in the collection
    #' @import ClustOfVar
    #' @return Update current [`PLNblockbisfit`] with smart starting values
    initialize = function(nb_blocks, sparsity, responses, covariates, offsets, weights, formula, control) {

      ## Initialize fields shared by the super class
      super$initialize(responses, covariates, offsets, weights, control)
      private$params  <- nb_blocks
      if (length(sparsity) == 1) sparsity <- rep(sparsity, length(nb_blocks))
      stopifnot(all.equal(length(sparsity),length(nb_blocks)))

      ## ==================================================
      ##
      ## Common PLN (use as starting point for B, M, S, Mu, Delta)
      if (isPLNfit(control$inception)) {
        if (control$trace > 1) cat("\n User defined inceptive PLN model")
        myPLN_init <- control$inception
      } else {
        control_init <- control
        control_init$config_optim <- config_default_nlopt
        control_init$backend <- "nlopt"
        control_init$covariance <- "diagonal"
        myPLN_init <- PLNfit$new(responses, covariates, offsets, weights, formula, control_init)
        myPLN_init$optimize(responses, covariates, offsets, weights, control_init$config_optim)
        control$inception <- myPLN_init
      }

      ## ==================================================
      ## Initial clustering
      ##
      ## Either user defined or obtained (the default)
      ## by kmeans on the variational parameters of the means of a fully parametrized PLN
      if (is.list(control$init_cl)) {
        stopifnot(length(control$init_cl) == length(nb_blocks),
                  all(sapply(control$init_cl, length) == private$p))
        blocks <- control$init_cl
      } else {
        blocks <- switch(control$init_cl,
          "kmeans" = {
            Means <- t(myPLN_init$var_par$M)
            blocks <- lapply(nb_blocks, function(k) {
              if (k == private$p) res <- 1:private$p
              else res <- kmeans(Means, centers = k, nstart = 30)$cl
              res
            })
          },
          "clustofvar" = hclustvar(myPLN_init$var_par$M) %>% cutree(nb_blocks) %>% as.data.frame() %>% as.list(),
          "hclust" = hclust(as.dist(1 - cov2cor(myPLN_init$model_par$Sigma)), method = "complete") %>% cutree(nb_blocks) %>% as.data.frame() %>% as.list()
        )
      }
      blocks <- blocks %>% map(as_indicator) %>% map(.check_boundaries)

      ## instantiate as many models as clusterings
      self$models <-
        map2(blocks, sparsity, function(blocks_, sparsity_) {
          if (sparsity_ > 0) {
            model <- PLNblockbisfit_sparse$new(blocks_, sparsity_, responses, covariates, offsets, weights, formula, control)
          } else {
            model <- PLNblockbisfit$new(blocks_, responses, covariates, offsets, weights, formula, control)
          }
        })
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Optimization ----------------------
    #' @description Call to the C++ optimizer on all models of the collection
    #' @param config a list for controlling the optimization.
    optimize_sequentially = function(config) {
      ## go along modes by decrezasing group sizes
      models_order <- order(self$nb_blocks, decreasing = TRUE)
      for (m in seq_along(models_order))  {
        if (config$trace == 1) {
          cat("\tNumber of blocks =", self$models[[models_order[m]]]$nb_block, "\r")
          flush.console()
        }
        if (config$trace > 1) {
          cat("\tNumber of blocks =", self$models[[models_order[m]]]$nb_block, "- iteration:")
        }
        self$models[[models_order[m]]]$optimize(self$responses, self$covariates, self$offsets, self$weights, config)
        ## Save time by starting the optimization of model m + 1  with optimal parameters of model m
        if (m < length(self$nb_blocks)) {
          blocks <-
            self$models[[models_order[m]]]$var_par$M  %>%
            kmeansvar(init = self$models[[models_order[m+1]]]$nb_block, nstart = 30) %>%
            pluck("cluster") %>% as_indicator() %>% .check_boundaries()
          self$models[[models_order[m+1]]]$update(
            B = self$models[[models_order[m]]]$model_par$B,
            M = self$models[[models_order[m]]]$var_par$M %*% blocks ,
            S = self$models[[models_order[m]]]$var_par$S %*% blocks,

            ###
            Mu = self$models[[models_order[m]]]$var_par$M,
            Delta = self$models[[models_order[m]]]$var_par$S,

          )
        }
        if (config$trace > 1) {
          cat("\r                                                                                    \r")
          flush.console()
        }
      }
    },
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Optimization ----------------------
    #' @description Call to the C++ optimizer on all models of the collection
    #' @param config a list for controlling the optimization.
    #' @importFrom purrr pluck
    #' @importFrom furrr future_map
    optimize = function(config) {
      self$models <- lapply(self$models, function(model) {
        if (config$trace >= 1) {
          cat("\tnumber of blocks =", model$nb_block, "\r")
          flush.console()
        }
        model$optimize(self$responses, self$covariates, self$offsets, self$weights, config)
        model
      })
    },

    #' @description Extract best model in the collection
    #' @param crit a character for the criterion used to performed the selection. Either
    #' "BIC", "ICL" or "loglik". Default is `ICL`
    #' @return a [`PLNblockbisfit`] object
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
    #' @description Display various outputs (goodness-of-fit criteria, robustness, diagnostic) associated with a collection of PLNblockbiss fits (a [`PLNblockbisfamily`])
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


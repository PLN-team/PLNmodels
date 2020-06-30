#' An R6 Class to represent a collection of PLNmixturefit
#'
#' @description The function [PLNmixture()] produces an instance of this class.
#'
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for [getBestModel()], [getModel()] and [`plot()`][plot.PLNPCAfamily()].
#'
#' @param responses the matrix of responses common to every models
#' @param covariates the matrix of covariates common to every models
#' @param offsets the matrix of offsets common to every models
#' @param control a list for controlling the optimization. See details.
#' @param clusters the dimensions of the successively fitted models
#' @param model model used for fitting, extracted from the formula in the upper-level call
#' @param control a list for controlling the optimization. See details.
#' @param xlevels named listed of factor levels included in the models, extracted from the formula in the upper-level call #'
#' @include PLNfamily-class.R
#' @importFrom R6 R6Class
#' @import ggplot2
#' @seealso The function \code{\link{PLNmixture}}, the class \code{\link[=PLNmixturefit]{PLNmixturefit}}
PLNmixturefamily <-
  R6Class(classname = "PLNmixturefamily",
    inherit = PLNfamily,
    active = list(
      clusters = function() private$params
    ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    public = list(
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Creation -----------------------
    #' @description Initialize all models in the collection.
      initialize = function(clusters, responses, covariates, offsets, model, xlevels, control) {

        ## initialize the required fields
        super$initialize(responses, covariates, offsets, rep(1, nrow(responses)), control)
        private$params <- clusters

        if (control$trace > 0) cat("\n Perform kmeans on the latent layer of the inceptive model...")
        myPLN <- PLNfit$new(responses, covariates, offsets, rep(1, nrow(responses)), model, xlevels, control)
        myPLN$optimize(responses, covariates, offsets, rep(1, nrow(responses)), control)

        if(control$covariance == 'spherical')
          Sbar <- c(myPLN$var_par$S2) * myPLN$p
        else
          Sbar <- rowSums(myPLN$var_par$S2)
        D <- sqrt(as.matrix(dist(myPLN$var_par$M)^2) + outer(Sbar,rep(1,myPLN$n)) + outer(rep(1, myPLN$n), Sbar))

        if (is.numeric(control$init_cl)) {
          clusterings <- control$init_cl
        } else if (is.character(control$init_cl)) {
          clusterings <-switch(control$init_cl,
            "kmeans"  = lapply(clusters, function(k) kmeans(D, centers = k, nstart = 30)$cl),
            "ward.D2" = D %>% as.dist() %>% hclust(method = "ward.D2") %>% cutree(clusters) %>% as.data.frame() %>% as.list()
          )
        }

        self$models <-
            clusterings %>%
            lapply(as_indicator) %>%
            lapply(.check_boundaries) %>%
            lapply(function(Z) PLNmixturefit$new(responses, covariates, offsets, Z, model, xlevels, control))
      },
      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Optimization ----------------------
      #' @description Call to the optimizer on all models of the collection
      optimize = function(control) {
        ## go along the number of clusters (i.e the models)
        for (m in seq_along(self$models))  {
          if (control$trace == 1) {
            cat("\tnumber of cluster =", self$models[[m]]$k, "\r")
            flush.console()
          }
          if (control$trace > 1) {
            cat("\tnumber of cluster =", self$models[[m]]$k, "- iteration:")
          }

          self$models[[m]]$optimize(self$responses, self$covariates, self$offsets, control)

          if (control$trace > 1) {
            cat("\r                                                                                    \r")
            flush.console()
          }

        }
      },
      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Graphical methods -------------
      #' @description
      #' Lineplot of selected criteria for all models in the collection
      #' @param criteria A valid model selection criteria for the collection of models. Any of "loglik", "BIC" or "ICL" (all).
      #' @param annotate Logical. Should R2 be added to the plot (defaults to `TRUE`)
      #' @return A [`ggplot2`] object
      plot = function(criteria = c("loglik", "BIC", "ICL"), annotate = FALSE) {
        # vlines <- sapply(intersect(criteria, c("BIC", "ICL")), function(crit) self$getBestModel(crit)$k)
        p <- super$plot(criteria, annotate) + xlab("# of clusters") # + geom_vline(xintercept = vlines, linetype = "dashed", alpha = 0.25)
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
          ggtitle("Objective along the alternate algorithm") + xlab("iteration (+ changes of model)") + theme_bw()
        p
      },

      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Extractors   -------------------
      #' @description Extract best model in the collection
      #' @param crit a character for the criterion used to performed the selection. Either
      #' "BIC", "ICL", or "R_squared". Default is `BIC`
      #' @return a [`PLNmixturefit`] object
      getBestModel = function(crit = c("BIC", "ICL", "R_squared")){
        crit <- match.arg(crit)
        stopifnot(!anyNA(self$criteria[[crit]]))
        id <- 1
        if (length(self$criteria[[crit]]) > 1) {
          id <- which.max(self$criteria[[crit]])
        }
        model <- self$models[[id]]$clone()
        model
      },
      show = function() {
        super$show()
        cat(" Task: Mixture Model \n")
        cat("========================================================\n")
        cat(" - Number of clusters considered: from", min(self$clusters), "to", max(self$clusters),"\n")
        cat(" - Best model (regarding BIC): cluster =", self$getBestModel("BIC")$k, "- R2 =", round(self$getBestModel("BIC")$R_squared, 2), "\n")
        cat(" - Best model (regarding ICL): cluster =", self$getBestModel("ICL")$k, "- R2 =", round(self$getBestModel("ICL")$R_squared, 2), "\n")
      },
      print = function() self$show()
    )
)


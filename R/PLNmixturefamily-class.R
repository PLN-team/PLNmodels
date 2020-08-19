#' An R6 Class to represent a collection of PLNmixturefit
#'
#' @description The function [PLNmixture()] produces an instance of this class.
#'
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for [getBestModel()], [getModel()] and [`plot()`][plot.PLNmixturefamily()].
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
#' @importFrom purrr map map_dbl map_int
#' @import ggplot2
#' @seealso The function \code{\link{PLNmixture}}, the class \code{\link[=PLNmixturefit]{PLNmixturefit}}
PLNmixturefamily <-
  R6Class(classname = "PLNmixturefamily",
    inherit = PLNfamily,
    active = list(
      clusters = function() private$params
    ),
    private = list(
      model = NULL,
      xlevels = NULL
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
        private$params  <- clusters
        private$model   <- model
        private$xlevels <- xlevels

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
            map(as_indicator) %>%
            map(.check_boundaries) %>%
            map(function(Z) PLNmixturefit$new(responses, covariates, offsets, Z, model, xlevels, control))
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
    smooth = function(control) {
        if (control$trace > 0) control$trace <- TRUE else control$trace <- FALSE
        for (i in 1:control$iterates) {
          if (control$smoothing %in% c('forward' , 'both')) self$smooth_forward(control)
          if (control$smoothing %in% c('backward', 'both')) self$smooth_backward(control)
        }
      },
    smooth_forward = function(control) {
      trace <- control$trace > 0; control$trace <- FALSE
      if (trace) cat("   Going forward ")
      for (i in self$clusters[-length(self$clusters)]) {
        if (trace) cat("+")
        cl0 <- self$models[[i]]$memberships
        if (length(unique(cl0)) == i) { # when would this not happens ?
          candidates <- mclapply(1:i, function(j) {
            cl <- cl0
            J  <- which(cl == j)
            if (length(J) > 1) {
              J1 <- base::sample(J, floor(length(J)/2))
              J2 <- setdiff(J, J1)
              cl[J1] <- j; cl[J2] <- i + 1
              # model <- self$models[[i + 1]]$clone()
              # model$posteriorProb <- as_indicator(cl)
              model <- PLNmixturefit$new(self$responses, self$covariates, self$offsets, as_indicator(cl), private$model, private$xlevels, control)
              model$optimize(self$responses, self$covariates, self$offsets, control)
            } else {
              model <- self$models[[i + 1]]$clone()
            }
            model
          }, mc.cores = control$cores)
          best_one <- candidates[[which.max(map_dbl(candidates, 'loglik'))]]
          if (best_one$loglik > self$models[[i + 1]]$loglik)
            self$models[[i + 1]] <- best_one
        }
      }
  if (trace) cat("\r                                                                                                    \r")
      },
      smooth_backward = function(control) {
        trace <- control$trace > 0; control$trace <- FALSE
        if (trace) cat("   Going backward ")
        for (i in rev(self$clusters[-1])) {
          if (trace) cat('+')
          cl0 <- factor(self$models[[i]]$memberships)
          if (nlevels(cl0) == i) {
            candidates <- mclapply(combn(i, 2, simplify = FALSE), function(couple) {
              cl_fusion <- cl0
              levels(cl_fusion)[which(levels(cl_fusion) == paste(couple[1]))] <- paste(couple[2])
              levels(cl_fusion) <- as.character(1:(i - 1))
              model <- PLNmixturefit$new(self$responses, self$covariates, self$offsets, as_indicator(cl_fusion), private$model, private$xlevels, control)
              # model <- self$models[[i - 1]]$clone()
              # model$posteriorProb <- as_indicator(cl_fusion)
              model$optimize(self$responses, self$covariates, self$offsets, control)
              model
            }, mc.cores = control$cores)
            best_one <- candidates[[which.max(map_dbl(candidates, 'loglik'))]]
            if (best_one$loglik > self$models[[i - 1]]$loglik)
              self$models[[i - 1]] <- best_one
          }
        }
        if (trace) cat("\r                                                                                                    \r")
      },
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Graphical methods -------------
      #' @description
      #' Lineplot of selected criteria for all models in the collection
      #' @param criteria A valid model selection criteria for the collection of models. Any of "loglik", "BIC" or "ICL" (all).
      #' @param annotate Logical. Should R2 be added to the plot (defaults to `FALSE`)
      #' @return A [`ggplot2`] object
      plot = function(criteria = c("loglik", "BIC", "ICL"), annotate = FALSE) {
        vlines <- map_int(intersect(criteria, c("BIC", "ICL")), function(crit) self$getBestModel(crit)$k)
        p <- super$plot(criteria, annotate) + xlab("# of clusters") + geom_vline(xintercept = vlines, linetype = "dashed", alpha = 0.25)
        p
       },
      #' @description Plot objective value of the optimization problem along the penalty path
      #' @return a [`ggplot`] graph
      plot_objective = function() {
        objective <- self$models %>% map('optim_par') %>% map('objective') %>% unlist
        changes   <- self$models %>% map('optim_par') %>% map('outer_iterations') %>% unlist %>% cumsum
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
      #' "BIC", "ICL", or "loglik". Default is `BIC`
      #' @return a [`PLNmixturefit`] object
      getBestModel = function(crit = c("BIC", "ICL")){
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
        cat(" - Best model (regarding BIC): cluster =", self$getBestModel("BIC")$k, "\n")
        cat(" - Best model (regarding ICL): cluster =", self$getBestModel("ICL")$k, "\n")
      },
      print = function() self$show()
    )
)


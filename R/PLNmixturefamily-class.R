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
#' @param formula model formula used for fitting, extracted from the formula in the upper-level call
#' @param control a list for controlling the optimization. See details.
#' @include PLNfamily-class.R
#' @importFrom R6 R6Class
#' @importFrom purrr map map_dbl map_int
#' @import ggplot2
#' @seealso The function \code{\link{PLNmixture}}, the class \code{\link[=PLNmixturefit]{PLNmixturefit}}
PLNmixturefamily <-
  R6Class(classname = "PLNmixturefamily",
    inherit = PLNfamily,
    active = list(
      #' @field clusters vector indicating the number of clusters considered is the successively fitted models
      clusters = function() private$params
    ),
    private = list(
      formula = NULL,
      #' @description helper function for forward smoothing: split a group
      add_one_cluster = function(model, k = NULL, control) {
        ## Control options
        control$trace <- FALSE
        config_fast <- control$config_optim
        config_fast$maxit_out <- 2

        ## Effective number of clusters (remove empty classes) and current clustering with clusters numbered in 1:k (with no gaps)
        cl  <- model$memberships
        k <- length(unique(cl))
        cl <- factor(cl) %>% as.integer() ## hacky way of getting memberships in 1:k
        ## all best split according to kmeans
        data_split <- model$latent_pos %>% as.data.frame() %>% split(cl)
        cl_splitable <- (1:k)[tabulate(cl) >= 3]
        cl_split <- vector("list", k)
        cl_split[cl_splitable] <- data_split[cl_splitable] %>% map(kmeans, 2, nstart = 10) %>% map("cluster")

        ## Reformating into indicator of clusters
        tau_candidates <- map(cl_splitable, function(k_)  {
          split <- cl_split[[k_]]
          split[cl_split[[k_]] == 1] <- k_
          split[cl_split[[k_]] == 2] <- k + 1
          candidate <- cl
          candidate[candidate == k_] <- split
          candidate
        }) %>% map(as_indicator)

        loglik_candidates <- future.apply::future_lapply(tau_candidates, function(tau_) {
          model <- PLNmixturefit$new(self$responses, self$covariates, self$offsets, tau_, private$formula, control)
          model$optimize(self$responses, self$covariates, self$offsets, config_fast)
          model$loglik
        }, future.seed = TRUE, future.scheduling = structure(TRUE, ordering = "random")) %>% unlist()

        best_one <- PLNmixturefit$new(self$responses, self$covariates, self$offsets, tau_candidates[[which.max(loglik_candidates)]], private$formula, control)
        best_one$optimize(self$responses, self$covariates, self$offsets, control$config_optim)
        best_one
      },
      smooth_forward  = function(control) {
        ## setup and verbosity levels
        trace <- control$trace > 0

        if (trace) cat("   Going forward ")
        for (model_index in head(seq_along(self$clusters), -1)) {
          if (trace) cat("+")
          ## effective and target number of clusters
          effective_k <- length(unique(self$models[[model_index]]$memberships)) ## instead of self$clusters[model_index] to protect against empty classes
          target_k <- self$clusters[model_index+1]
          candidate <- private$add_one_cluster(self$models[[model_index]],
                                               control = control)
          candidate_k <- length(unique(candidate$memberships))

          ## While the target number of clusters has not been reached and the number of clusters increases, keep going
          while (candidate_k > effective_k && candidate_k < target_k) {
            effective_k <- candidate_k
            candidate <- private$add_one_cluster(candidate, control = control)
            candidate_k <- length(unique(candidate$memberships))
          }

          ## Sanity check: candidate is viable only if it has the correct number of clusters
          if (candidate_k != target_k) next
          if (candidate$loglik > self$models[[model_index + 1]]$loglik) {
            self$models[[model_index + 1]] <- candidate
            # cat("found one")
          }

      }
      if (trace) cat("\r                                                                                                    \r")
      },
      remove_one_cluster = function(model, k = NULL, control) {
        ## Control options
        control$trace <- FALSE
        config_fast <- control$config_optim
        config_fast$maxit_out <- 2

        ## number of clusters
        if (is.null(k)) k <- length(model$components)

        tau <- model$posteriorProb
        tau_candidates <- lapply(combn(k, 2, simplify = FALSE), function(couple) {
          i <- min(couple); j <- max(couple)
          tau_merged <- tau[, -j, drop = FALSE]
          tau_merged[, i] <- rowSums(tau[, c(i,j)])
          tau_merged
        })

        loglik_candidates <- future.apply::future_lapply(tau_candidates, function(tau_) {
          model <- PLNmixturefit$new(self$responses, self$covariates, self$offsets, tau_, private$formula, control)
          model$optimize(self$responses, self$covariates, self$offsets, config_fast)
          model$loglik
        }, future.seed = TRUE, future.scheduling = structure(TRUE, ordering = "random")) %>% unlist()

        best_one <- PLNmixturefit$new(self$responses, self$covariates, self$offsets, tau_candidates[[which.max(loglik_candidates)]], private$formula, control)
        best_one$optimize(self$responses, self$covariates, self$offsets, control$config_optim)
        best_one
      },
      smooth_backward = function(control) {
        trace <- control$trace > 0
        if (trace) cat("   Going backward ")
        for (model_index in rev(seq_along(self$clusters)[-1])) {
          if (trace) cat('+')

          ## current and target number of clusters
          target_k <- self$clusters[model_index - 1]
          current_k <- self$clusters[model_index]
          candidate <- private$remove_one_cluster(self$models[[model_index]], k = current_k, control = control)
          current_k <- current_k - 1

          ## The number of clusters always decreases by one after merging (unlike). Keep going until the target number of clusters has been reached
          while (current_k > target_k) {
            candidate <- private$remove_one_cluster(candidate, k = current_k, control = control)
            current_k <- current_k - 1
          }

          ## Sanity check: candidate is viable only if it has the correct number of clusters
          if (current_k != target_k) next ## should never happen
          if (candidate$loglik > self$models[[model_index - 1]]$loglik) {
            self$models[[model_index - 1]] <- candidate
            # cat("found one")
          }
        }
        if (trace) cat("\r                                                                                                    \r")
      }

    ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    public = list(
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Creation -----------------------
    #' @description Initialize all models in the collection.
      initialize = function(clusters, responses, covariates, offsets, formula, control) {

        ## initialize the required fields
        super$initialize(responses, covariates, offsets, rep(1, nrow(responses)), control)
        private$params  <- clusters
        private$formula <- formula

        ## Default clustering is obtained by performing kmeans or CAH on the variational parameters of the means of a fully parametrized PLN
        if (is.numeric(control$config_optim$init_cl)) {
          clusterings <- control$init_cl
        } else if (is.character(control$init_cl)) {
          myPLN <- PLNfit$new(responses, covariates, offsets, rep(1, nrow(responses)), formula, control)
          myPLN$optimize(responses, covariates, offsets, rep(1, nrow(responses)), control$config_optim)
          Sbar <- rowSums(myPLN$var_par$S2)
          D <- sqrt(as.matrix(dist(myPLN$var_par$M)^2) + outer(Sbar,rep(1,myPLN$n)) + outer(rep(1, myPLN$n), Sbar))
          clusterings <-switch(control$init_cl,
            "kmeans"  = lapply(clusters, function(k) kmeans(D, centers = k, nstart = 30)$cl),
            "ward.D2" = D %>% as.dist() %>% hclust(method = "ward.D2") %>% cutree(clusters) %>% as.data.frame() %>% as.list()
          )
        }
        self$models <-
          clusterings %>%
            map(as_indicator) %>%
            map(.check_boundaries) %>%
            map(function(Z) {
              PLNmixturefit$new(responses, covariates, offsets, Z, formula, control)}
            )
      },
      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Optimization ----------------------
      #' @description Call to the optimizer on all models of the collection
      #' @param config a list for controlling the optimization
      optimize = function(config) {
        ## go along the number of clusters (i.e the models)
         self$models <- future.apply::future_lapply(self$models, function(model) {
          if (config$trace == 1) {
            cat("\tnumber of cluster =", model$k, "\r")
            flush.console()
          }
          model$optimize(self$responses, self$covariates, self$offsets, config)
          if (config$trace > 1) {
            cat("\r                                                                                    \r")
            flush.console()
          }
          model
        }, future.seed = TRUE)
      },
      #' @description
      #' function to restart clustering to avoid local minima by smoothing the loglikelihood values as a function of the number of clusters
      #' @param control a list to control the smoothing process
      smooth = function(control) {
        if (control$trace > 0) control$trace <- TRUE else control$trace <- FALSE
        for (i in seq_len(control$config_optim$it_smooth)) {
          if (control$smoothing %in% c('backward', 'both')) private$smooth_backward(control)
          if (control$smoothing %in% c('forward' , 'both')) private$smooth_forward(control)
        }
      },
      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Graphical methods -------------
      #' @description
      #' Lineplot of selected criteria for all models in the collection
      #' @param criteria A valid model selection criteria for the collection of models. Any of "loglik", "BIC" or "ICL" (all).
      #' @param reverse A logical indicating whether to plot the value of the criteria in the "natural" direction
      #' (loglik - 0.5 penalty) or in the "reverse" direction (-2 loglik + penalty). Default to FALSE, i.e use the
      #' natural direction, on the same scale as the log-likelihood..
      #' @return A [`ggplot2`] object
      plot = function(criteria = c("loglik", "BIC", "ICL"), reverse = FALSE) {
        vlines <- map_int(intersect(criteria, c("BIC", "ICL")), function(crit) self$getBestModel(crit)$k)
        p <- super$plot(criteria, reverse) + xlab("# of clusters") + geom_vline(xintercept = vlines, linetype = "dashed", alpha = 0.25)
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
      #' "BIC", "ICL" or "loglik". Default is `ICL`
      #' @return a [`PLNmixturefit`] object
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
      ## Print methods ---------------------
      #' @description User friendly print method
      show = function() {
        super$show()
        cat(" Task: Mixture Model \n")
        cat("========================================================\n")
        cat(" - Number of clusters considered: from", min(self$clusters), "to", max(self$clusters),"\n")
        cat(" - Best model (regarding BIC): cluster =", self$getBestModel("BIC")$k, "\n")
        cat(" - Best model (regarding ICL): cluster =", self$getBestModel("ICL")$k, "\n")
      },
      #' @description User friendly print method
      print = function() self$show()
    )
)


#' An R6 Class to represent a collection of PLNnetworkfit
#'
#' @description The function [PLNnetwork()] produces an instance of this class.
#'
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for [getBestModel()],
#' [getModel()] and [plot()][plot.PLNnetworkfamily()]
#'
## Parameters shared by many methods
#' @param penalties a vector of positive real number controlling the level of sparsity of the underlying network.
#' @param responses the matrix of responses common to every models
#' @param covariates the matrix of covariates common to every models
#' @param offsets the matrix of offsets common to every models
#' @param weights the vector of observation weights
#' @param formula model formula used for fitting, extracted from the formula in the upper-level call
#' @param control a list for controlling the optimization.
#' @param var value of the parameter (`rank` for PLNPCA, `sparsity` for PLNnetwork) that identifies the model to be extracted from the collection. If no exact match is found, the model with closest parameter value is returned with a warning.
#' @param index Integer index of the model to be returned. Only the first value is taken into account
#'
#' @include PLNfamily-class.R
#' @importFrom R6 R6Class
#' @importFrom glassoFast glassoFast
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' fits <- PLNnetwork(Abundance ~ 1, data = trichoptera)
#' class(fits)
#' @seealso The function [PLNnetwork()], the class [`PLNnetworkfit`]
PLNnetworkfamily <- R6Class(
  classname = "PLNnetworkfamily",
  inherit = PLNfamily,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ------
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Creation functions ----------------
    #' @description Initialize all models in the collection
    #' @return Update current [`PLNnetworkfit`] with smart starting values
    initialize = function(penalties, responses, covariates, offsets, weights, formula, control) {

      ## Initialize fields shared by the super class
      super$initialize(responses, covariates, offsets, weights, control)

      ## A basic model for inception, useless one is defined by the user
### TODO check if it is useful
      if (is.null(control$inception)) {
        myPLN <- PLNfit$new(responses, covariates, offsets, weights, formula, control)
        myPLN$optimize(responses, covariates, offsets, weights, control$config_optim)
        control$inception <- myPLN
      }

      if (is.null(control$penalty_weights))
        control$penalty_weights <- matrix(1, ncol(responses), ncol(responses))
      ## Get the number of penalty
      if (is.null(penalties)) {
        if (is.list(control$penalty_weights))
          control$n_penalties <- length(control$penalty_weights)
      } else {
        control$n_penalties <- length(penalties)
      }
      ## Define a matrix of weights for each penalty
      if (!is.list(control$penalty_weights))
        list_penalty_weights <- rep(list(control$penalty_weights), control$n_penalties)
      else
        list_penalty_weights <- control$penalty_weights

      ## Get an appropriate grid of penalties
      if (is.null(penalties)) {
        if (control$trace > 1) cat("\n Recovering an appropriate grid of penalties.")
        max_pen <- list_penalty_weights %>%
          map(~ myPLN$model_par$Sigma / .x) %>%
          map_dbl(~ max(abs(.x[upper.tri(.x, diag = control$penalize_diagonal)]))) %>%
          max()
        penalties <- 10^seq(log10(max_pen), log10(max_pen*control$min_ratio), len = control$n_penalties)
      } else {
        if (control$trace > 1) cat("\nPenalties already set by the user")
        stopifnot(all(penalties > 0))
      }
      ## Sort eh penalty in decreasing order
      o <- order(penalties, decreasing = TRUE)
      private$params <- penalties[o]
      list_penalty_weights <- list_penalty_weights[o]

      ## instantiate as many models as penalties
      control$trace <- 0
      self$models <- map2(private$params, list_penalty_weights, function(penalty, penalty_weights) {
        PLNnetworkfit$new(penalty, penalty_weights, responses, covariates, offsets, weights, formula, control)
      })

    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Optimization ----------------------
    #' @description Call to the C++ optimizer on all models of the collection
    #' @param config a list for controlling the optimization.
    optimize = function(config) {
      ## Go along the penalty grid (i.e the models)
      for (m in seq_along(self$models))  {

        if (config$trace == 1) {
          cat("\tsparsifying penalty =", self$models[[m]]$penalty, "\r")
          flush.console()
        }
        if (config$trace > 1) {
          cat("\tsparsifying penalty =", self$models[[m]]$penalty, "- iteration:")
        }
        self$models[[m]]$optimize(self$responses, self$covariates, self$offsets, self$weights, config)
        ## Save time by starting the optimization of model m + 1  with optimal parameters of model m
        if (m < length(self$penalties))
          self$models[[m + 1]]$update(
            B = self$models[[m]]$model_par$B,
            M = self$models[[m]]$var_par$M,
            S = self$models[[m]]$var_par$S
          )

        if (config$trace > 1) {
          cat("\r                                                                                    \r")
          flush.console()
        }

      }

    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Stability -------------------------
    #' @description Compute the stability path by stability selection
    #' @param subsamples a list of vectors describing the subsamples. The number of vectors (or list length) determines the number of subsamples used in the stability selection. Automatically set to 20 subsamples with size \code{10*sqrt(n)} if \code{n >= 144} and \code{0.8*n} otherwise following Liu et al. (2010) recommendations.
    #' @param control a list controlling the main optimization process in each call to PLNnetwork. See [PLNnetwork()] for details.
    stability_selection = function(subsamples = NULL, control = PLNnetwork_param()) {

      ## select default subsamples according
      if (is.null(subsamples)) {
        subsample.size <- round(ifelse(private$n >= 144, 10*sqrt(private$n), 0.8*private$n))
        subsamples <- replicate(20, sample.int(private$n, subsample.size), simplify = FALSE)
      }

      ## got for stability selection
      cat("\nStability Selection for PLNnetwork: ")
      cat("\nsubsampling: ")

      stabs_out <- future.apply::future_lapply(subsamples, function(subsample) {
        cat("+")
        inception_ <- self$getModel(self$penalties[1])
        inception_$update(
          M  = inception_$var_par$M[subsample, ],
          S  = inception_$var_par$S[subsample, ]
        )

        ## force some control parameters
        control$inception = inception_
        control$penalty_weights = map(self$models, "penalty_weights")
        control$penalize_diagonal = (sum(diag(inception_$penalty_weights)) != 0)
        control$trace <- 0
        control$config_optim$trace <- 0

        myPLN <- PLNnetworkfamily$new(penalties  = self$penalties,
                                      responses  = self$responses [subsample, , drop = FALSE],
                                      covariates = self$covariates[subsample, , drop = FALSE],
                                      offsets    = self$offsets   [subsample, , drop = FALSE],
                                      formula    = private$formula,
                                      weights    = self$weights   [subsample], control = control)

        myPLN$optimize(control$config_optim)
        nets <- do.call(cbind, lapply(myPLN$models, function(model) {
          as.matrix(model$latent_network("support"))[upper.tri(diag(private$p))]
        }))
        nets
      }, future.seed = TRUE, future.scheduling = structure(TRUE, ordering = "random"))

      prob <- Reduce("+", stabs_out, accumulate = FALSE) / length(subsamples)
      ## formatting/tyding
      node_set <- colnames(self$getModel(index = 1)$model_par$B)
      colnames(prob) <- self$penalties
      private$stab_path <- prob %>%
        as.data.frame() %>%
        mutate(Edge = 1:n()) %>%
        gather(key = "Penalty", value = "Prob", -Edge) %>%
        mutate(Penalty = as.numeric(Penalty),
               Node1   = node_set[edge_to_node(Edge)$node1],
               Node2   = node_set[edge_to_node(Edge)$node2],
               Edge    = paste0(Node1, "|", Node2))

      invisible(subsamples)
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Extractors ------------------------
    #' @description Extract the regularization path of a [`PLNnetworkfamily`]
    #' @param precision Logical. Should the regularization path be extracted from the precision matrix Omega (`TRUE`, default) or from the variance matrix Sigma (`FALSE`)
    #' @param corr Logical. Should the matrix be transformed to (partial) correlation matrix before extraction? Defaults to `TRUE`
    coefficient_path = function(precision = TRUE, corr = TRUE) {
      lapply(self$penalties, function(x) {
        if (precision) {
          G <- self$getModel(x)$model_par$Omega
        } else {
          G <- self$getModel(x)$model_par$Sigma
          dimnames(G) <- dimnames(self$getModel(x)$model_par$Omega)
        }
        if (corr) {
          G <- ifelse(precision, -1, 1) * G / tcrossprod(sqrt(diag(G)))
        }
        setNames(
          cbind(
            expand.grid(colnames(G), rownames(G)),
            as.vector(G)), c("Node1", "Node2", "Coeff")
        ) %>%
          mutate(Penalty = x,
                 Node1   = as.character(Node1),
                 Node2   = as.character(Node2),
                 Edge    = paste0(Node1, "|", Node2)) %>%
          filter(Node1 < Node2)
      }) %>% bind_rows()
    },

    #' @description Extract the best network in the family according to some criteria
    #' @param crit character. Criterion used to perform the selection. Is "StARS" is chosen but `$stability` field is empty, will compute stability path.
    #' @param stability Only used for "StARS" criterion. A scalar indicating the target stability (= 1 - 2 beta) at which the network is selected. Default is `0.9`.
    getBestModel = function(crit = c("BIC", "EBIC", "StARS"), stability = 0.9){
      crit <- match.arg(crit)
      if (crit == "StARS") {
        if (is.null(private$stab_path)) self$stability_selection()
        id_stars <- self$criteria %>%
          select(param, stability) %>% rename(Stability = stability) %>%
          filter(Stability > stability) %>%
          pull(param) %>% min() %>% match(self$penalties)
        model <- self$models[[id_stars]]$clone()
      } else {
        stopifnot(!anyNA(self$criteria[[crit]]))
        id <- 1
        if (length(self$criteria[[crit]]) > 1) {
          id <- which.max(self$criteria[[crit]])
        }
        model <- self$models[[id]]$clone()
      }
      model
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Graphical methods -----------------
    #' @description Display various outputs (goodness-of-fit criteria, robustness, diagnostic) associated with a collection of PLNnetwork fits (a [`PLNnetworkfamily`])
    #' @param criteria vector of characters. The criteria to plot in `c("loglik", "pen_loglik", "BIC", "EBIC")`. Defaults to all of them.
    #' @param reverse A logical indicating whether to plot the value of the criteria in the "natural" direction
    #' (loglik - 0.5 penalty) or in the "reverse" direction (-2 loglik + penalty). Default to FALSE, i.e use the
    #' natural direction, on the same scale as the log-likelihood..
    #' @param log.x logical: should the x-axis be represented in log-scale? Default is `TRUE`.
    #' @return a [`ggplot`] graph
    plot = function(criteria = c("loglik", "pen_loglik", "BIC", "EBIC"), reverse = FALSE, log.x = TRUE) {
      vlines <- sapply(intersect(criteria, c("BIC", "EBIC")) , function(crit) self$getBestModel(crit)$penalty)
      p <- super$plot(criteria, reverse) + xlab("penalty") + geom_vline(xintercept = vlines, linetype = "dashed", alpha = 0.25)
      if (log.x) p <- p + ggplot2::coord_trans(x = "log10")
      p
    },

  #' @description Plot stability path
  #' @param stability scalar: the targeted level of stability in stability plot. Default is `0.9`.
  #' @param log.x logical: should the x-axis be represented in log-scale? Default is `TRUE`.
  #' @return a [`ggplot`] graph
  plot_stars = function(stability = 0.9, log.x = TRUE) {
    if (anyNA(self$stability)) stop("stability selection has not yet been performed! Use stability_selection()")
    dplot <- self$criteria %>% select(param, density, stability) %>%
      rename(Penalty = param) %>%
      gather(key = "Metric", value = "Value", stability:density)
    penalty_stars <- dplot %>% filter(Metric == "stability" & Value >= stability) %>%
      pull(Penalty) %>% min()

    p <- ggplot(dplot, aes(x = Penalty, y = Value, group = Metric, color = Metric)) +
      geom_point() +  geom_line() + theme_bw() +
      ## Add information correspinding to best lambda
      geom_vline(xintercept = penalty_stars, linetype = 2) +
      geom_hline(yintercept = stability, linetype = 2) +
      annotate(x = penalty_stars, y = 0,
               label = paste("lambda == ", round(penalty_stars, 5)),
               parse = TRUE, hjust = -0.05, vjust = 0, geom = "text") +
      annotate(x = penalty_stars, y = stability,
               label = paste("stability == ", stability),
               parse = TRUE, hjust = -0.05, vjust = 1.5, geom = "text")
    if (log.x) p <- p + ggplot2::scale_x_log10() + annotation_logticks(sides = "b")
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
               label = paste("penalty=",format(self$criteria$param, digits = 1)), hjust = -.1, size = 3, alpha = 0.7) + theme_bw()
    p
  },


  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## Print methods ---------------------
  #' @description User friendly print method
  show = function() {
    super$show()
    cat(" Task: Network Inference \n")
    cat("========================================================\n")
    cat(" -", length(self$penalties) , "penalties considered: from", min(self$penalties), "to", max(self$penalties), "\n")
    cat(" - Best model (greater BIC): lambda =", format(self$getBestModel("BIC")$penalty, digits = 3), "\n")
    cat(" - Best model (greater EBIC): lambda =", format(self$getBestModel("EBIC")$penalty, digits = 3), "\n")
    if (!anyNA(self$criteria$stability))
    cat(" - Best model (regarding StARS): lambda =", format(self$getBestModel("StARS")$penalty, digits = 3), "\n")
  }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ------
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  private = list(
    stab_path = NULL # a field to store the stability path,
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## ACTIVE BINDINGS ------
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field penalties the sparsity level of the network in the successively fitted models
    penalties = function() private$params,
    #' @field stability_path the stability path of each edge as returned by the stars procedure
    stability_path = function() private$stab_path,
    #' @field stability mean edge stability along the penalty path
    stability = function() {
      if (!is.null(private$stab_path)) {
        stability <- self$stability_path %>%
          dplyr::select(Penalty, Prob) %>%
          group_by(Penalty) %>%
          summarize(Stability = 1 - mean(4 * Prob * (1 - Prob))) %>%
          arrange(desc(Penalty)) %>%
          pull(Stability)
      } else {
        stability <- rep(NA, length(self$penalties))
      }
      stability
    },
    #' @field criteria a data frame with the values of some criteria (approximated log-likelihood, (E)BIC, ICL and R2, stability) for the collection of models / fits
    #' BIC, ICL and EBIC are defined so that they are on the same scale as the model log-likelihood, i.e. with the form, loglik - 0.5 penalty
    criteria = function() {mutate(super$criteria, stability = self$stability)}
  )

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## END OF CLASS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
)


#' An R6 Class to represent a PLNfit in a sparse inverse covariance framework
#'
#' @description The function [PLNnetwork()] produces a collection of models which are instances of object with class [`PLNnetworkfit`].
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for [`plot()`][plot.PLNnetworkfit()] and methods inherited from [`PLNfit`].
#'
## Parameters common to all PLN-xx-fit methods (shared with PLNfit but inheritance does not work)
#' @param data a named list used internally to carry the data matrices
#' @param control a list for controlling the optimization.
#' @param nullModel null model used for approximate R2 computations. Defaults to a GLM model with same design matrix but no latent variable.
#' @param B matrix of regression coefficients
#' @param Sigma variance-covariance matrix of the latent variables
#' @param Omega precision matrix of the latent variables. Inverse of Sigma.
#'
## Parameters specific to PLNnetwork-fit methods
#' @param penalty a positive real number controlling the level of sparsity of the underlying network.
#' @param penalty_weights either a single or a list of p x p matrix of weights (default: all weights equal to 1) to adapt the amount of shrinkage to each pair of node. Must be symmetric with positive values.
#'
#' @include PLNnetworkfit-class.R
#' @examples
#' \dontrun{
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' nets <- PLNnetwork(Abundance ~ 1, data = trichoptera)
#' myPLNnet <- getBestModel(nets)
#' class(myPLNnet)
#' print(myPLNnet)
#' }
#' @seealso The function [PLNnetwork()], the class [`PLNnetworkfamily`]
PLNnetworkfit <- R6Class(
  classname = "PLNnetworkfit",
  inherit = PLNfit_fixedcov,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public  = list(
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Creation functions ----------------
    #' @description Initialize a [`PLNnetworkfit`] object
    initialize = function(data, control) {
      super$initialize(data$Y, data$X, data$O, data$w, data$formula, control)
      ## Default for penalty weights (if not already set)
      if (is.null(control$penalty_weights)) control$penalty_weights <- matrix(1, self$p, self$p)
      stopifnot(isSymmetric(control$penalty_weights), all(control$penalty_weights >= 0))
      if (!control$penalize_diagonal) diag(control$penalty_weights) <- 0
      private$lambda <- control$penalty
      private$rho    <- control$penalty_weights
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Optimization ----------------------
    #' @description Call to the C++ optimizer and update of the relevant fields
    #' @param config a list for controlling the optimization
    optimize = function(data, config) {
      cond <- FALSE; iter <- 0
      objective   <- numeric(config$maxit_out)
      convergence <- numeric(config$maxit_out)
      ## start from the standard PLN at initialization
      objective.old <- -self$loglik
      args <- list(data   = data,
                   params = list(B = private$B, M = private$M, S = private$S),
                   config = config)
      while (!cond) {
        iter <- iter + 1
        if (config$trace > 1) cat("", iter)
        ## CALL TO GLASSO TO UPDATE Omega
        glasso_out <- glassoFast::glassoFast(private$Sigma, rho = self$penalty * self$penalty_weights)
        if (anyNA(glasso_out$wi)) break
        private$Omega <- args$params$Omega <- Matrix::symmpart(glasso_out$wi)

        ## CALL TO NLOPT OPTIMIZATION TO UPDATE OTHER PARMAETERS
        optim_out <- do.call(private$optimizer$main, args)
        do.call(self$update, optim_out)

        ## Check convergence
        objective[iter]   <- -self$loglik # + self$penalty * sum(abs(private$Omega))
        convergence[iter] <- abs(objective[iter] - objective.old)/abs(objective[iter])
        if ((convergence[iter] < config$ftol_out) | (iter >= config$maxit_out)) cond <- TRUE

        ## Prepare next iterate
        args$params <- list(B = private$B, M = private$M, S = private$S)
        objective.old <- objective[iter]
      }

      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## OUTPUT
      private$Sigma <- Matrix::symmpart(glasso_out$w)
      private$monitoring$objective   <- objective[1:iter]
      private$monitoring$convergence <- convergence[1:iter]
      private$monitoring$iterations  <- iter
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Extractors ------------------------
    #' @description Extract interaction network in the latent space
    #' @param type edge value in the network. Can be "support" (binary edges), "precision" (coefficient of the precision matrix) or "partial_cor" (partial correlation between species)
    #' @importFrom Matrix Matrix
    #' @return a square matrix of size `PLNnetworkfit$n`
    latent_network = function(type = c("partial_cor", "support", "precision")) {
      net <- switch(
        match.arg(type),
        "support"     = 1 * (private$Omega != 0 & !diag(TRUE, ncol(private$Omega))),
        "precision"   = private$Omega,
        "partial_cor" = {
          tmp <- -private$Omega / tcrossprod(sqrt(diag(private$Omega))); diag(tmp) <- 1
          tmp
        }
      )
      ## Enforce sparse Matrix encoding to avoid downstream problems with igraph::graph_from_adjacency_matrix
      ## as it fails when given dsyMatrix objects
      Matrix(net, sparse = TRUE)
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Graphical methods------------------
    #' @description plot the latent network.
    #' @param type edge value in the network. Either "precision" (coefficient of the precision matrix) or "partial_cor" (partial correlation between species).
    #' @param output Output type. Either `igraph` (for the network) or `corrplot` (for the adjacency matrix)
    #' @param edge.color Length 2 color vector. Color for positive/negative edges. Default is `c("#F8766D", "#00BFC4")`. Only relevant for igraph output.
    #' @param node.labels vector of character. The labels of the nodes. The default will use the column names ot the response matrix.
    #' @param remove.isolated if `TRUE`, isolated node are remove before plotting. Only relevant for igraph output.
    #' @param layout an optional igraph layout. Only relevant for igraph output.
    #' @param plot logical. Should the final network be displayed or only sent back to the user. Default is `TRUE`.
    plot_network = function(type            = c("partial_cor", "support"),
                            output          = c("igraph", "corrplot"),
                            edge.color      = c("#F8766D", "#00BFC4"),
                            remove.isolated = FALSE,
                            node.labels     = NULL,
                            layout          = layout_in_circle,
                            plot = TRUE) {
      .plot_network(self$latent_network(match.arg(type)),
                   type            = match.arg(type),
                   output          = match.arg(output),
                   edge.color      = edge.color,
                   remove.isolated = remove.isolated,
                   node.labels     = node.labels,
                   layout          = layout,
                   plot            = plot)
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Print methods ---------------------
    #' @description User friendly print method
    show = function() {
      super$show(paste0("Poisson Lognormal with sparse inverse covariance (penalty = ", format(self$penalty,digits = 3),")\n"))
      cat("* Additional fields for sparse network\n")
      cat("    $EBIC, $density, $penalty \n")
      cat("* Additional S3 methods for network\n")
      cat("    plot.PLNnetworkfit() \n")
    }

  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    lambda = NA, # the sparsity tuning parameter
    rho    = NA  # the p x p penalty weight
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field vcov_model character: the model used for the residual covariance
    vcov_model = function() {"sparse"},
    #' @field penalty the global level of sparsity in the current model
    penalty         = function() {private$lambda},
    #' @field penalty_weights a matrix of weights controlling the amount of penalty element-wise.
    penalty_weights = function() {private$rho},
    #' @field n_edges number of edges if the network (non null coefficient of the sparse precision matrix)
    n_edges         = function() {sum(private$Omega[upper.tri(private$Omega, diag = FALSE)] != 0)},
    #' @field nb_param number of parameters in the current PLN model
    nb_param        = function() {self$p * self$d + self$p + self$n_edges},
    #' @field pen_loglik variational lower bound of the l1-penalized loglikelihood
    pen_loglik      = function() {self$loglik - private$lambda * sum(abs(private$Omega))},
    #' @field EBIC variational lower bound of the EBIC
    EBIC      = function() {self$BIC - .5 * ifelse(self$n_edges > 0, self$n_edges * log(.5 * self$p*(self$p - 1)/self$n_edges), 0)},
    #' @field density proportion of non-null edges in the network
    density   = function() {mean(self$latent_network("support"))},
    #' @field criteria a vector with loglik, penalized loglik, BIC, EBIC, ICL, R_squared, number of parameters, number of edges and graph density
    criteria  = function() {data.frame(super$criteria, n_edges = self$n_edges, EBIC = self$EBIC, pen_loglik = self$pen_loglik, density = self$density)}
  )
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## END OF CLASS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
)

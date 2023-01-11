#' An R6 Class to represent a PLNfit in a sparse inverse covariance framework
#'
#' @description The function [PLNnetwork()] produces a collection of models which are instances of object with class [`PLNnetworkfit`].
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for [`plot()`][plot.PLNnetworkfit()] and methods inherited from [`PLNfit`].
#'
## Parameters common to all PLN-xx-fit methods (shared with PLNfit but inheritance does not work)
#' @param responses the matrix of responses (called Y in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param covariates design matrix (called X in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param offsets offset matrix (called O in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param weights an optional vector of observation weights to be used in the fitting process.
#' @param formula model formula used for fitting, extracted from the formula in the upper-level call
#' @param control a list for controlling the optimization.
#' @param nullModel null model used for approximate R2 computations. Defaults to a GLM model with same design matrix but not latent variable.
#' @param B matrix of regression matrix
#' @param Sigma variance-covariance matrix of the latent variables
#' @param Omega precision matrix of the latent variables. Inverse of Sigma.
#'
## Parameters specific to PLNnetwork-fit methods
#' @param penalty a positive real number controlling the level of sparsity of the underlying network.
#' @param penalty_weights either a single or a list of p x p matrix of weights (default filled with 1) to adapt the amount of shrinkage to each pairs of node. Must be symmetric with positive values.
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
    initialize = function(penalty, penalty_weights, responses, covariates, offsets, weights, formula, control) {
      stopifnot(isSymmetric(penalty_weights), all(penalty_weights >= 0))
      super$initialize(responses, covariates, offsets, weights, formula, control)
      private$lambda <- penalty
      private$rho    <- penalty_weights
      if (!control$penalize_diagonal) diag(private$rho) <- 0
    },
    #' @description Update fields of a [`PLNnetworkfit`] object
    #' @param B matrix of regression matrix
    #' @param Sigma variance-covariance matrix of the latent variables
    #' @param Omega precision matrix of the latent variables. Inverse of Sigma.
    #' @param M     matrix of mean vectors for the variational approximation
    #' @param S     matrix of variance vectors for the variational approximation
    #' @param Z     matrix of latent vectors (includes covariates and offset effects)
    #' @param A     matrix of fitted values
    #' @param Ji    vector of variational lower bounds of the log-likelihoods (one value per sample)
    #' @param R2    approximate R^2 goodness-of-fit criterion
    #' @param monitoring a list with optimization monitoring quantities
    update = function(penalty=NA, B=NA, Sigma=NA, Omega=NA, M=NA, S=NA, Z=NA, A=NA, Ji=NA, R2=NA, monitoring=NA) {
      super$update(B = B, Sigma = Sigma, Omega = Omega, M, S = S, Z = Z, A = A, Ji = Ji, R2 = R2, monitoring = monitoring)
      if (!anyNA(penalty)) private$lambda <- penalty
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Optimization ----------------------
    #' @description Call to the C++ optimizer and update of the relevant fields
    #' @param config a list for controlling the optimization
    optimize = function(responses, covariates, offsets, weights, config) {
      cond <- FALSE; iter <- 0
      objective   <- numeric(config$maxit_out)
      convergence <- numeric(config$maxit_out)
      ## start from the standard PLN at initialization
      objective.old <- -self$loglik
      args <- list(data   = list(Y = responses, X = covariates, O = offsets, w = weights),
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
        objective[iter]   <- -self$loglik + self$penalty * sum(abs(private$Omega))
        convergence[iter] <- abs(objective[iter] - objective.old)/abs(objective[iter])
        if ((convergence[iter] < config$ftol_out) | (iter >= config$maxit_out)) cond <- TRUE

        ## Prepare next iterate
        args$params <- list(B = private$B, M = private$M, S = private$S)
        objective.old <- objective[iter]
      }

      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## OUTPUT
      private$Sigma <- Matrix::symmpart(glasso_out$w)
      private$monitoring$objective        <- objective[1:iter]
      private$monitoring$convergence      <- convergence[1:iter]
      private$monitoring$outer_iterations <- iter
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

      type <- match.arg(type)
      output <- match.arg(output)

      net <- self$latent_network(type)

      if (output == "igraph") {

        G <-  graph_from_adjacency_matrix(net, mode = "undirected", weighted = TRUE, diag = FALSE)

        if (!is.null(node.labels)) {
          igraph::V(G)$label <- node.labels
        } else {
          igraph::V(G)$label <- colnames(net)
        }
        ## Nice nodes
        V.deg <- degree(G)/sum(degree(G))
        igraph::V(G)$label.cex <- V.deg / max(V.deg) + .5
        igraph::V(G)$size <- V.deg * 100
        igraph::V(G)$label.color <- rgb(0, 0, .2, .8)
        igraph::V(G)$frame.color <- NA
        ## Nice edges
        igraph::E(G)$color <- ifelse(igraph::E(G)$weight > 0, edge.color[1], edge.color[2])
        if (type == "support")
          igraph::E(G)$width <- abs(igraph::E(G)$weight)
        else
          igraph::E(G)$width <- 15*abs(igraph::E(G)$weight)

        if (remove.isolated) {
          G <- delete.vertices(G, which(degree(G) == 0))
        }
        if (plot) plot(G, layout = layout)
      }
      if (output == "corrplot") {
        if (plot) {
          if (ncol(net) > 100)
            colnames(net) <- rownames(net) <- rep(" ", ncol(net))
          G <- net
          diag(net) <- 0
          corrplot(as.matrix(net), method = "color", is.corr = FALSE, tl.pos = "td", cl.pos = "n", tl.cex = 0.5, type = "upper")
        } else  {
          G <- net
        }
      }
      invisible(G)
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

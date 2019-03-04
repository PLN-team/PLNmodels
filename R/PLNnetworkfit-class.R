#' An R6 Class to represent a PLNfit in a sparse inverse covariance framework
#'
#' @description The function \code{\link{PLNnetwork}} produces a collection of models which are instances of object with class \code{PLNnetworkfit}.
#'
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for \code{\link[=plot.PLNnetworkfit]{plot.PLNnetworkfit}} + methods inherited from PLNfit.
#'
#' @field penalty the level of sparsity in the current model
#' @field model_par a list with the matrices associated with the estimated parameters of the pPCA model: Theta (covariates), Sigma (latent covariance) and Theta (latent precision matrix). Note Omega and Sigma are inverse of each other.
#' @field var_par a list with two matrices, M and S, which are the estimated parameters in the variational approximation
#' @field latent a matrix: values of the latent vector (Z in the model)
#' @field optim_par a list with parameters useful for monitoring the optimization
#' @field loglik variational lower bound of the loglikelihood
#' @field pen_loglik variational lower bound of the l1-penalized loglikelihood
#' @field BIC variational lower bound of the BIC
#' @field EBIC variational lower bound of the EBIC
#' @field ICL variational lower bound of the ICL
#' @field R_squared approximated goodness-of-fit criterion
#' @field nb_param number of parameters in the current PLN model
#' @field density proportion of non-null edges in the network
#' @field criteria a vector with loglik, BIC, ICL, R_squared and number of parameters
#' @include PLNnetworkfit-class.R
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' nets <- PLNnetwork(Abundance ~ 1, data = trichoptera)
#' myPLNnet <- getBestModel(nets)
#' class(myPLNnet)
#' print(myPLNnet)
#' @seealso The function \code{\link{PLNnetwork}}, the class \code{\link[=PLNnetworkfamily]{PLNnetworkfamily}}
PLNnetworkfit <-
  R6Class(classname = "PLNnetworkfit",
    inherit = PLNfit,
    public  = list(
      initialize = function(penalty, responses, covariates, offsets, weights, model, control) {
        super$initialize(responses, covariates, offsets, weights, model, control)
        private$lambda <- penalty
      },
      update = function(penalty=NA, Theta=NA, Sigma=NA, Omega=NA, M=NA, S=NA, Z=NA, A=NA, Ji=NA, R2=NA, monitoring=NA) {
        super$update(Theta = Theta, Sigma = Sigma, M, S = S, Z = Z, A = A, Ji = Ji, R2 = R2, monitoring = monitoring)
        if (!anyNA(penalty)) private$lambda <- penalty
        if (!anyNA(Omega))   private$Omega  <- Omega
      }
    ),
    private = list(
      Omega  = NA, # the p x p precision matrix
      lambda = NA  # the sparsity tuning parameter
    ),
    active = list(
      penalty         = function() {private$lambda},
      n_edges    = function() {sum(private$Omega[upper.tri(private$Omega, diag = FALSE)] != 0)},
      nb_param   = function() {self$p * self$d + self$n_edges},
      pen_loglik = function() {self$loglik - private$lambda * sum(abs(private$Omega))},
      model_par  = function() {
        par <- super$model_par
        par$Omega <- private$Omega
        par
      },
      EBIC      = function() {
        self$BIC - .5 * ifelse(self$n_edges > 0, self$n_edges * log(.5 * self$p*(self$p - 1)/self$n_edges), 0)
      },
      density   = function() {mean(self$latent_network("support"))},
      criteria  = function() {data.frame(super$criteria, n_edges = self$n_edges, EBIC = self$EBIC, pen_loglik = self$pen_loglik, density = self$density)}
    )
)


PLNnetworkfit$set("public", "optimize",
function(responses, covariates, offsets, weights, control) {

  ## shall we penalize the diagonal? in glassoFast
  rho <- matrix(self$penalty, self$p, self$p)
  if (!control$penalize_diagonal) diag(rho) <- 0

  cond <- FALSE; iter <- 0
  objective   <- numeric(control$maxit_out)
  convergence <- numeric(control$maxit_out)
  ## start from the standard PLN at initialization
  par0  <- c(private$Theta, private$M, private$S)
  Sigma <- private$Sigma
  objective.old <- -self$loglik
  while (!cond) {
    iter <- iter + 1
    if (control$trace > 1) cat("", iter)

    ## CALL TO GLASSO TO UPDATE Omega/Sigma
    glasso_out <- glassoFast::glassoFast(Sigma, rho = rho)
    if (anyNA(glasso_out$wi)) break
    Omega  <- glasso_out$wi ; if (!isSymmetric(Omega)) Omega <- Matrix::symmpart(Omega)

    ## CALL TO NLOPT OPTIMIZATION WITH BOX CONSTRAINT
    control$Omega <- Omega
    optim.out <- optim_sparse(par0, responses, covariates, offsets, weights, control)

    ## Check convergence
    objective[iter]   <- -sum(weights * optim.out$loglik) + self$penalty * sum(abs(Omega))
    convergence[iter] <- abs(objective[iter] - objective.old)/abs(objective[iter])
    if ((convergence[iter] < control$ftol_out) | (iter >= control$maxit_out)) cond <- TRUE

    ## Prepare next iterate
    Sigma <- optim.out$Sigma
    par0  <- c(optim.out$Theta, optim.out$M, optim.out$S)
    objective.old <- objective[iter]
  }

  ## ===========================================
  ## OUTPUT
  self$update(
    Theta = optim.out$Theta,
    Omega = Omega,
    Sigma = optim.out$Sigma,
    M = optim.out$M,
    S = optim.out$S,
    Z = optim.out$Z,
    A = optim.out$A,
    Ji = optim.out$loglik,
    monitoring = list(objective        = objective[1:iter],
                      convergence      = convergence[1:iter],
                      outer_iterations = iter,
                      inner_iterations = optim.out$iterations,
                      inner_status     = optim.out$status,
                      inner_message    = statusToMessage(optim.out$status)))

})

PLNnetworkfit$set("public", "postTreatment",
function(responses, covariates, offsets, weights, nullModel) {
  super$postTreatment(responses, covariates, offsets, weights, nullModel = nullModel)
  dimnames(private$Omega) <- dimnames(private$Sigma)
  colnames(private$S) <- 1:self$p
})

#' @importFrom Matrix Matrix
PLNnetworkfit$set("public", "latent_network",
  function(type = c("partial_cor", "support", "precision")) {
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
  }
)

# Plot the network (support of the inverse covariance) for a \code{PLNnetworkfit} object
PLNnetworkfit$set("public", "plot_network",
  function(type            = c("partial_cor", "support"),
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
})

PLNnetworkfit$set("public", "show",
function() {
  super$show(paste0("Poisson Lognormal with sparse inverse covariance (penalty = ", format(self$penalty,digits = 3),")\n"))
  cat("* Additional fields for sparse network\n")
  cat("    $EBIC, $density, $penalty \n")
  cat("* Additional S3 methods for network\n")
  cat("    plot.PLNnetworkfit() \n")
})

#' An R6 Class to represent a collection of PLNMMfit
#'
#' @description The function \code{\link{PLNMM}} produces an instance of this class.
#'
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for \code{\link[=PLNfamily_getBestModel]{getBestModel}},
#' \code{\link[=PLNfamily_getModel]{getModel}}, \code{\link[=plot.PLNfamily]{plot}}
#' and \code{\link[=predict.PLNfit]{predict}}.
#'
#' @field responses the matrix of responses common to every models
#' @field covariates the matrix of covariates common to every models
#' @field offsets the matrix of offsets common to every models
#' @field clusters the dimensions of the successively fitted models
#' @field models a list of \code{\link[=PLNMMfit]{PLNMMfit}} object, one for each number of clusters considered.
#' @field inception a \code{\link[=PLNfit]{PLNfit}} object, obtained when one single class is considered.
#' @field criteria a data frame with the value of some criteria (variational lower bound J, BIC, ICL and R2) for the different models.
#' @include PLNfamily-class.R
#' @importFrom R6 R6Class
#' @import ggplot2
#' @seealso The function \code{\link{PLNMM}}, the class \code{\link[=PLNMMfit]{PLNMMfit}}
PLNMMfamily <-
  R6Class(classname = "PLNMMfamily",
    inherit = PLNfamily,
     active = list(
      clusters = function() private$params
    )
)

PLNMMfamily$set("public", "initialize",
  function(clusters, responses, covariates, offsets, control) {

  ## initialize the required fields
  super$initialize(responses, covariates, offsets, control)
  private$params <- clusters

  if (control$trace > 0) cat("\n Perform kmeans on the inceptive model for initializing...")
  M     <- self$inception$var_par$M
  S     <- self$inception$var_par$S
  S_bar <- diag(colMeans(S))

  ## instantiate as many models as clusters
  self$models <- lapply(clusters, function(k){
    out_kmeans <- kmeans(M, centers = k)
    clusters <- out_kmeans$cl
    M_tilde <- M - out_kmeans$centers[clusters, ]
    tau <- matrix(0, private$n, k)
    tau[cbind(1:private$n, clusters)] <- 1

    model <- PLNMMfit$new(
      Theta = self$inception$model_par$Theta,
      Sigma = crossprod(M_tilde)/private$n + S_bar,
      mu    = out_kmeans$centers ,
      pi    = colMeans(tau),
      M     = M_tilde,
      S     = S,
      tau   = tau,
      J     = self$inception$loglik
    )
    return(model)
  })

})

PLNMMfamily$set("public", "optimize",
  function(control) {

  self$models <- mclapply(self$models, function(model) {
    ## initial parameters (model + variational)
    par0 <- c(
      model$model_par$Theta,
      model$model_par$mu,
      model$var_par$M,
      pmax(model$var_par$S,10*control$lbvar),
    )

    ## ===========================================
    ## OPTIMISATION
    if (control$trace == 1) {
      cat("\t Number of clusters =", model$k, "\r")
      flush.console()
    }

    if (control$trace > 1) {
      cat(" Number of clusters =", model$k)
      cat("\n\t conservative convex separable approximation for gradient descent")
    }

    ## CALL TO NLOPT OPTIMIZATION WITH BOX CONSTRAINT
    opts <- list(
      "algorithm"   = control$method,
      "maxeval"     = control$maxeval,
      "ftol_rel"    = control$ftol_rel,
      "ftol_abs"    = control$ftol_abs,
      "xtol_rel"    = control$xtol_rel,
      "xtol_abs"    = c(rep(0, private$p*(private$d + model$rank) + private$n * model$rank),
                        rep(control$xtol_abs, private$n*model$rank)),
      "lower_bound" = c(rep(-Inf, private$p*private$d), # Theta
                        rep(-Inf, private$p*model$k  ), # mu
                        rep(-Inf, private$n*private$p), # M
                        rep(control$lbvar,private$n*model$p)) # S
    )
    # optim.out <- optimization_PLNMM(par0, self$responses, self$covariates, self$offsets, model$k, opts)
    # optim.out$message <- statusToMessage(optim.out$status)

    ## ===========================================
    ## OUTPUT

    ### TODO
    ### COMPUTE Tau and pi from the optimized parameters

    ## formating parameters for output
    Theta <- matrix(optim.out$solution[1:(private$p*private$d)                                                        ], private$p, private$d)
    mu    <- matrix(optim.out$solution[private$p*private$d              + 1:(private$p*model$k)                    ], private$p, model$k)
    M     <- matrix(optim.out$solution[private$p*(private$d+model$k)    + 1:(private$n*private$p)                    ], private$n, model$rank)
    S     <- matrix(optim.out$solution[private$p*(private$d+model$rank) + private$n*model$rank + 1:(private$n*model$rank)], private$n, model$rank)
    Sigma <- B %*% (crossprod(M)/private$n + diag(colMeans(S), nrow = model$rank, ncol = model$rank)) %*% t(B)
    rownames(Theta) <- colnames(self$responses); colnames(Theta) <- colnames(self$covariates)
    rownames(B)     <- colnames(self$responses); colnames(B) <- 1:model$rank
    rownames(M)     <- rownames(self$responses); colnames(M) <- 1:model$rank
    rownames(Sigma) <- colnames(Sigma) <- colnames(self$responses)

    model$update(B = B, Theta = Theta, Sigma = Sigma, M = M, S = S, J = -optim.out$objective,
                 monitoring = list(objective = optim.out$objective,
                                   iterations = optim.out$iterations,
                                   status = optim.out$status,
                                   message = optim.out$message))
    return(model)
  }, mc.cores = control$cores, mc.allow.recursive = FALSE)
})

PLNMMfamily$set("public", "plot",
function(criteria = c("loglik", "BIC", "ICL")) {
  vlines <- sapply(criteria, function(crit) self$getBestModel(crit)$cluster)
  p <- super$plot(criteria) + xlab("# of clusters") + geom_vline(xintercept = vlines, linetype = "dashed", alpha = 0.25)
  p
})

PLNMMfamily$set("public", "show",
function() {
  super$show()
  cat(" Task: Mixture Model \n")
  cat("========================================================\n")
  cat(" - Number of clusters considered: from", min(self$clusters), "to", max(self$clusters),"\n")
  cat(" - Best model (regarding BIC): rank =", self$getBestModel("BIC")$cluster, "- R2 =", round(self$getBestModel("BIC")$R_squared, 2), "\n")
  cat(" - Best model (regarding ICL): rank =", self$getBestModel("ICL")$cluster, "- R2 =", round(self$getBestModel("ICL")$R_squared, 2), "\n")
})
PLNMMfamily$set("public", "print", function() self$show())

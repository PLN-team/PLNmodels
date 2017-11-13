#' An R6 Class to represent a collection of PLNPCAfit
#'
#' @description The function \code{\link{PLNPCA}} produces an instance of this class.
#'
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for \code{\link[=PLNPCAfamily_getBestModel]{getBestModel}},
#' \code{\link[=PLNPCAfamily_getModel]{getModel}} and \code{\link[=PLNPCAfamily_plot]{plot}}. Other methods
#'  should not be called as they are designed to be used during the optimization process.
#'
#' @field responses the matrix of responses common to every models
#' @field covariates the matrix of covariates common to every models
#' @field offsets the matrix of offsets common to every models
#' @field ranks the dimensions of the successively fitted models
#' @field models a list of \code{\link[=PLNPCAfit]{PLNPCAfit}} object, one per rank.
#' @field inception a \code{\link[=PLNfit-class]{PLNfit}} object, obtained when full rank is considered.
#' @field criteria a data frame with the value of some criteria (variational lower bound J, BIC, ICL and R2) for the different models.
#' @include PLNfamily-class.R
#' @importFrom R6 R6Class
#' @importFrom nloptr nloptr
#' @import ggplot2
#' @seealso The function \code{\link{PLNPCA}}, the class \code{\link[=PLNPCAfit-class]{PLNPCAfit}}
PLNPCAfamily <-
  R6Class(classname = "PLNPCAfamily",
    inherit = PLNfamily,
     active = list(
      ranks = function() self$params
    )
)

PLNPCAfamily$set("public", "initialize",
  function(ranks, responses, covariates, offsets, control) {

  ## initialize the required fields
  super$initialize(responses, covariates, offsets, control)
  self$params <- ranks

  if (control$trace > 0) cat("\n Perform SVD to initialize other parameters...")
  svdM     <- svd(self$inception$var_par$M, nu=max(ranks), nv=max(ranks))
  svdS     <- svd(self$inception$var_par$S, nu=max(ranks), nv=max(ranks))
  svdSigma <- svd(self$inception$model_par$Sigma, nu=max(ranks), nv=0)

  ## instantiate as many models as ranks
  self$models <- lapply(ranks, function(q){
    model <- PLNPCAfit$new(
      Theta = self$inception$model_par$Theta,
      Sigma = svdSigma$u[, 1:q, drop=FALSE] %*% diag(svdSigma$d[1:q],nrow=q, ncol=q) %*% t(svdSigma$u[, 1:q, drop=FALSE]),
      B = svdSigma$u[, 1:q, drop=FALSE] %*% sqrt(diag(svdSigma$d[1:q],nrow=q, ncol=q)),
      M = svdM$u[, 1:q, drop=FALSE] %*% diag(svdM$d[1:q], nrow=q, ncol=q) %*% t(svdM$v[1:q, 1:q, drop=FALSE]),
      S = svdS$u[, 1:q, drop=FALSE] %*% diag(svdS$d[1:q], nrow=q, ncol=q) %*% t(svdS$v[1:q, 1:q, drop=FALSE])
    )
    return(model)
  })

  ## declare the objective and gradient functions for optimization
  private$fn_optim <- fn_optim_PLNPCA_Cpp

})

PLNPCAfamily$set("public", "optimize",
  function(control) {

  KY <- sum(.logfactorial(self$responses))

  self$models <- mclapply(self$models, function(model) {
    ## initial parameters (model + variational)
    par0 <- c(model$model_par$Theta  , model$model_par$B,
              model$var_par$M, pmax(model$var_par$S,10*control$lbvar))

    ## ===========================================
    ## OPTIMISATION
    if (control$trace == 1) {
      cat("\t Rank approximation =",model$rank, "\r")
      flush.console()
    }

    if (control$trace > 1) {
      cat(" Rank approximation =",model$rank)
      cat("\n\t conservative convex separable approximation for gradient descent")
    }

    ## CALL TO NLOPT OPTIMIZATION WITH BOX CONSTRAINT
    lower.bound <- c(rep(-Inf, private$p*private$d)     , # Theta
                     rep(-Inf, private$p*model$rank) , # B
                     rep(-Inf, private$n*model$rank) , # M
                     rep(control$lbvar,private$n*model$rank)) # S
    opts <- list("algorithm"   = paste("NLOPT_LD",control$method, sep="_"),
                 "maxeval"     = control$maxeval,
                 "ftol_rel"    = control$ftol_rel,
                 "ftol_abs"    = control$ftol_abs,
                 "xtol_rel"    = control$xtol_rel,
                 "xtol_abs"    = control$xtol_abs,
                 "print_level" = max(0,control$trace-1))

    optim.out <- nloptr(par0, eval_f = private$fn_optim, lb = lower.bound, opts = opts,
                        q=model$rank, Y=self$responses, X=self$covariates, O=self$offsets, KY=KY)

    ## ===========================================
    ## OUTPUT

    ## formating parameters for output
    Theta <- matrix(optim.out$solution[1:(private$p*private$d)                                                     ], private$p, private$d)
    B     <- matrix(optim.out$solution[private$p*private$d              + 1:(private$p*model$rank)                    ], private$p, model$rank)
    M     <- matrix(optim.out$solution[private$p*(private$d+model$rank) + 1:(private$n*model$rank)                    ], private$n, model$rank)
    S     <- matrix(optim.out$solution[private$p*(private$d+model$rank) + private$n*model$rank + 1:(private$n*model$rank)], private$n, model$rank)
    Sigma <- B %*% (crossprod(M)/private$n + diag(colMeans(S), nrow = model$rank, ncol = model$rank)) %*% t(B)
    rownames(Theta) <- colnames(self$responses); colnames(Theta) <- colnames(self$covariates)
    rownames(B)     <- colnames(self$responses); colnames(B) <- 1:model$rank
    rownames(M)     <- rownames(self$responses); colnames(M) <- 1:model$rank
    rownames(Sigma) <- colnames(Sigma) <- colnames(self$responses)

    model$update(B = B, Theta = Theta, Sigma = Sigma, M = M, S = S,
                 monitoring = list(objective = optim.out$objective,
                                   iterations = optim.out$iterations,
                                   status = optim.out$status,
                                   message = optim.out$message))
    return(model)
  }, mc.cores = control$cores, mc.allow.recursive = FALSE)
})

PLNPCAfamily$set("public", "postTreatment",
function() {
  super$postTreatment()
  for (model in self$models) {
    model$setVisualization()
  }
})

#' A plot method for a collection of PLNPCAfit
#'
#' @name PLNPCAfamily_plot
#' @return Produces a plot  representing the evolution of the criteria of the different models considered,
#' highlighting the best model in terms of ICL. If not smooth, you may consider running PLNPCA again with
#' a smaller tolerance for convergence.
NULL
PLNPCAfamily$set("public", "plot",
function() {
  p <- super$plot() + xlab("rank")
  p <- p + annotate("text", x=self$params, y=min(self$criteria$loglik), angle=90, label=paste("R2 =", round(self$criteria$R_squared, 2)), size=3, alpha=0.7) +
    geom_vline(xintercept=self$getBestModel("ICL")$rank, linetype="dashed", alpha=0.5)
  p
})

PLNPCAfamily$set("public", "show",
function() {
  super$show()
  cat(" Task: Principal Component Analysis\n")
  cat("========================================================\n")
  cat(" - Ranks considered: from", min(self$params), "to", max(self$params),"\n")
  cat(" - Best model (regarding ICL): rank =", self$getBestModel("ICL")$rank, "- R2 =", round(self$getBestModel("ICL")$R_squared, 2), "\n")
})

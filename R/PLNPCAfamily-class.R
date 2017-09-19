#' An R6 Class to represent a collection of PLNPCAfit
#'
#' @description The function \code{\link{PLNPCA}} produces an instance of this class.
#'
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for \code{\link[=PLNPCAfamily_getBestModel]{getBestModel}},
#' \code{\link[=PLNPCAfamily_getModel]{getModel}} and \code{\link[=PLNPCAfamily_plot]{plot}}. Other methods
#'  should not be called as they are designed to be used during the optimization process.
#'
#' Fields should not be changed or manipulated by the user as they are updated internally
#' during the estimation process.
#'
#' @field ranks the dimensions of the successively fitted models
#' @field models a list of \code{\link[=PLNPCAfit]{PLNPCAfit}} object, one per rank.
#' @field inception a \code{\link[=PLNfit-class]{PLNfit}} object, obtained when full rank is considered.
#' @field criteria a data frame with the value of some criteria (variational lower bound J, BIC, ICL and R2) for the different models.
#' @field responses the matrix of responses common to every models
#' @field covariates the matrix of covariates common to every models
#' @field offsets the matrix of offsets common to every models
#' @field fn_optim the R functions used to compute the model's objective and gradient during the optimization process
#' @include PLNfamily-class.R
#' @importFrom R6 R6Class
#' @importFrom nloptr nloptr
#' @import ggplot2
#' @seealso The function \code{\link{PLNPCA}}, the class \code{\link[=PLNPCAfit-class]{PLNPCAfit}}
PLNPCAfamily <-
  R6Class(classname = "PLNPCAfamily",
    inherit = PLNfamily,
     public = list(
      ranks = "numeric"
    )
)

PLNPCAfamily$set("public", "initialize",
  function(ranks, responses, covariates, offsets, control) {

  ## initialize the required fields
  super$initialize(responses, covariates, offsets, control)
  self$ranks <- ranks

  ## instantiate as many models as ranks
  fit <- PLNPCAfit$new(model.par = self$inception$model.par, variational.par = list())

  svdSigma <- svd(self$inception$model.par$Sigma  , nu=max(ranks), nv=0)
  svdM     <- svd(self$inception$variational.par$M)
  svdS     <- svd(self$inception$variational.par$S)

  self$models <- lapply(ranks, function(q){
    model <- fit$clone()
    model$rank <- q
    model$model.par$B       <- svdSigma$u[, 1:q, drop=FALSE] %*% diag(sqrt(svdSigma$d[1:q]),nrow=q, ncol=q)
    model$variational.par$M <- svdM$u[, 1:q, drop=FALSE] %*% diag(svdM$d[1:q], nrow=q, ncol=q) %*% t(svdM$v[1:q, 1:q, drop=FALSE])
    model$variational.par$S <- svdS$u[, 1:q, drop=FALSE] %*% diag(svdS$d[1:q], nrow=q, ncol=q) %*% t(svdS$v[1:q, 1:q, drop=FALSE])
    return(model)
  })
  names(self$models) <- as.character(ranks)

  ## declare the objective and gradient functions for optimization
  self$fn_optim <- fn_optim_PLNPCA_Cpp

})

PLNPCAfamily$set("public", "optimize",
  function(control) {

  KY <-sum(.logfactorial(self$responses))

  self$models <- mclapply(self$models, function(model) {
    ## initial parameters (model + variational)
    par0 <- c(model$model.par$Theta  , model$model.par$B,
              model$variational.par$M, pmax(model$variational.par$S,10*control$lbvar))

    ## ===========================================
    ## OPTIMISATION
    if (control$trace > 0) cat("\n Rank approximation =",model$rank)
    if (control$trace > 1) cat("\n\t conservative convex separable approximation for gradient descent")

    ## CALL TO NLOPT OPTIMIZATION WITH BOX CONSTRAINT
    lower.bound <- c(rep(-Inf, self$p*self$d)     , # Theta
                     rep(-Inf, self$p*model$rank) , # B
                     rep(-Inf, self$n*model$rank) , # M
                     rep(control$lbvar,self$n*model$rank)) # S
    ## CALL TO NLOPT OPTIMIZATION WITH BOX CONSTRAINT
    opts <- list("algorithm"   = paste("NLOPT_LD",control$method, sep="_"),
                 "maxeval"     = control$maxeval,
                 "ftol_rel"    = control$ftol_rel,
                 "ftol_abs"    = control$ftol_abs,
                 "xtol_rel"    = control$xtol_rel,
                 "xtol_abs"    = control$xtol_abs,
                 "print_level" = max(0,control$trace-1))

    optim.out <- nloptr(par0, eval_f = self$fn_optim, lb = lower.bound, opts = opts,
                        q=model$rank, Y=self$responses, X=self$covariates, O=self$offsets, KY=KY)

    ## ===========================================
    ## OUTPUT

    ## formating parameters for output
    Theta <- matrix(optim.out$solution[1:(self$p*self$d)                                                     ], self$p, self$d)
    B     <- matrix(optim.out$solution[self$p*self$d              + 1:(self$p*model$rank)                    ], self$p, model$rank)
    M     <- matrix(optim.out$solution[self$p*(self$d+model$rank) + 1:(self$n*model$rank)                    ], self$n, model$rank)
    S     <- matrix(optim.out$solution[self$p*(self$d+model$rank) + self$n*model$rank + 1:(self$n*model$rank)], self$n, model$rank)
    Sigma <- B %*% (crossprod(M)/self$n + diag(colMeans(S), nrow = model$rank)) %*% t(B)
    rownames(Theta) <- colnames(self$responses); colnames(Theta) <- colnames(self$covariates)
    rownames(B)     <- colnames(self$responses); colnames(B) <- 1:model$rank
    rownames(M)     <- rownames(self$responses); colnames(M) <- 1:model$rank
    rownames(Sigma) <- colnames(Sigma) <- colnames(self$responses)

    ## compute some criteria for model selection
    J   <- -optim.out$objective
    BIC <- J - self$p * (self$d + model$rank) * log(self$n)
    ICL <- BIC - .5*self$n*model$rank*log(2*pi*exp(1)) - .5 * sum(log(S))

    model$model.par       <- list(B = B, Theta = Theta, Sigma = Sigma)
    model$variational.par <- list(M = M, S = S)
    model$criteria        <- c(J = J, BIC = BIC, ICL = ICL)
    model$convergence     <- data.frame(status = optim.out$message, objective = optim.out$objective, iterations=optim.out$iterations)
    return(model)
  }, mc.cores = control$cores, mc.allow.recursive = FALSE)
})

PLNPCAfamily$set("public", "postTreatment",
function() {
  self$computeR2()
  self$setCriteria()
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
  p <- p + annotate("text", x=self$ranks, y=min(self$criteria$J), angle=90, label=paste("R2 =", round(self$criteria$R2, 2)), size=3, alpha=0.7) +
    geom_vline(xintercept=self$getBestModel("ICL")$rank, linetype="dashed", alpha=0.5)
  p
})

PLNPCAfamily$set("public", "show",
function() {
  super$show()
  cat(" Task: Principal Component Analysis\n")
  cat("======================================================\n")
  cat(" - Ranks considered: from", min(self$ranks), "to", max(self$ranks),"\n")
  cat(" - Best model (regardings ICL): rank =", self$getBestModel("ICL")$rank, "- R2 =", round(self$getBestModel("ICL")$criteria["R2"], 2), "\n")
})

# PLNfamily$set("public", "print", function() self$show())

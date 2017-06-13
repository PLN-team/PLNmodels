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
  KY <-sum(.logfactorial(self$responses))

  ## instantiate as many models as ranks
  fit <- PLNPCAfit$new(model.par       = self$inception$model.par,
                       variational.par = self$inception$variational.par)
  svdSigma <- svd(self$inception$model.par$Sigma, nu=max(ranks), nv=0)
  self$models <- lapply(ranks, function(q){
    model <- fit$clone()
    model$rank <- q
    model$model.par$B <- svdSigma$u[, 1:q, drop=FALSE] %*% diag(sqrt(svdSigma$d[1:q]),nrow=q, ncol=q)
    return(model)
  })
  names(self$models) <- as.character(ranks)

  ## declare the objective and gradient functions for optimization
  self$fn_optim <- function(par,q) {

    Theta <- matrix(par[1:(self$p*self$d)]                          , self$p,self$d)
    B     <- matrix(par[self$p*self$d                + 1:(self$p*q)], self$p,q)
    M     <- matrix(par[self$p*(self$d+q)            + 1:(self$n*q)], self$n,q)
    S     <- matrix(par[self$p*(self$d+q)+(self$n*q) + 1:(self$n*q)], self$n,q)

    Z     <- tcrossprod(M,B) + tcrossprod(self$covariates, Theta) + self$offsets
    A <- exp (.trunc(Z + .5*tcrossprod(S^2, B^2)))

    gr.Theta <- crossprod(A-self$responses, self$covariates)
    gr.B     <- crossprod(A-self$responses, M) + crossprod(A,S^2) * B
    gr.M     <- (A-self$responses) %*% B + M
    gr.S     <- S - 1/S  + (A %*% (B^2)) * S

    return(list(
      "objective" = sum(as.numeric(A -self$responses*Z)) +.5*sum(M^2 + S^2 - 2*log(S)-1) + KY,
      "gradient"  = c(gr.Theta,gr.B,gr.M,gr.S)
    ))
  }

  # self$objective <- function(par,q,KY) {
  #   Theta <- matrix(par[1:(self$p*self$d)]                          , self$p,self$d)
  #   B     <- matrix(par[self$p*self$d                + 1:(self$p*q)], self$p,q)
  #   M     <- matrix(par[self$p*(self$d+q)            + 1:(self$n*q)], self$n,q)
  #   S     <- matrix(par[self$p*(self$d+q)+(self$n*q) + 1:(self$n*q)], self$n,q)
  #   Z     <- tcrossprod(M,B) + tcrossprod(self$covariates, Theta) + self$offsets
  #   return(sum(as.numeric(exp(.trunc(Z +.5*tcrossprod(S^2, B^2)))-self$responses*Z)) +.5*sum(M^2 + S^2 - 2*log(S)-1) + KY)
  # }
  #
  # self$gradient <- function(par,q,KY) {
  #   Theta <- matrix(par[1:(self$p*self$d)]                          , self$p,self$d)
  #   B     <- matrix(par[self$p*self$d                + 1:(self$p*q)], self$p,q)
  #   M     <- matrix(par[self$p*(self$d+q)            + 1:(self$n*q)], self$n,q)
  #   S     <- matrix(par[self$p*(self$d+q)+(self$n*q) + 1:(self$n*q)], self$n,q)
  #   A <- exp (.trunc(tcrossprod(M,B) + tcrossprod(self$covariates,Theta) + self$offsets + .5*tcrossprod(S^2, B^2)))
  #   gr.Theta <- crossprod(A-self$responses, self$covariates)
  #   gr.B     <- crossprod(A-self$responses, M) + crossprod(A,S^2) * B
  #   gr.M     <- (A-self$responses) %*% B + M
  #   gr.S     <- S - 1/S  + (A %*% (B^2)) * S
  #   return(c(gr.Theta,gr.B,gr.M,gr.S))
  # }
})

PLNPCAfamily$set("public", "optimize",
  function(control) {

  self$models <- mclapply(self$models, function(model) {
    ## initial parameters (model + variational)
    par0 <- c(model$model.par$Theta  , model$model.par$B,
              #model$variational.par$M, model$variational.par$S)
              rep(0, self$n*model$rank) , # M
              rep(10*control$lbvar,self$n*model$rank))
              # model$variational.par$M, model$variational.par$S)
    ## ===========================================
    ## OPTIMISATION
    if (control$trace > 0) cat("\n Rank approximation =",model$rank)
    if (control$trace > 1) cat("\n\t conservative convex separable approximation for gradient descent")

    ## CALL TO L-BFGS OPTIMIZATION WITH BOX CONSTRAINT
    lower.bound <- c(rep(-Inf, self$p*self$d)     , # Theta
                     rep(-Inf, self$p*model$rank) , # B
                     rep(-Inf, self$n*model$rank) , # M
                     rep(control$lbvar,self$n*model$rank)) # S
    ## CALL TO NLOPT OPTIMIZATION WITH BOX CONSTRAINT
    opts <- list("algorithm" = "NLOPT_LD_TNEWTON_PRECOND",
                 "maxeval"   = control$maxit,
                 "xtol_rel"  = control$xtol,
                 "ftol_rel"  = control$ftol,
                 "print_level" = max(0, control$trace-1))
    optim.out <- nloptr(par0, eval_f = self$fn_optim, lb = lower.bound, opts = opts, q=model$rank)

    ## ===========================================
    ## OUTPUT
    ## formating parameters for output
    Theta <- matrix(optim.out$solution[1:(self$p*self$d)]           , self$p,self$d)
    B     <- matrix(optim.out$solution[self$p*self$d      + 1:(self$p*model$rank)], self$p,model$rank)
    M     <- matrix(optim.out$solution[self$p*(self$d+model$rank)  + 1:(self$n*model$rank)], self$n,model$rank)
    S     <- matrix(optim.out$solution[self$p*(self$d+model$rank)+self$n*model$rank + 1:(self$n*model$rank)],self$n,model$rank)
    rownames(Theta) <- colnames(self$responses); colnames(Theta) <- colnames(self$covariates)
    rownames(B)     <- colnames(self$responses); colnames(B) <- 1:model$rank
    rownames(M)     <- rownames(self$responses); colnames(M) <- 1:model$rank

    ## compute some criteria for evaluation
    J   <- -optim.out$objective
    BIC <- J - self$p * (self$d + model$rank) * log(self$n)
    ICL <- BIC - .5*self$n*model$rank *log(2*pi*exp(1)) - sum(log(S))
    loglik <- sapply(1:model$rank, function(q_) {
      Z <- tcrossprod(M[,1:q_,drop=FALSE],B[,1:q_,drop=FALSE]) + tcrossprod(self$covariates, Theta) + self$offsets
      return(sum(self$responses * (Z)) - sum(as.numeric(self$responses)))
    })
    lmin <- sum(self$responses * (self$offsets + tcrossprod(self$covariates, Theta))) - sum(as.numeric(self$responses))
    lmax <- sum(self$responses * (log(self$responses + 1*(self$responses == 0)) - 1))
    R2  <- (loglik[model$rank] - lmin)  / (lmax - lmin)

    model$model.par       <- list(B = B, Theta = Theta, Sigma = tcrossprod(B))
    model$variational.par <- list(M = M, S = S)
    model$criteria        <- c(J = J, R2 = R2, BIC = BIC, ICL = ICL)
    model$convergence     <- data.frame(status = optim.out$message, objective = optim.out$objective, iterations=optim.out$iterations)
    return(model)
  }, mc.cores = control$cores, mc.allow.recursive = FALSE)
})

PLNPCAfamily$set("public", "postTreatment",
function() {
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

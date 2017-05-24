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
#' @field criteria a data frame with the value of some criteria (variational lower bound J, BIC, ICL and R2) for the different models.
#' @field responses the matrix of responses common to every models
#' @field covariates the matrix of covariates common to every models
#' @field offsets the matrix of offsets common to every models
#' @field objective the R function used to compute the model's objective during the optimization process
#' @field gradient the R function to compute the model's gradient during the optimization process
#' @include PLNfamily-class.R
#' @importFrom R6 R6Class
#' @seealso The function \code{\link{PLNPCA}}, the class \code{\link[=PLNPCAfit]{PLNPCAfit}}
PLNPCAfamily <-
  R6Class(classname = "PLNPCAfamily",
    inherit = PLNfamily,
     public = list(
      ranks = "numeric"
    )
)

PLNPCAfamily$set("public", "initialize",
  function(ranks, responses, covariates, offsets) {

  ## initialize the required fields
  super$initialize(responses, covariates, offsets)
  self$ranks <- ranks

  ## recover the initial model for each rank with glm Poisson models
  fit <- PLNPCAfit$new()
  svdSigma <- svd(self$init.par$Sigma, nu=max(ranks), nv=0)
  self$models <- lapply(ranks, function(q){
    B <- svdSigma$u[, 1:q, drop=FALSE] %*% diag(sqrt(svdSigma$d[1:q]),nrow=q, ncol=q)
    model <- fit$clone()
    model$model.par <- list(B = B, Theta = self$init.par$Theta)
    model$rank <- q
    return(model)
  })
  names(self$models) <- as.character(ranks)

  ## declare the objective and gradient functions for optimization

  self$objective <- function(par,n,p,q,d,KY) {
    Theta <- matrix(par[1:(p*d)]                , p,d)
    B     <- matrix(par[p*d           + 1:(p*q)], p,q)
    M     <- matrix(par[p*(d+q)       + 1:(n*q)], n,q)
    S     <- matrix(par[p*(d+q)+(n*q) + 1:(n*q)], n,q)
    Z     <- tcrossprod(M,B) + tcrossprod(self$covariates, Theta) + self$offsets
    return(sum(as.numeric(exp(.trunc(Z +.5*tcrossprod(S^2, B^2)))-self$responses*Z)) +.5*sum(M^2 + S^2 - 2*log(S)-1) + KY)
  }

  self$gradient <- function(par,n,p,q,d,KY) {
    Theta <- matrix(par[1:(p*d)]                , p,d)
    B     <- matrix(par[p*d           + 1:(p*q)], p,q)
    M     <- matrix(par[p*(d+q)       + 1:(n*q)], n,q)
    S     <- matrix(par[p*(d+q)+(n*q) + 1:(n*q)], n,q)
    A <- exp (.trunc(tcrossprod(M,B) + tcrossprod(self$covariates,Theta) + self$offsets + .5*tcrossprod(S^2, B^2)))
    gr.Theta <- crossprod(A-self$responses, self$covariates)
    gr.B     <- crossprod(A-self$responses, M) + crossprod(A,S^2) * B
    gr.M     <- (A-self$responses) %*% B + M
    gr.S     <- S - 1/S  + (A %*% (B^2)) * S
    return(c(gr.Theta,gr.B,gr.M,gr.S))
  }
})

PLNPCAfamily$set("public", "optimize",
  function(control) {

  ## ===========================================
  ## INITIALISATION
  n  <- nrow(self$responses); p <- ncol(self$responses); d <- ncol(self$covariates) # problem dimension
  KY <- sum(.logfactorial(self$responses)) ## constant quantity in the objective
  control.lbfgsb <- list(maxit = control$maxit, pgtol=control$pgtol, factr=control$factr,trace=ifelse(control$trace>1,control$trace,0))

  self$models <- mclapply(self$models, function(model) {
    ## initial parameters (model + variational)
    q <- model$rank
    par0 <- c(model$model.par$Theta, model$model.par$B, matrix(0, n, q), matrix(control$lbvar*10,n,q))

    ## ===========================================
    ## OPTIMISATION
    if (control$trace > 0) cat("\n Rank approximation =",q)
    if (control$trace > 1) cat("\n\t L-BFGS-B for gradient descent...")

    ## CALL TO L-BFGS OPTIMIZATION WITH BOX CONSTRAINT
    lower.bound <- c(rep(-Inf, p*d), ## Theta
                     rep(-Inf, p*q), ## B
                     rep(-Inf, n*q), ## M
                     rep(control$lbvar,n*q)) # S
    optim.out <- optim(par0, fn=self$objective, gr=self$gradient, lower=lower.bound, control = control.lbfgsb, method="L-BFGS-B",n,p,q,d,KY)
    ## C++ implementation: TOO FRAGILE for the moment...
    ##    par <- solvePLN(q, c(Theta,B,M,S), X, Y, O, min.var*10, max.iter, 1e-2, 0, 1e-2)$par

    ## ===========================================
    ## OUTPUT
    ## formating parameters for output
    Theta <- matrix(optim.out$par[1:(p*d)]           , p,d)
    B     <- matrix(optim.out$par[p*d      + 1:(p*q)], p,q)
    M     <- matrix(optim.out$par[p*(d+q)  + 1:(n*q)], n,q)
    S     <- matrix(optim.out$par[p*(d+q)+n*q + 1:(n*q)],n,q)
    rownames(Theta) <- colnames(self$responses); colnames(Theta) <- colnames(self$covariates)
    rownames(B)     <- colnames(self$responses); colnames(B) <- 1:q
    rownames(M)     <- rownames(self$responses); colnames(M) <- 1:q

    ## compute some criteria for evaluation
    loglik <- sapply(1:q, function(q_) {
      Z <- tcrossprod(M[,1:q_,drop=FALSE],B[,1:q_,drop=FALSE]) + tcrossprod(self$covariates, Theta) + self$offsets
      return(sum(self$responses * (Z)) - sum(as.numeric(self$responses)) - KY)
    })
    J   <- -optim.out$value
    lmin <- sum(self$responses * (self$offsets + tcrossprod(self$covariates, model$model.par$Theta))) - sum(as.numeric(self$responses)) - KY
    lmax <- sum(self$responses * (log(self$responses + 1*(self$responses == 0)) - 1)) - KY
    R2  <- (loglik[q] - lmin)  / (lmax - lmin)
    BIC <- J - p * (d + q) * log(n)
    ICL <- BIC - .5*n*q *log(2*pi*exp(1)) - sum(log(S))

    model$model.par       <- list(B = B, Theta = Theta)
    model$variational.par <- list(M = M, S = S)
    model$criteria        <- c(J = J, R2 = R2, BIC = BIC, ICL = ICL, lmin = lmin, lmax = lmax, lq = loglik[q])
    model$loglik          <- loglik
    model$convergence     <- optim.out$convergence
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

#' An R6 Class to represent a collection of PLNnetworkfit
#'
#' @description The function \code{\link{PLNnetwork}} produces an instance of this class.
#'
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for \code{\link[=PLNnetworkfamily_getBestModel]{getBestModel}},
#' \code{\link[=PLNnetworkfamily_getModel]{getModel}} and \code{\link[=PLNnetworkfamily_plot]{plot}}. Other methods
#'  should not be called as they are designed to be used during the optimization process.
#'
#' Fields should not be changed or manipulated by the user as they are updated internally
#' during the estimation process.
#'
#' @field type a character indicating the model used for the covariance matrix in the variational Gaussian approximation. Either "diagonal" or "spherical".
#' @field penatlies the sparsity level of the network in the successively fitted models
#' @field models a list of \code{\link[=PLNnetworkfit-class]{PLNnetworkfit}} object, one per rank.
#' @field criteria a data frame with the value of some criteria (variational lower bound J, BIC, ICL and R2) for the different models.
#' @field responses the matrix of responses common to every models
#' @field covariates the matrix of covariates common to every models
#' @field offsets the matrix of offsets common to every models
#' @field objective the R function used to compute the model's objective during the optimization process
#' @field gradient the R function to compute the model's gradient during the optimization process
#' @include PLNfamily-class.R
#' @importFrom R6 R6Class
#' @seealso The function \code{\link{PLNnetwork}}, the class \code{\link[=PLNnetworkfit-class]{PLNnetworkfit}}
PLNnetworkfamily <-
  R6Class(classname = "PLNnetworkfamily",
    inherit = PLNfamily,
     public = list(
      penalties = "numeric"
    )
)

PLNnetworkfamily$set("public", "initialize",
  function(type, penalties, responses, covariates, offsets) {

  ## initialize the required fields
  super$initialize(type, responses, covariates, offsets)
  stopifnot(all(penalties >= 0))
  self$penalties <- penalties

  ## recover the initial model for each rank with glm Poisson models
  fit <- PLNnetworkfit$new(type = self$type)
  self$models <- lapply(penalties, function(penalty){
    # glasso.out <- glasso::glasso(self$init.par$Sigma, rho = penalty, penalize.diagonal = FALSE, trace=TRUE, approx=TRUE)
    # Omega <- glasso.out$wi
    # Sigma <- glasso.out$w
    model <- fit$clone()
    model$model.par <- list(Sigma = self$init.par$Sigma, Theta = self$init.par$Theta)
    model$penalty <- penalty
    return(model)
  })

  ## declare the objective and gradient functions for optimization
  self$objective <- function(par,n,p,d,KY,logDetOmega,Omega) {
      Theta <- matrix(par[1:(d*p)]             , d,p)
      M     <- matrix(par[d*p        + 1:(n*d)], n,d)
      S     <- matrix(par[(n+p)*d    + 1:(n*d)], n,d)

      Z <- self$offsets + tcrossprod(self$covariates, Theta) + M
      logP.Z  <- n/2 * (logDetOmega - sum(diag(Omega)*colMeans(S))) - .5*sum(diag(Omega %*% crossprod(M)))

      # we return - J
      return(sum(as.numeric(exp(.trunc(Z + .5*S)) - self$responses*Z)) - logP.Z - .5*sum(log(S)+1)+ KY)
  }

  self$gradient <- function(par,n,p,d,KY,logDetOmega,Omega) {
      Theta <- matrix(par[1:(p*d)]             , d,p)
      M     <- matrix(par[p*d        + 1:(n*d)], n,d)
      S     <- matrix(par[(n+p)*d    + 1:(n*d)], n,d)

      A <- exp (.trunc(self$offsets + tcrossprod(self$covariates, Theta) + M + .5*S))
      gr.Theta <- crossprod(self$covariates, A - self$responses)
      gr.M  <- M %*% Omega + A - self$responses
      gr.S  <- .5 * (matrix(rep(diag(Omega),n), n,d, byrow = TRUE) + A - 1/S)

      return(c(gr.Theta,gr.M,gr.S))
    }

})

PLNnetworkfamily$set("public", "optimize",
  function(control) {

  ## ===========================================
  ## INITIALISATION
  n  <- nrow(self$responses); d <- ncol(self$responses); p <- ncol(self$covariates) # problem dimension
  KY <- sum(.logfactorial(self$responses)) ## constant quantity in the objective
  control.lbfgsb <- list(maxit = control$maxit, pgtol=control$pgtol, factr=control$factr,trace=ifelse(control$trace>1,control$trace,0))
  control.MM     <- list(maxit = control$MMmaxit, tol=control$MMtol, trace=control$trace)

  self$models <- mclapply(self$models, function(model) {

    penalty <- model$penalty

    ## initial parameters (model + variational)
    Omega <- solve(model$model.par$Sigma + diag(rep(1e-6,d)))
    logDetOmega <- determinant(Omega, logarithm=TRUE)$modulus
    par0 <- c(model$model.par$Theta, matrix(0,n,d), matrix(control$lbvar*10,n,d))
    lower.bound <- c(rep(-Inf, p*d), ## Theta
                     rep(-Inf, n*d), ## M
                     rep(control$lbvar, n*d)) # S

    ## ===========================================
    ## OPTIMISATION
    if (control$trace > 0) cat("\n sparsifying penalty =",penalty)
    if (control$trace > 1) cat("\n\t L-BFGS-B for gradient descent")
    if (control$trace > 1) cat("\n\t graphical-Lasso for sparse covariance estimation")
    # browser()

    cond <- FALSE; iter <- 0
    convergence <- numeric(control.MM$maxit)
    objective <- numeric(control.MM$maxit)
    if (control$trace > 0) cat("\n\titeration: ")
    while(!cond) {
      iter <- iter + 1
      if (control$trace > 0) cat("",iter)

      ## CALL TO L-BFGS OPTIMIZATION WITH BOX CONSTRAINT
      optim.out <- optim(par0, fn=self$objective, gr=self$gradient, lower=lower.bound, control = control.lbfgsb, method="L-BFGS-B",n,p,d,KY,logDetOmega,Omega)
      objective[iter] <- optim.out$value
      convergence[iter] <- sqrt(sum((optim.out$par - par0)^2)/sum(par0^2))
      if ((convergence[iter] < control.MM$tol) | (iter > control.MM$maxit)) {
        cond <- TRUE
      } else {
        M <- matrix(optim.out$par[p*d        + 1:(n*d)], n,d)
        S <- matrix(optim.out$par[(n+p)*d    + 1:(n*d)], n,d)
        Sigma <- crossprod(M)/n + diag(colMeans(S))
        if (penalty <=1e-5) {
          Omega <- solve(Sigma)
        } else {
          Omega <- glasso::glasso(Sigma, rho=penalty, penalize.diagonal = FALSE)$wi
        }
        logDetOmega <- determinant(Omega, logarithm=TRUE)$modulus
        par0 <- optim.out$par
      }
    }
    ## ===========================================
    ## OUTPUT
    ## formating parameters for output
    Theta <- matrix(optim.out$par[1:(p*d)] ,d,p)
    rownames(Theta) <- colnames(self$responses); colnames(Theta) <- colnames(self$covariates)
    dimnames(S)     <- dimnames(self$responses)
    dimnames(M)     <- dimnames(self$responses)

    ## compute some criteria for evaluation
    ## loglik <- sum(self$responses * (Z)) - sum(as.numeric(self$responses)) - KY
    J   <- -optim.out$value
    # lmin <- sum(self$responses * (self$offsets + tcrossprod(self$covariates, model$model.par$Theta))) - sum(as.numeric(self$responses)) - KY
    # lmax <- sum(self$responses * (log(self$responses + 1*(self$responses == 0)) - 1)) - KY
    # R2  <- (loglik[q] - lmin)  / (lmax - lmin)
    ## BIC <- J - p * (d + q) * log(n)
    ## ICL <- BIC - .5*n*q *log(2*pi*exp(1)) - sum(log(S)) * switch(self$type, "diagonal"=1, "spherical"=q)

    model$model.par       <- list(Omega = Omega, Sigma = solve(Sigma), Theta = Theta)
    model$variational.par <- list(M = M, S = S)
    model$criteria        <- c(J = J)
    model$loglik          <- NA #loglik
    model$convergence     <- data.frame(convergence = convergence[1:iter], objective = objective[1:iter])
    return(model)
  }, mc.cores=control$cores)
})

#' A Reference Class to represent a collection of PLNfit
#'
#' @description The function \code{\link{PLNPCA}} produces an instance of this class.
#'
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for \code{\link[=PLNfamily_getBestModel]{getBestModel}},
#' \code{\link[=PLNfamily_getModel]{getModel}} and \code{\link[=PLNfamily_plot]{plot}}. Other methods
#'  should not be called as they are designed to be used during the optimization process.
#'
#' Fields should not be changed or manipulated by the user as they are updated internally
#' during the estimation process.
#'
#' @field type a character indicating the model used for the covariance matrix in the variational Gaussian approximation. Either "diagonal" or "spherical".
#' @field ranks the dimensions of the successively fitted models
#' @field models a list of \code{\link[=PLNfit.PCA-class]{PLNfit}} object, one per rank.
#' @field criteria a data frame with the value of some criteria (variational lower bound J, BIC, ICL and R2) for the different models.
#' @field responses the matrix of responses common to every models
#' @field covariates the matrix of covariates common to every models
#' @field offsets the matrix of offsets common to every models
#' @field objective the R function used to compute the model's objective during the optimization process
#' @field gradient the R function to compute the model's gradient during the optimization process
#' @seealso The function \code{\link{PLNPCA}}, the class \code{\link[=PLNfit.PCA-class]{PLNfit.PCA}}
PLNfamily <-
  setRefClass(Class="PLNfamily",
    fields = list(
      type       = "character",
      ranks      = "numeric",
      models     = "list",
      criteria   = "data.frame",
      responses  = "matrix",
      covariates = "matrix",
      offsets    = "matrix",
      objective  = "function",
      gradient   = "function"
    )
)

PLNfamily$methods(initialize = function(type,ranks,responses,covariates,offsets) {
  ## initialize the required fields
  type       <<- type
  ranks      <<- ranks
  responses  <<- responses
  covariates <<- covariates
  offsets    <<- offsets

  ## set names of the data matrices
  if (is.null(rownames(responses)))  rownames(responses)  <<- 1:nrow(responses)
  if (is.null(colnames(responses)))  colnames(responses)  <<- 1:ncol(responses)
  if (is.null(rownames(covariates))) rownames(covariates) <<- 1:nrow(covariates)
  if (is.null(colnames(covariates))) colnames(covariates) <<- 1:ncol(covariates)

  ## recover the initial model for each rank with glm Poisson models
  glmP  <- lapply(1:ncol(responses), function(j) glm.fit(covariates, responses[, j], offset = offsets[,j], family = poisson()))
  Theta <- matrix(Reduce(rbind, lapply(glmP, coefficients)), ncol=ncol(covariates))
  Sigma <- cov(sapply(glmP, residuals.glm, type="pearson"))
  svdSigma <- svd(Sigma, nu=max(ranks), nv=0)

  models <<- lapply(ranks, function(q){
    B <- svdSigma$u[, 1:q, drop=FALSE] %*% diag(sqrt(svdSigma$d[1:q]),nrow=q, ncol=q)
    return(PLNfit$new(type = "poisson-GLM", rank  = q, model.par = list(B = B, Theta = Theta))
    )
  })

  ## declare the objective and gradient functions for optimization
  objective <<- switch(type,
      "diagonal" = function(par,n,p,q,d,KY) {
        Theta <- matrix(par[1:(p*d)]                , p,d)
        B     <- matrix(par[p*d           + 1:(p*q)], p,q)
        M     <- matrix(par[p*(d+q)       + 1:(n*q)], n,q)
        S     <- matrix(par[p*(d+q)+(n*q) + 1:(n*q)], n,q)
        Z     <- tcrossprod(M,B) + tcrossprod(covariates, Theta) + offsets
      return(sum(as.numeric(exp(.trunc(Z +.5*tcrossprod(S^2, B^2)))-responses*Z)) +.5*sum(M^2 + S^2 - 2*log(S)-1) + KY)
    },
      "spherical" = function(par,n,p,q,d,KY) {
        Theta <- matrix(par[1:(p*d)]           , p,d)
        B     <- matrix(par[p*d      + 1:(p*q)], p,q)
        M     <- matrix(par[p*(d+q)  + 1:(n*q)], n,q)
        S     <- c(par[p*(q+d) + (n*q) + 1:n])
        Z <- tcrossprod(M,B) + tcrossprod(covariates, Theta) + offsets
      return(sum(as.numeric(exp(.trunc(Z + .5*outer(S^2, rowSums(B^2))))-responses*Z)) + .5*sum(M^2) - .5*q*sum(2*log(S)-S^2+1) + KY)
    })

  gradient <<- switch(type,
        "diagonal" = function(par,n,p,q,d,KY) {
          Theta <- matrix(par[1:(p*d)]                , p,d)
          B     <- matrix(par[p*d           + 1:(p*q)], p,q)
          M     <- matrix(par[p*(d+q)       + 1:(n*q)], n,q)
          S     <- matrix(par[p*(d+q)+(n*q) + 1:(n*q)], n,q)
          A <- exp (.trunc(tcrossprod(M,B) + tcrossprod(covariates,Theta) + offsets + .5*tcrossprod(S^2, B^2)))
          gr.Theta <- crossprod(A-responses, covariates)
          gr.B     <- crossprod(A-responses, M) + crossprod(A,S^2) * B
          gr.M     <- (A-responses) %*% B + M
          gr.S     <- S - 1/S  + (A %*% (B^2)) * S
        return(c(gr.Theta,gr.B,gr.M,gr.S))
      },
        "spherical" = function(par,n,p,q,d,KY) {
          Theta <- matrix(par[1:(p*d)]           , p,d)
          B     <- matrix(par[p*d      + 1:(p*q)], p,q)
          M     <- matrix(par[p*(d+q)  + 1:(n*q)], n,q)
          S     <- c(par[p*(q+d) + (n*q) + 1:n])
          A <- exp (.trunc(tcrossprod(M,B) + tcrossprod(covariates, Theta) + offsets + .5 * outer(S^2, rowSums(B^2))))
          gr.M  <- (A-responses) %*% B + M
          gr.S  <- q*(S - 1/S) + rowSums(A %*% B^2 * S)
          gr.Theta <- crossprod(A-responses,covariates)
          ## using R recycling trick
          gr.B <-  crossprod(A-responses,M) + B * colSums(A * S^2)
        return(c(gr.Theta,gr.B,gr.M,gr.S))
      })

})

PLNfamily$methods(optimize = function(control) {

  ## ===========================================
  ## INITIALISATION
  n  <- nrow(responses); p <- ncol(responses); d <- ncol(covariates) # problem dimension
  KY <- sum(.logfactorial(responses)) ## constant quantity in the objective
  control.lbfgsb <- list(maxit = control$maxit, pgtol=control$pgtol, factr=control$factr,trace=ifelse(control$trace>1,control$trace,0))

  models <<- mclapply(models, function(model) {
    ## initial parameters (model + variational)
    q <- model$rank
    par0 <- c(model$model.par$Theta, model$model.par$B, matrix(0, n, q),
              switch(type, "diagonal" = matrix(control$lbvar*10,n,q), "spherical" = rep(control$lbvar*10,n)))

    ## ===========================================
    ## OPTIMISATION
    if (control$trace > 0) cat("\n Rank approximation =",q)
    if (control$trace > 1) cat("\n\t L-BFGS-B for gradient descent...")

    ## CALL TO L-BFGS OPTIMIZATION WITH BOX CONSTRAINT
    lower.bound <- c(rep(-Inf, p*d), ## Theta
                     rep(-Inf, p*q), ## B
                     rep(-Inf, n*q), ## M
                     rep(control$lbvar, switch(type, "diagonal" = n*q, "spherical" = n))) # S
    optim.out <- optim(par0, fn=objective, gr=gradient, lower=lower.bound, control = control.lbfgsb, method="L-BFGS-B",n,p,q,d,KY)
    ## C++ implementation: TOO FRAGILE for the moment...
    ##    par <- solvePLN(q, c(Theta,B,M,S), X, Y, O, min.var*10, max.iter, 1e-2, 0, 1e-2)$par

    ## ===========================================
    ## OUTPUT

    ## formating parameters for output
    Theta <- matrix(optim.out$par[1:(p*d)]           , p,d); colnames(Theta) <- colnames(covariates)
    B     <- matrix(optim.out$par[p*d      + 1:(p*q)], p,q); rownames(B) <- colnames(responses)
    M     <- matrix(optim.out$par[p*(d+q)  + 1:(n*q)], n,q); rownames(M) <- rownames(responses)
    S     <- switch(type,
                    "diagonal"  = matrix(optim.out$par[p*(d+q)+n*q + 1:(n*q)],n,q),
                    "spherical" = c(optim.out$par[p*(d+q)+n*q + 1:n]))
    rownames(Theta) <- colnames(responses); colnames(Theta) <- colnames(covariates)
    rownames(B)     <- colnames(responses); colnames(B) <- 1:q
    rownames(M)     <- rownames(responses); colnames(M) <- 1:q


    ## compute some criteria for evaluation
    loglik <- sapply(1:q, function(q_) {
      Z <- tcrossprod(M[,1:q_,drop=FALSE],B[,1:q_,drop=FALSE]) + tcrossprod(covariates, Theta) + offsets
      return(sum(responses * (Z)) - sum(as.numeric(responses)) - KY)
    })
    J   <- -optim.out$value # objective(optim.out$par,n,p,q,d,KY)
    lmin <- sum(responses * (offsets + tcrossprod(covariates, model$model.par$Theta))) - sum(as.numeric(responses)) - KY
    lmax <- sum(responses * (log(responses + 1*(responses == 0)) - 1)) - KY
    R2  <- (loglik[q] - lmin)  / (lmax - lmin)
    BIC <- J - p * (d + q) * log(n)
    ICL <- BIC - .5*n*q *log(2*pi*exp(1)) - sum(log(S)) * switch(type, "diagonal"=1, "spherical"=q)

    return(PLNfit$new(
      type            = type,
      rank            = q,
      model.par       = list(B = B, Theta = Theta),
      variational.par = list(M = M, S = S),
      criteria        = c(J = J, R2 = R2, BIC = BIC, ICL = ICL, lmin = lmin, lmax = lmax, lq = loglik[q]),
      loglik          = loglik,
      convergence     = optim.out$convergence))
  }, mc.cores = control$cores, mc.allow.recursive = FALSE)
})

#' Best model extraction from a collection of PLNfit
#'
#' @name PLNfamily_getBestModel
#'
#' @param crit a character for the criterion used to performed the selection. Either
#' "ICL", "BIC", "J" or "R2". Default is "ICL.
#' @return  Send back a object with class \code{\link[=PLNfit.PCA-class]{PLNfit}}.
NULL
PLNfamily$methods(getBestModel = function(crit=c("ICL", "BIC", "J", "R2")){
  crit <- match.arg(crit)
  if(length(criteria$BIC) >1) {
    id <- switch(crit,
    "BIC" = which.max(criteria$BIC),
    "ICL" = which.max(criteria$ICL),
    "J"   = which.max(criteria$J),
    "R2"  = which.max(criteria$R2))
  } else {id <- 1}
    model <- models[[id]]$copy()
    return(model)
})

#' Model extraction from a collection of PLNfit
#'
#' @name PLNfamily_getModel
#'
#' @param rank an integer given the rank of the model to be extracted from the collection.
#' @return Send back a object with class \code{\link[=PLNfit.PCA-class]{PLNfit}}.
NULL
PLNfamily$methods(getModel = function(rank){
  id <- match(rank, ranks)
  if (!is.na(id)){
    return(models[[id]]$copy())
  } else {
    stop("No model with such a rank available.")
  }
})

PLNfamily$methods(setCriteria =function() {
  criteria <<- data.frame(rank = ranks, t(sapply(models, function(model) model$criteria)))
})

PLNfamily$methods(setPCA =function() {
  models <<- lapply(models, function(model) {
    m <- PLNfit.PCA$new(model)
    m$setVisualization()
    return(m)
  })
})

#' A plot method for a collection of PLNfit
#'
#' @name PLNfamily_plot
#' @return Produces a plot  representing the evolution of the criteria of the different models considered,
#' highlighting the best model in terms of ICL. If not smooth, you may consider running PLNPCA again with
#' a smaller tolerance for convergence.
NULL
PLNfamily$methods(plot = function() {
  dplot <- melt(criteria[, c("rank", "J", "BIC", "ICL")], id.vars = 1, variable.name = "criterion")
  p <- ggplot(dplot, aes(x=rank, y=value, group=criterion, colour=criterion)) +
        geom_line() + geom_point() + ggtitle("Model selection criteria")
  p <- p + annotate("text", x=criteria$rank, y=min(criteria$J), angle=90, label=paste("R2 =", round(criteria$R2, 2)), size=3, alpha=0.7) +
    geom_vline(xintercept=getBestModel("ICL")$rank, linetype="dashed", alpha=0.5)
  p
})
NULL
PLNfamily$methods(show = function() {
  cat("Poisson-log normal models\n")
  cat("======================================================\n")
  cat(" - Variational approximation: Gaussian with ",type," covariances\n")
  cat(" - Ranks considered: from", min(ranks), "to", max(ranks),"\n")
  cat(" - Best model (regardings ICL): rank =", getBestModel("ICL")$rank, "- R2 =", round(getBestModel("ICL")$criteria["R2"], 2), "\n")
})

PLNfamily$methods(print = function() .self$show())


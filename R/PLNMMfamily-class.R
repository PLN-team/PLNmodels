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
#' @import mclust
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

  if (control$trace > 0) cat("\n Perform GMM on the inceptive model for initializing...")
  M     <- self$inception$var_par$M
  S     <- self$inception$var_par$S
  S_bar <- diag(colMeans(S))

  ## instantiate as many models as clusters
  self$models <- lapply(clusters, function(k) {
    mclust_out <- mclust::Mclust(M, k, modelNames = "EII", verbose = FALSE)

    Tau <- mclust_out$z

    Mu  <- mclust_out$parameters$mean

  ##  Mk <- array(M, c(private$n, private$p, k)) %>% sweep(c(2, 3), Mu, "-")

    Mk <- array(0, c(private$n, private$p, k))

    Sigma <- matrix(0, private$p, private$p)
    for (k_ in 1:k) Sigma <- Sigma + t(Mk[ , ,k_]) %*% diag(Tau[, k_]) %*% Mk[ , ,k_]
    Sigma <- diag(rep(1, private$p)) / private$n

    model <- PLNMMfit$new(
      Theta = self$inception$model_par$Theta,
      Sigma = Sigma,
      Mu    = Mu,
      M     = Mk,
      S     = array(10 * control$lbvar, c(private$n, private$p, k)),
      Tau   = Tau
    )
    return(model)
  })

})

PLNMMfamily$set("public", "optimize",
  function(control) {

  ## set option for call to NLOPT: optim typ, lower bound, tolerance...
  opts <- list("algorithm"   = control$method,
               "maxeval"     = control$maxeval,
               "ftol_rel"    = control$ftol_rel,
               "ftol_abs"    = control$ftol_abs,
               "xtol_rel"    = control$xtol_rel,
               "xtol_abs"    = c(rep(0, private$p*(private$d + 1)), # Theta + Muk
                                 rep(0, private$n*private$p), # Mk
                                 rep(control$xtol_abs, private$n*private$p)), #Sk
               "lower_bound" = c(rep(-Inf, private$p*(private$d + 1)), # Theta + Muk
                                 rep(-Inf, private$n*private$p), # Mk
                                 rep(control$lbvar, private$n*private$p)), # Sk
               "weighted"    = TRUE
  )



  ## ===========================================
  ## GO ALONG THE NUMBER OF CLUSTER (i.e the models)
  for (m in seq_along(self$models))  {
    ## k is the total number of cluster
    k <- self$models[[m]]$k

    ## ===========================================
    ## INITIALISATION
    objective_old <- -Inf
    tau_old <- matrix(NA, private$n, k); par0    <- list()
    for (k_ in 1:k) {
      tau_old[, k_] <- self$models[[m]]$posteriorProb[, k_]
      par0[[k_]] <- c(
        cbind(self$models[[m]]$model_par$Mu[,k_], self$models[[m]]$model_par$Theta),
        self$models[[m]]$var_par$M[, , k_],
        self$models[[m]]$var_par$S[, , k_]
      )
    }

    if (control$trace == 1) {
      cat("\tnumber of cluster =", k, "\r")
      flush.console()
    }
    if (control$trace > 1) {
      cat("\tnumber of cluster =", k, "- iteration:")
    }

    ## ===========================================
    ## OPTIMISATION
    cond <- FALSE; iter <- 0
    objective   <- numeric(control$maxit_out)
    objective_component <- vector("numeric", k)
    convergence <- numeric(control$maxit_out)
    while (!cond) {
      iter <- iter + 1
      if (control$trace > 1) cat("", iter)
      ## UPDATE COMPONENTS OF THE MIXTURE
      for (k_ in 1:k) {
        ## UPDATE COMPONENTS OF THE MIXTURE
        ## weighted-PLN model - probabilities of cluster belonging
        optim_out <- optimization_PLN(
          par0[[k_]],
          self$responses,
          cbind(rep(1,private$n), self$covariates),
          self$offsets,
          tau_old[, k_], opts
        )
        par0[[k_]] <- optim_out$solution
        objective_component[k_] <- optim_out$objective
      }
      ## UPDATE THE POSTERIOR PROBABILITIES
      kappa <- get_kappa(self$responses, cbind(rep(1,private$n),self$covariates), self$offsets, tau_old, par0)
      tau   <- matrix(t(apply(kappa, 1, .softmax)), ncol = k)
      objective[iter]   <- sum(objective_component * colMeans(tau)) - sum(.xlogx(tau))
      # convergence[iter] <- sum(abs(tau_old - tau)) / private$n
      convergence[iter] <- abs(objective_old - objective[iter])/abs(objective[iter])

      ## Check convergence
      if ((convergence[iter] < control$ftol_out) | (iter >= control$maxit_out)) cond <- TRUE
      objective_old <- objective[iter]
      tau_old <- tau
    }

    ## ===========================================
    ## OUTPUT
    ## formating parameters for output

    Theta <- list()
    Mu    <- list()
    M     <- list()
    S     <- list()
    Sigma <- list()

    for (k_ in seq.int(k)) {
      Theta[[k_]] <- matrix(par0[[k_]][1:(private$p*(private$d + 1))], private$p, private$d + 1)
      Mu[[k_]]    <- Theta[[k_]][, 1, drop = FALSE]
      Theta[[k_]] <- Theta[[k_]][, -1, drop = FALSE]
      M[[k_]]     <- matrix(par0[[k_]][private$p*(private$d + 1)   + 1:(private$n*private$p)], private$n,private$p)
      S[[k_]]     <- matrix(par0[[k_]][private$p*(private$d + 1 + private$n) + 1:(private$n*private$p)], private$n, private$p)
      Sigma[[k_]] <- t(M[[k_]]) %*% diag(tau[, k_]) %*% M[[k_]]/ sum(tau[, k_] + diag(tau[, k_] %*% S[[k_]]))
    }

    Theta <- Reduce("+", Theta)
    rownames(Theta) <- colnames(self$responses); colnames(Theta) <- colnames(self$covariates)
    Sigma <- Reduce("+", Sigma)
    rownames(Sigma) <- colnames(Sigma) <- colnames(self$responses)
    M <- array(unlist(M), c(private$n, private$p, k) )
    S <- array(unlist(S), c(private$n, private$p, k) )

    self$models[[m]]$update(
      Sigma = Sigma, Theta = Theta, M = M, S = S, J = -objective[iter], Tau = tau,
          monitoring = list(objective        = objective[1:iter],
                            convergence      = convergence[1:iter],
                            outer_iterations = iter))
    if (control$trace > 1) {
      cat("\r                                                                                    \r")
      flush.console()
    }

  }
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
  cat(" - Best model (regarding BIC): cluster =", self$getBestModel("BIC")$cluster, "- R2 =", round(self$getBestModel("BIC")$R_squared, 2), "\n")
  cat(" - Best model (regarding ICL): cluster =", self$getBestModel("ICL")$cluster, "- R2 =", round(self$getBestModel("ICL")$R_squared, 2), "\n")
})
PLNMMfamily$set("public", "print", function() self$show())

get_kappa <- function(Y, X, O, Tau, par) {

  n <- nrow(Y); p <- ncol(Y); d <- ncol(X)

  Theta <- list()
  M     <- list()
  S     <- list()
  Sigma <- list()
  A     <- list()
  Z     <- list()
  kappa <- matrix(NA, n, ncol(Tau))

  for (k_ in seq.int(ncol(Tau))) {
    Theta[[k_]] <- matrix(par[[k_]][1:(p*d)   ], p, d)
    M[[k_]]     <- matrix(par[[k_]][p*d   + 1:(n*p)], n,p)
    S[[k_]]     <- matrix(par[[k_]][p*(d + n) + 1:(n*p)], n,p)
    Sigma[[k_]] <- t(M[[k_]]) %*% diag(Tau[, k_]) %*% M[[k_]]/ sum(Tau[, k_] + diag(Tau[, k_] %*% S[[k_]]))
  }

  Sigma_bar <- Reduce("+", Sigma)
  Omega_bar <- solve(Sigma_bar)
  pi <- colMeans(Tau)
  log_det_Sigma <- determinant(Sigma_bar, logarithm = TRUE)$modulus

  for (k_ in seq.int(ncol(Tau))) {
    Z[[k_]] <- O +  X %*% t(Theta[[k_]]) + M[[k_]]
    A[[k_]] <- exp(Z[[k_]] + .5 * S[[k_]])
    kappa[, k_] <- rowSums(Y * Z[[k_]]) - rowSums(A[[k_]]) - .5 * rowSums(log(S[[k_]])) -
      .5 * (log_det_Sigma + diag(M[[k_]] %*% Omega_bar %*% t(M[[k_]])) + S[[k_]] %*% diag(Omega_bar)) + log(pi[k_])
  }

  kappa
}

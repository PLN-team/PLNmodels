# ZI optimization:
#
# init_parameters: named list(B, B0, Omega, M, S2, R)
# Y, X, O: observed data, constant
# configuration: named list(algorithm, xtol_rel, xtol_abs, ftol_rel, ftol_abs, maxeval, maxtime)
#
# Returns named list(parameters: named list, nb_iter: int, stop_reason: str)
#
# Configuration :
# Forwarded to nloptr iterative steps (c,e,f).
# Except for algorithm, the keys are optional termination criteria. Not all of them are needed.
# A missing key or a key set to a negative number means disabled.
# xtol_abs: double for every parameter, or named list with values for each parameter.
#
# xtol_abs, xtol_rel, maxeval are also used for overall loop.
# ftol_abs and ftol_rel could maybe be used if defined on J(params), but I'm not sure.
#
# Dimensions are checked only on C++ side.
optimize_zi <- function(init_parameters, Y, X, O, configuration) {

    n <- nrow(Y); p <- ncol(Y); d <- ncol(X)
    
    # Link to the approximate function to optimize Omega ,depending on the target structure
    optim_zipln_Omega <- switch(
      configuration$covariance,
      "spherical" = optim_zipln_Omega_spherical,
      "diagonal"  = optim_zipln_Omega_diagonal,
      "full"      = optim_zipln_Omega_full,
      "sparse"    = function(M, X, B, S) optim_zipln_Omega_sparse(M, X, B, S, rho = configuration$rho)
    )

    optim_zipln_zipar <- switch(
      configuration$ziparam,
      "single" = function(init_B0, X, R, config) list(Pi = matrix(    mean(R), n, p)              , B0 = matrix(NA, d, p)),
      "row"    = function(init_B0, X, R, config) list(Pi = matrix(rowMeans(R), n, p)              , B0 = matrix(NA, d, p)),
      "col"    = function(init_B0, X, R, config) list(Pi = matrix(colMeans(R), n, p, byrow = TRUE), B0 = matrix(NA, d, p)),
      "covar"  = optim_zipln_zipar_covar
    )
    
    maxit_out <- if("maxit_out" %in% names(configuration)) { configuration$maxit_out } else { 50 }

    # Main loop
    nb_iter <- 0
    parameters <- init_parameters
    criterion <- vector("numeric", maxit_out)
    vloglik <- -Inf; objective <- Inf
    repeat {
        # Check maxeval
        if(maxit_out >= 0 && nb_iter >= maxit_out) {
            return(list(
                parameters = parameters,
                nb_iter = nb_iter,
                stop_reason = "maximum number of iterations reached",
                criterion = criterion[1:nb_iter],
                vloglik = vloglik
            ))
        }

        # M Step
        new_Omega <- optim_zipln_Omega(
            M = parameters$M, X = X, B = parameters$B, S = parameters$S
        )
        new_B <- optim_zipln_B(
          M = parameters$M, X = X, Omega = new_Omega, configuration
        )
        
        optim_new_zipar <- optim_zipln_zipar(
          init_B0 = parameters$B0, X = X, R = parameters$R, config = configuration
        )
        new_B0 <- optim_new_zipar$B0
        new_Pi <- optim_new_zipar$Pi

        # VE Step
        new_R <- optim_zipln_R(
          Y = Y, X = X, O = O, M = parameters$M, S = parameters$S, Pi = new_Pi
        )
        optim_new_M <- optim_zipln_M(
          init_M = parameters$M,
          Y = Y, X = X, O = O, R = new_R, S = parameters$S, B = new_B, Omega = new_Omega,
          configuration = configuration
        )
        new_M <- optim_new_M$M
        optim_new_S <- optim_zipln_S(
          init_S = parameters$S,
          O = O, M = new_M, R = new_R, B = new_B, diag_Omega = diag(new_Omega),
          configuration = configuration
        )
        new_S <- optim_new_S$S
        # print(optim_new_S$status)

        # Check convergence
        new_parameters <- list(
            Omega = new_Omega, B = new_B, B0 = new_B0, Pi = new_Pi,
            R = new_R, M = new_M, S = new_S
        )
        nb_iter <- nb_iter + 1

        vloglik <- zipln_vloglik(
          Y, X, O, new_Pi, new_Omega, new_B, new_R, new_M, new_S
        )

        criterion[nb_iter] <- new_objective <- -sum(vloglik)

        objective_converged <-
            (objective - new_objective) < configuration$ftol_out |
            (objective - new_objective)/abs(new_objective) < configuration$ftol_out

        parameters_converged <- parameter_list_converged(
            parameters, new_parameters,
            xtol_abs = configuration$xtol_abs, xtol_rel = configuration$xtol_rel
        )

        if (parameters_converged | objective_converged) {
            return(list(
                parameters = new_parameters,
                nb_iter = nb_iter,
                stop_reason = "converged",
                criterion = criterion[1:nb_iter],
                vloglik = vloglik
            ))
        }

        parameters <- new_parameters
        objective  <- new_objective
    }
}

#' @importFrom glassoFast glassoFast
optim_zipln_Omega_sparse <- function(M, X, B, S, rho) {
  n <- nrow(M); p <- ncol(M)
  glassoFast::glassoFast( crossprod(M - X %*% B)/n + diag(colMeans(S * S), p, p), rho = rho )$wi
}

#' @importFrom glmnet glmnet
optim_zipln_B <- function(M, X, Omega, config) {

  if(config$lambda > 0) {
    if (!is.null(config$ind_intercept)) {
      m_bar <- colMeans(M)
      x_bar <- colMeans(X[, -config$ind_intercept])
      X <- scale(X[, -config$ind_intercept], x_bar, FALSE)
      M <- scale(M, m_bar, FALSE)
    }
    p <- ncol(M); d <- ncol(X)
    if (d > 0) {
      Omega12 <- chol(Omega)
      y <- as.vector(M %*% t(Omega12))
      x <- kronecker(Omega12, X)
      glmnet_out <- glmnet(x, y, lambda = config$lambda, intercept = FALSE, standardize = FALSE)
      B <- matrix(as.numeric(glmnet_out$beta), nrow = d, ncol = p)
    } else {
      B <- matrix(0, nrow = d, ncol = p)
    }

    if (!is.null(config$ind_intercept)) {
      mu0 <- m_bar - as.vector(crossprod(B, x_bar))
      B <- rbind(mu0, B)
    }

  } else {
    B <- optim_zipln_B_dense(M, X)
  }
  B
}

# Test convergence for a named list of parameters
# oldp, newp: named list of parameters
# xtol_rel: double ; negative or NULL = disabled
# xtol_abs: double ; negative or NULL = disabled
# Returns boolean
parameter_list_converged <- function(oldp, newp, xtol_abs = NULL, xtol_rel = NULL) {
    # Strategy is to compare each pair of list elements with matching names.
    # Named lists are just vectors (T,str) using order of insertion.
    # mapply() is handy to do the pair tests, but it works on the underlying vector order (ignoring names).
    # So reorder lists by their names to use mapply.
    stopifnot(is.list(oldp), is.list(newp))
    oldp <- oldp[order(names(oldp))]
    newp <- newp[order(names(newp))]
    stopifnot(all(names(oldp) == names(newp)))

    # Check convergence with xtol_rel if enabled
    if(is.double(xtol_rel) && xtol_rel > 0) {
        if(all(mapply(function(o, n) { all(abs(n - o) <= xtol_rel * abs(o)) }, oldp, newp))) {
            return(TRUE)
        }
    }

    # Check convergence with xtol_abs (homogeneous) if enabled
    if(is.double(xtol_abs) && xtol_abs > 0) {
        if(all(mapply(function(o, n) { all(abs(n - o) <= xtol_abs) }, oldp, newp))) {
            return(TRUE)
        }
    }

    # Check convergence with xtol_abs as list(xtol_abs for each param_name)
    if(is.list(xtol_abs)) {
        xtol_abs <- xtol_abs[order(names(xtol_abs))]
        stopifnot(all(names(oldp) == names(xtol_abs)))
        # Due to the possible presence of NULLs, mapply may return a list. unlist allows all() to operate anyway.
        if(all(unlist(mapply(
            function(o, n, tol) {
                if((is.double(tol) && tol > 0) || is.matrix(tol)) {
                    all(abs(n - o) <= tol)
                } else {
                    NULL # Ignore comparison in outer all()
                }
            },
            oldp, newp, xtol_abs
        )))) {
            return(TRUE)
        }
    }

    # If no criteria has triggered, indicate no convergence
    FALSE
}


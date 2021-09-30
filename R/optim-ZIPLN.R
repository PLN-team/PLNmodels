# ZI optimization:
#
# init_parameters: named list(Theta, Theta0, Omega, M, S2, Pi)
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

    # Prepare configuration structs for steps with nloptr
    config_for <- function(parameter_name) {
        # This relies on "copy on write" semantics in R ; configuration is copied when modified
        if("xtol_abs" %in% names(configuration)) {
            xtol_abs <- configuration$xtol_abs
            if(is.double(xtol_abs)) {
                # double means applying the value for everything.
                # nloptr wrapper already treats it as such, keep it unchanged.
            } else if(is.list(xtol_abs)) {
                # list(param_names = values) means using specific values (double / matrix) for each parameter.
                # extract value for specific parameter ; nlopt wrapper will check dimensions
                configuration$xtol_abs <- xtol_abs[[parameter_name]]
            } else {
                warning("optimize_zi: configuration$xtol_abs is neither a list or a double")
            }
        }
        configuration
    }
    configuration_step_c <- config_for("Theta0")
    configuration_step_e <- config_for("M")
    configuration_step_f <- config_for("S")

    # Follow the convention of maxeval = -1 => disabled
    maxeval <- if("maxeval" %in% names(configuration)) { configuration$maxeval } else { -1 }

    # Main loop
    nb_iter <- 0
    parameters <- init_parameters
    criterion <- vector("numeric", 100)
    objective <- Inf
    repeat {
        # Check maxeval
        if(maxeval >= 0 && nb_iter >= maxeval) {
            return(list(
                parameters = parameters,
                nb_iter = nb_iter,
                stop_reason = "maxeval"
            ))
        }
        # Steps
        new_Omega <- cpp_optimize_zi_Omega(
            M = parameters$M, X = X, Theta = parameters$Theta, S = parameters$S
        )
        new_Theta <- cpp_optimize_zi_Theta(
            M = parameters$M, X = X
        )
        new_Theta0 <- cpp_optimize_zi_Theta0(
            init_Theta0 = parameters$Theta0,
            X = X, Pi = parameters$Pi,
            configuration = configuration_step_c
        )$Theta0
        new_Pi <- cpp_optimize_zi_Pi(
            Y = Y, X = X, O = O, M = parameters$M, S = parameters$S, Theta0 = new_Theta0
        )
        new_M <- cpp_optimize_zi_M(
            init_M = parameters$M,
            Y = Y, X = X, O = O, Pi = new_Pi, S = parameters$S, Theta = new_Theta, Omega = new_Omega,
            configuration = configuration_step_e
        )$M
        new_S <- cpp_optimize_zi_S(
            init_S = parameters$S,
            O = O, M = new_M, Pi = new_Pi, Theta = new_Theta, Omega = new_Omega,
            configuration = configuration_step_f
        )$S


        # Check convergence
        new_parameters <- list(
            Omega = new_Omega, Theta = new_Theta, Theta0 = new_Theta0,
            Pi = new_Pi, M = new_M, S = new_S
        )
        nb_iter <- nb_iter + 1

        criterion[nb_iter] <- new_objective <- get_objective(Y, X, O, new_parameters)

        objective_converged <-
            (objective - new_objective) < configuration$ftol_out |
            (objective - new_objective)/abs(new_objective) < configuration$ftol_out

        parameters_converged <- parameter_list_converged(
            parameters, new_parameters,
            xtol_abs = configuration$xtol_abs, xtol_rel = configuration$xtol_rel
        )

        if (parameters_converged | objective_converged) {
            return(list(
                parameters = parameters,
                nb_iter = nb_iter,
                stop_reason = "converged",
                criterion = criterion[1:nb_iter],
                vloglik = get_vloglik(Y, X, O, new_parameters)
            ))
        }

        if(nb_iter >= configuration$maxit_out) {
            return(list(
                parameters = parameters,
                nb_iter = nb_iter,
                stop_reason = "maximum number of iterations reached",
                criterion = criterion[1:nb_iter],
                vloglik = get_vloglik(Y, X, O, new_parameters)
            ))
        }

        parameters <- new_parameters
        objective  <- new_objective
    }
}

get_objective <- function(Y, X, O, par) {

    A <- exp(O + par$M + .5 * par$S^2)
    n <- nrow(Y)
    logdet_Omega <-determinant(par$Omega, logarithm = TRUE)$modulus
    res <- - sum ((1 - par$Pi) * ( Y * (O + par$M) - A - .logfactorial(Y)) ) -
      sum(par$Pi* X %*% par$Theta0) + sum( .xlogx(par$Pi) + .xlogx(1-par$Pi) ) -
        .5 * n * logdet_Omega - .5 * sum(log(par$S^2)) + sum(  log (1 + exp(X %*% par$Theta0) ))
    res
}

get_vloglik <- function(Y, X, O, par) {

    Z <- O + par$M
    A <- exp(Z + .5 * par$S^2)
    p <- ncol(Y)
    logdet_Omega <- determinant(par$Omega, logarithm = TRUE)$modulus
    res <- rowSums( (1 - par$Pi) * ( ( Y * Z - A  ) - .logfactorial(Y) ) ) +
      rowSums(par$Pi* X %*% par$Theta0) - rowSums( .xlogx(par$Pi) + .xlogx(1-par$Pi) ) +
        .5 * logdet_Omega + .5 * rowSums(log(par$S^2)) - rowSums( log (1 + exp(X %*% par$Theta0) ) ) + p
    res
}

# Test convergence for a named list of parameters
# oldp, newp: named list of parameters
# xtol_rel: double ; negative or NULL = disabled
# xtol_abs: double (negative or NULL = disabled) or list of double/matrices for each parameter (any negative is disabled)
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

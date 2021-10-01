context("test-plnzi")

epsilon = 1e-6


test_that("PLNzi: steps are dimensionally correct", {
    # Choose different characteristic dimensions to detect mismatches
    n = 20
    p = 6
    d = 4

    # Prepare test values
    M = S = Pi = Y = O = matrix(1., nrow=n, ncol=p)
    X = matrix(runif(n*d), nrow=n, ncol=d)
    Theta = Theta0 = matrix(1., nrow=d, ncol=p)
    Omega = diag(p)

    # Steps
    expect_identical(
        dim(Omega),
        dim(PLNmodels:::cpp_optimize_zi_Omega_full(M = M, X = X, Theta = Theta, S = S))
    )
    expect_identical(
        dim(Omega),
        dim(PLNmodels:::cpp_optimize_zi_Omega_spherical(M = M, X = X, Theta = Theta, S = S))
    )
    expect_identical(
        dim(Omega),
        dim(PLNmodels:::cpp_optimize_zi_Omega_diagonal(M = M, X = X, Theta = Theta, S = S))
    )

    expect_identical(
        dim(Theta),
        dim(PLNmodels:::cpp_optimize_zi_Theta(M = M, X = X))
    )

    step_c = PLNmodels:::cpp_optimize_zi_Theta0(
        init_Theta0 = Theta0, X = X, Pi = Pi,
        configuration = list(
            algorithm = "CCSAQ",
            xtol_rel = epsilon, xtol_abs = epsilon,
            ftol_rel = epsilon, ftol_abs = epsilon,
            maxeval = 3 # make it stop fast
        )
    )
    expect_identical(dim(Theta0), dim(step_c$Theta0))

    expect_identical(
        dim(Pi),
        dim(PLNmodels:::cpp_optimize_zi_Pi(
            Y = Y, X = X, O = O, M = M, S = S, Theta0 = Theta0
        ))
    )

    step_e = PLNmodels:::cpp_optimize_zi_M(
        init_M = M, Y = Y, X = X, O = O, Pi = Pi, S = S, Theta = Theta, Omega = Omega,
        configuration = list(
            algorithm = "CCSAQ",
            xtol_rel = epsilon, xtol_abs = M * epsilon, # demonstrate per-cell xtol_abs
            ftol_rel = epsilon, ftol_abs = epsilon,
            maxeval = 3 # make it stop fast
        )
    )
    expect_identical(dim(M), dim(step_e$M))

    step_f = PLNmodels:::cpp_optimize_zi_S(
        init_S = S, O = O, M = M, Pi = Pi, Theta = Theta, diag_Omega = diag(Omega),
        configuration = list(
            algorithm = "CCSAQ",
            xtol_rel = epsilon, xtol_abs = epsilon,
            ftol_rel = epsilon, ftol_abs = epsilon,
            maxeval = 3 # make it stop fast
        )
    )
    expect_identical(dim(S), dim(step_f$S))
})

test_that("PLNzi: parameter_list_convergence", {
    p = list(d = 4., mat = diag(10) + 1)
    p_abs = list(d = p$d + 2 * epsilon, mat = p$mat - 2 * epsilon) # to check with epsilon
    p_rel = list(d = p$d + 0.1, mat = p$mat - 0.1) # to check with xtol_rel = 0.1
    p_switched = list(mat = p$mat, d = p$d) # Check order does not matter

    # With no criterion set, no convergence
    expect_false(PLNmodels:::parameter_list_converged(p, p))
    expect_false(PLNmodels:::parameter_list_converged(p, p_switched))
    expect_false(PLNmodels:::parameter_list_converged(p, p, xtol_rel = -1))
    expect_false(PLNmodels:::parameter_list_converged(p, p, xtol_abs = -1))
    expect_false(PLNmodels:::parameter_list_converged(p, p, xtol_abs = -1, xtol_rel = -1))

    # Any criterion should lead to convergence for (p,p)
    expect_true(PLNmodels:::parameter_list_converged(p, p, xtol_rel = epsilon))
    expect_true(PLNmodels:::parameter_list_converged(p, p_switched, xtol_rel = epsilon))
    expect_true(PLNmodels:::parameter_list_converged(p, p, xtol_abs = epsilon))
    expect_true(PLNmodels:::parameter_list_converged(p, p_switched, xtol_abs = epsilon))
    expect_true(PLNmodels:::parameter_list_converged(p, p, xtol_abs = list(d = epsilon, mat = epsilon)))
    expect_true(PLNmodels:::parameter_list_converged(p, p, xtol_rel = epsilon, xtol_abs = epsilon))

    # (p,p_abs) cases, test xtol_abs
    expect_false(PLNmodels:::parameter_list_converged(p, p_abs, xtol_abs = epsilon))
    expect_true(PLNmodels:::parameter_list_converged(p, p_abs, xtol_abs = 3 * epsilon))
    expect_false(PLNmodels:::parameter_list_converged(p, p_abs, xtol_abs = list(d = 3 * epsilon, mat = epsilon)))
    expect_false(PLNmodels:::parameter_list_converged(p, p_abs, xtol_abs = list(d = epsilon, mat = 3 * epsilon)))
    expect_true(PLNmodels:::parameter_list_converged(p, p_abs, xtol_abs = list(d = 3 * epsilon, mat = 3 * epsilon)))
    expect_true(PLNmodels:::parameter_list_converged(p, p_abs, xtol_abs = list(d = 3 * epsilon, mat = -1)))
    expect_true(PLNmodels:::parameter_list_converged(p, p_abs, xtol_abs = list(d = 3 * epsilon, mat = NULL)))

    # (p,prel) cases, test xtol_rel
    expect_false(PLNmodels:::parameter_list_converged(p, p_rel, xtol_rel = epsilon))
    expect_true(PLNmodels:::parameter_list_converged(p, p_rel, xtol_rel = 0.1))

    # Check error handling in case of mismatched lists
    expect_error(PLNmodels:::parameter_list_converged(p, list(d = 3), xtol_abs = epsilon), "names")
    expect_error(PLNmodels:::parameter_list_converged(p, p_abs, xtol_abs = list(d = 3 * epsilon)), "names")
})

test_that("PLNzi: optimize_zi is dimensionally correct", {
    # Choose different caracteristic dimensions to detect mismatches
    n = 20
    p = 6
    d = 3

    # Prepare test values
    M = S = Pi = Y = O = matrix(1., nrow=n, ncol=p)
    X = matrix(runif(n*d), nrow=n, ncol=d)
    Theta = Theta0 = matrix(1., nrow=d, ncol=p)
    Omega = diag(p)

    parameters = list(
        Theta = Theta, Theta0 = Theta0, Omega = Omega,
        M = M, S = S, Pi = Pi
    )

    # Test dimensions (and also example of how to use optimize_zi)
    result = PLNmodels:::optimize_zi(
        init_parameters = parameters, Y = Y, X = X, O = O,
        configuration = list(
            covariance = "full",
            maxit_out  = 3,
            ftol_out = 1e-4,
            algorithm = "CCSAQ",
            xtol_rel = epsilon, xtol_abs = epsilon,
            ftol_rel = epsilon, ftol_abs = epsilon
        )
    )
    expect_identical(result$stop_reason, "maximum number of iterations reached")
    expect_true(setequal(names(parameters), names(result$parameters)))
    expect_identical(result$nb_iter, 3)
})

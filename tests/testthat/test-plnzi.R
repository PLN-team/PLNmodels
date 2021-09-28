context("test-plnzi")

epsilon = 1e-6

test_that("PLNzi: steps are dimensionally correct", {
    # Choose different characteristic dimensions to detect mismatches
    n = 3
    p = 9
    d = 6

    # Prepare test values
    M = S2 = Pi = Y = O = matrix(1., nrow=n, ncol=p)
    X = matrix(runif(n*d), nrow=n, ncol=d)
    Theta = Theta0 = matrix(1., nrow=d, ncol=p)
    Omega = diag(p)

    # Steps
    expect_identical(
        dim(Omega),
        dim(cpp_optimize_zi_step_a(M = M, X = X, Theta = Theta, S2 = S2))
    )

    expect_identical(
        dim(Theta),
        dim(cpp_optimize_zi_step_b(M = M, X = X))
    )

    step_c = cpp_optimize_zi_step_c(
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
        dim(cpp_optimize_zi_step_d(
            Y = Y, X = X, O = O, M = M, S2 = S2, Theta0 = Theta0
        ))
    )

    step_e = cpp_optimize_zi_step_e(
        init_M = M, Y = Y, X = X, O = O, Pi = Pi, S2 = S2, Theta = Theta, Omega = Omega,
        configuration = list(
            algorithm = "CCSAQ",
            xtol_rel = epsilon, xtol_abs = M * epsilon, # demonstrate per-cell xtol_abs
            ftol_rel = epsilon, ftol_abs = epsilon,
            maxeval = 3 # make it stop fast
        )
    )
    expect_identical(dim(M), dim(step_e$M))

    step_f = cpp_optimize_zi_step_f(
        init_S2 = S2, O = O, M = M, Pi = Pi, Theta = Theta, Omega = Omega,
        configuration = list(
            algorithm = "CCSAQ",
            xtol_rel = epsilon, xtol_abs = epsilon,
            ftol_rel = epsilon, ftol_abs = epsilon,
            maxeval = 3 # make it stop fast
        )
    )
    expect_identical(dim(S2), dim(step_f$S2))
})

test_that("PLNzi: parameter_list_convergence", {
    p = list(d = 42., mat = diag(3) + 3)
    p_abs = list(d = p$d + 2 * epsilon, mat = p$mat - 2 * epsilon) # to check with epsilon
    p_rel = list(d = p$d + 0.1, mat = p$mat - 0.1) # to check with xtol_rel = 0.1
    p_switched = list(mat = p$mat, d = p$d) # Check order does not matter

    # With no criterion set, no convergence
    expect_false(parameter_list_converged(p, p))
    expect_false(parameter_list_converged(p, p_switched))
    expect_false(parameter_list_converged(p, p, xtol_rel = -1))
    expect_false(parameter_list_converged(p, p, xtol_abs = -1))
    expect_false(parameter_list_converged(p, p, xtol_abs = -1, xtol_rel = -1))

    # Any criterion should lead to convergence for (p,p)
    expect_true(parameter_list_converged(p, p, xtol_rel = epsilon))
    expect_true(parameter_list_converged(p, p_switched, xtol_rel = epsilon))
    expect_true(parameter_list_converged(p, p, xtol_abs = epsilon))
    expect_true(parameter_list_converged(p, p_switched, xtol_abs = epsilon))
    expect_true(parameter_list_converged(p, p, xtol_abs = list(d = epsilon, mat = epsilon)))
    expect_true(parameter_list_converged(p, p, xtol_rel = epsilon, xtol_abs = epsilon))

    # (p,p_abs) cases, test xtol_abs
    expect_false(parameter_list_converged(p, p_abs, xtol_abs = epsilon))
    expect_true(parameter_list_converged(p, p_abs, xtol_abs = 3 * epsilon))
    expect_false(parameter_list_converged(p, p_abs, xtol_abs = list(d = 3 * epsilon, mat = epsilon)))
    expect_false(parameter_list_converged(p, p_abs, xtol_abs = list(d = epsilon, mat = 3 * epsilon)))
    expect_true(parameter_list_converged(p, p_abs, xtol_abs = list(d = 3 * epsilon, mat = 3 * epsilon)))
    expect_true(parameter_list_converged(p, p_abs, xtol_abs = list(d = 3 * epsilon, mat = -1)))
    expect_true(parameter_list_converged(p, p_abs, xtol_abs = list(d = 3 * epsilon, mat = NULL)))

    # (p,prel) cases, test xtol_rel
    expect_false(parameter_list_converged(p, p_rel, xtol_rel = epsilon))
    expect_true(parameter_list_converged(p, p_rel, xtol_rel = 0.1))

    # Check error handling in case of mismatched lists
    expect_error(parameter_list_converged(p, list(d = 3), xtol_abs = epsilon), "names")
    expect_error(parameter_list_converged(p, p_abs, xtol_abs = list(d = 3 * epsilon)), "names")
})

test_that("PLNzi: optimize_zi is dimensionally correct", {
    # Choose different caracteristic dimensions to detect mismatches
    n = 3
    p = 9
    d = 6

    # Prepare test values
    M = S2 = Pi = Y = O = matrix(1., nrow=n, ncol=p)
    X = matrix(runif(n*d), nrow=n, ncol=d)
    Theta = Theta0 = matrix(1., nrow=d, ncol=p)
    Omega = diag(p)

    parameters = list(
        Theta = Theta, Theta0 = Theta0, Omega = Omega,
        M = M, S2 = S2, Pi = Pi
    )

    # Test dimensions (and also example of how to use optimize_zi)
    result = optimize_zi(
        init_parameters = parameters, Y = Y, X = X, O = O,
        configuration = list(
            algorithm = "CCSAQ",
            xtol_rel = epsilon, xtol_abs = epsilon,
            ftol_rel = epsilon, ftol_abs = epsilon,
            maxeval = 3 # make it stop fast
        )
    )
    expect_identical(result$stop_reason, "maxeval")
    expect_true(setequal(names(parameters), names(result$parameters)))
    expect_identical(result$nb_iter, 3)
})

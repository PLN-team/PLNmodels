context("test-plnzi")

epsilon = 1e-6

test_that("PLNzi: steps are dimensionally correct", {
    # Choose different caracteristic dimensions to detect mismatches
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
            xtol_rel = epsilon, xtol_abs = M * epsilon, # per-cell xtol_abs
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
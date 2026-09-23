###############################################################################
## The in-package graphical Lasso (src/graphical_lasso.h), which replaced
## glassoFast in PLNnetwork and ZIPLNnetwork. It is a port of glassoFast's
## Fortran: the equivalence is checked directly when glassoFast is installed,
## and otherwise pinned by the properties that matter -- the returned matrix
## solves the penalized problem, and the closed-form cases come out exact.
###############################################################################

## a well-conditioned covariance: 3n degrees of freedom plus a small ridge, so
## that the tests measure the solver rather than the conditioning of S
rand_S <- function(n, m = 3 * n) {
  A <- matrix(rnorm(n * m), n)
  tcrossprod(A) / m + diag(0.1, n)
}

## -log det(Theta) + tr(S Theta) + || rho o Theta ||_1
glasso_objective <- function(S, rho, Theta) {
  ld <- determinant(Theta, logarithm = TRUE)
  if (ld$sign <= 0) return(Inf)
  -as.numeric(ld$modulus) + sum(S * Theta) + sum(abs(rho * Theta))
}

off_diag_weights <- function(n) {
  w <- matrix(1, n, n)
  diag(w) <- 0
  w
}

## support of the off-diagonal part, as read along a penalty path
support <- function(Theta) { diag(Theta) <- 0; which(Theta != 0) }

## penalty path from the largest to 1e-4 times the largest off-diagonal entry
penalty_path <- function(S, length = 50) {
  x <- abs(S[upper.tri(S)]); hi <- max(x); lo <- min(max(min(x[x > 0]), hi * 1e-4), hi)
  10^seq(log10(hi), log10(lo), length.out = length)
}

test_that("the solution is symmetric, positive definite, and minimizes the objective", {
  set.seed(101)
  for (n in c(4, 12, 30)) {
    S   <- rand_S(n)
    rho <- 0.08 * off_diag_weights(n)
    Theta <- graphical_lasso(S, rho)$wi

    expect_true(isSymmetric(Theta))
    expect_gt(min(eigen(Theta, symmetric = TRUE, only.values = TRUE)$values), 0)

    ## no small symmetric perturbation staying in the PD cone does better
    base  <- glasso_objective(S, rho, Theta)
    worse <- 0
    for (k in 1:150) {
      P <- matrix(0, n, n)
      i <- sample(n, 1); j <- sample(n, 1)
      P[i, j] <- P[j, i] <- rnorm(1, sd = 1e-3)
      if (glasso_objective(S, rho, Theta + P) < base - 1e-10) worse <- worse + 1
    }
    expect_equal(worse, 0)
  }
})

test_that("the penalty actually sparsifies, monotonically in rho", {
  set.seed(102)
  n <- 20
  S <- rand_S(n)
  nnz <- vapply(c(0.001, 0.05, 0.2, 0.6, 2), function(pen) {
    Theta <- graphical_lasso(S, pen * off_diag_weights(n))$wi
    sum(abs(Theta[upper.tri(Theta)]) > 1e-10)
  }, numeric(1))

  expect_true(all(diff(nnz) <= 0))
  expect_gt(nnz[1], nnz[length(nnz)])
  expect_equal(nnz[length(nnz)], 0) # a big enough penalty empties the network
})

test_that("an unpenalized problem returns the plain inverse", {
  set.seed(103)
  n <- 10
  S <- rand_S(n)
  ## at the default thr = 1e-4 this is an iterative solve stopped early, not an
  ## exact inverse, hence the tightened threshold
  expect_equal(graphical_lasso(S, matrix(0, n, n), thr = 1e-12)$wi, solve(S),
               tolerance = 1e-6)
})

test_that("a separable problem is solved exactly rather than through the recursion", {
  ## No off-diagonal mass: Theta is diagonal, with entries 1 / (S_ii + rho_ii).
  ## This is the case glassoFast gets wrong -- it drops S_ii and returns
  ## 1 / max(rho_ii, eps), i.e. ~9.09e15 on an unpenalized diagonal.
  S <- diag(c(4, 9, 16))
  expect_equal(diag(graphical_lasso(S, 0.1 * off_diag_weights(3))$wi), 1 / diag(S))
  expect_equal(graphical_lasso(matrix(4, 1, 1), 0)$wi[1], 0.25)
  ## a penalized diagonal shifts it, and is still exact
  expect_equal(diag(graphical_lasso(S, diag(1, 3))$wi), 1 / (diag(S) + 1))
})

test_that("degenerate input is reported rather than hung on or silently accepted", {
  ## glassoFast loops forever on a non-finite input: its inner loop only exits
  ## on `dlx < thrLasso`, which is never true once dlx is NaN
  res <- graphical_lasso(matrix(c(1, NA, NA, 1), 2), 0.1)
  expect_true(all(is.na(res$wi)))
  expect_false(res$converged)
  expect_equal(res$status, "degenerate")

  ## a zero-variance coordinate is what would divide by zero inside the sweep
  res0 <- graphical_lasso(matrix(c(0, 0, 0, 1), 2), matrix(0, 2, 2))
  expect_true(all(is.na(res0$wi)))
  expect_false(res0$converged)
  expect_equal(res0$status, "degenerate")
})

###############################################################################
## Stagnation: on an ill-conditioned covariance the sweeps settle into a small
## limit cycle, so the stopping criterion is never met while the solution no
## longer moves. The solver detects this instead of spending its whole budget.
###############################################################################

test_that("converging problems are left alone by the stagnation detection", {
  ## the detector must never fire where the criterion is reachable: the result
  ## has to be exactly what it is with the detection switched off
  set.seed(311)
  for (n in c(10, 25)) {
    S <- rand_S(n)
    for (pen in c(0.3, 0.05, 0.005)) {
      rho <- pen * off_diag_weights(n)
      g <- graphical_lasso(S, rho)
      h <- graphical_lasso(S, rho, stall_patience = Inf)
      expect_equal(g$status, "converged")
      expect_true(g$converged)
      expect_lte(g$delta, 1)          # delta <= 1 is exactly the criterion
      expect_identical(g$wi, h$wi)
      expect_identical(g$niter, h$niter)
    }
  }
})

test_that("a cycling solve is detected, reported, and lands where grinding on would", {
  skip_on_cran()
  ## 80 most abundant oak species: dense enough, at a low penalty, for the
  ## sweeps to start cycling -- glassoFast burns all 10000 of them here
  data(oaks)
  Y <- as.matrix(oaks$Abundance)
  S <- cov(log1p(Y[, order(-colMeans(Y))[1:80]])) + diag(1e-3, 80)
  rho <- max(abs(S[upper.tri(S)])) * 1e-3

  stalled <- graphical_lasso(S, rho)
  expect_equal(stalled$status, "stalled")
  expect_false(stalled$converged)
  expect_lt(stalled$niter, 5000)  # stopped well short of maxit
  expect_gt(stalled$delta, 1)     # the criterion was not met, by definition
  expect_lt(stalled$delta, 10)    # but it came close: this is a cycle, not a failure
  expect_false(anyNA(stalled$wi))

  ## Grinding on does not buy accuracy: the iterates wander inside the cycle,
  ## so a much longer run is another point of it, not a better answer. What is
  ## pinned here is the amplitude -- same number of edges, values within the
  ## cycle -- not entry-by-entry agreement, which the cycle does not provide
  ## (a couple of borderline edges flip whichever sweep one stops at).
  ground <- graphical_lasso(S, rho, maxit = 10000, stall_patience = Inf)
  expect_equal(length(support(stalled$wi)), length(support(ground$wi)))
  expect_equal(stalled$wi, ground$wi, tolerance = 1e-2)
  expect_lt(max(abs(stalled$wi - ground$wi)) / max(abs(ground$wi)), 1e-2)
})

test_that("the per-sweep criterion can be traced", {
  set.seed(312)
  S <- rand_S(12)
  g <- graphical_lasso(S, 0.05 * off_diag_weights(12), trace = TRUE)
  expect_length(g$dw_trace, g$niter)
  expect_true(all(is.finite(g$dw_trace)))
  expect_length(graphical_lasso(S, 0.05 * off_diag_weights(12))$dw_trace, 0)
})

test_that("a nearly collapsed covariance terminates without hanging or corrupting output", {
  ## Rank 2, max|diag| ~ 2.4e-8: glassoFast processes the first ~45 penalties of
  ## this path (on Linux/glibc), then never returns from around the 46th. The
  ## exact cutoff is platform-dependent -- it sits at the edge of machine
  ## precision, where BLAS/LAPACK rounding differences (observed between glibc
  ## and macOS's Accelerate) shift it by a few penalties either way -- so only
  ## the well-conditioned (large penalty) end and the overall termination /
  ## honest reporting are pinned here; see recipe A in the handoff document
  ## this test is transcribed from for the exact glibc counts.
  set.seed(1); p <- 20
  for (sc in c(3e-4, 2e-4, 1.5e-4, 1e-4, 7e-5)) L <- matrix(rnorm(p * 2), p, 2) * sc
  L <- matrix(rnorm(p * 2), p, 2) * 5e-5
  S <- tcrossprod(L)
  x <- abs(S[upper.tri(S)]); hi <- max(x); lo <- max(min(x[x > 0]), hi * 1e-4)
  fits <- lapply(10^seq(log10(hi), log10(lo), length.out = 50),
                 function(rho) graphical_lasso(S, rho))
  converged <- vapply(fits, `[[`, logical(1), "converged")
  expect_true(all(converged[1:20]))
  expect_gt(mean(converged), 0.5)
  expect_false(any(vapply(fits, function(f) anyNA(f$wi), logical(1))))
})

test_that("input is validated, and a scalar rho and a constant matrix rho agree", {
  set.seed(104)
  S <- rand_S(8)
  expect_equal(graphical_lasso(S, 0.1)$wi, graphical_lasso(S, matrix(0.1, 8, 8))$wi)
  expect_error(graphical_lasso(S, matrix(0.1, 3, 3)), "same dimensions")
  expect_error(graphical_lasso(S, -1), "non-negative")
  expect_error(graphical_lasso(S[, 1:3], 0.1), "square")
})

test_that("per-pair weights are respected: a heavily penalized pair is zeroed", {
  set.seed(105)
  n <- 8
  S <- rand_S(n)
  rho <- 0.02 * off_diag_weights(n)
  rho[2, 5] <- rho[5, 2] <- 50
  Theta <- graphical_lasso(S, rho)$wi
  expect_equal(Theta[2, 5], 0)
  expect_gt(sum(abs(Theta[upper.tri(Theta)]) > 1e-10), 0) # the rest survives
})

test_that("a warm start lands on the same solution as a cold one", {
  set.seed(106)
  n <- 15
  S1 <- rand_S(n)
  S2 <- S1 + 0.02 * rand_S(n) # a nearby problem
  rho <- 0.05 * off_diag_weights(n)

  first <- graphical_lasso(S1, rho)
  cold  <- graphical_lasso(S2, rho)
  warm  <- graphical_lasso(S2, rho, w_init = first$w, wi_init = first$wi)

  ## the stopping rule measures per-sweep progress, so starting closer exits
  ## sooner and slightly short of where a cold start would stop
  expect_equal(cold$wi, warm$wi, tolerance = 1e-2)
  expect_lte(warm$niter, cold$niter)

  ## tightening the threshold collapses the difference
  tight_cold <- graphical_lasso(S2, rho, thr = 1e-10)
  tight_warm <- graphical_lasso(S2, rho, thr = 1e-10, w_init = first$w, wi_init = first$wi)
  expect_equal(tight_cold$wi, tight_warm$wi, tolerance = 1e-6)
})

test_that("w and wi are inverses of each other at convergence", {
  set.seed(107)
  n <- 10
  S <- rand_S(n)
  res <- graphical_lasso(S, 0.05 * off_diag_weights(n))
  expect_true(res$converged)
  expect_equal(res$w %*% res$wi, diag(n), tolerance = 1e-3)
})

test_that("an R time limit stops a long solve with R's own error", {
  skip_on_cran()
  set.seed(5); p <- 300
  X <- matrix(rnorm(2 * p * p), 2 * p, p); S <- crossprod(X) / (2 * p)
  on.exit(setTimeLimit())
  expect_error({
    setTimeLimit(elapsed = 0.5, transient = TRUE)
    graphical_lasso(S, rho = 1e-4, thr = 1e-12, maxit = 1e6)
  }, "time limit")
  setTimeLimit()
  ## the session is still usable afterwards
  expect_true(graphical_lasso(diag(2), 0.1)$converged)
})

###############################################################################
## Equivalence with glassoFast, whose supports along a penalty path downstream
## analyses may rely on.
###############################################################################

test_that("the result matches glassoFast on ordinary matrices, all along the path", {
  skip_if_not_installed("glassoFast")
  set.seed(7)
  for (k in 1:6) {
    X <- matrix(rnorm(50 * 20), 50, 20); S <- crossprod(X) / 50
    if (k %% 2 == 0) S <- S * 10^(-k)
    for (rho in penalty_path(S)) {
      ref <- glassoFast::glassoFast(S, rho = rho)$wi
      new <- graphical_lasso(S, rho)$wi
      expect_identical(support(new), support(ref))
      expect_equal(new, ref, tolerance = 1e-12)
    }
  }
})

test_that("the result matches glassoFast with weights and an unpenalized diagonal", {
  skip_if_not_installed("glassoFast")
  set.seed(11)
  X <- matrix(rnorm(60 * 15), 60, 15); S <- crossprod(X) / 60
  W <- matrix(runif(15 * 15, .5, 2), 15); W <- (W + t(W)) / 2; diag(W) <- 0
  for (lambda in c(0.3, 0.1, 0.03, 0.01)) {
    ref <- glassoFast::glassoFast(S, rho = lambda * W)
    new <- graphical_lasso(S, lambda * W)
    expect_identical(support(new$wi), support(ref$wi))
    expect_equal(new$wi, ref$wi, tolerance = 1e-12)
    expect_equal(new$w,  ref$w,  tolerance = 1e-12)
  }
})

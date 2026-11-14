context("test-pln")
require(purrr)

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance[1:20, 1:5], trichoptera$Covariate[1:20, ])

test_that("PLN: Check that PLN is running and robust",  {

  expect_is(PLN(Abundance ~ 1, data = trichoptera), "PLNfit")

  expect_is(PLN(Abundance ~ 0, data = trichoptera), "PLNfit")

  expect_is(PLN(trichoptera$Abundance ~ 1), "PLNfit")

  expect_equal(PLN(trichoptera$Abundance ~ 1 + trichoptera$Wind)$fitted,
               PLN(Abundance ~ Wind, data = trichoptera)$fitted)
})

test_that("PLN: Check consistency of initialization - fully parametrized covariance",  {

  ## use default initialization (LM), run enough EM iters to ensure convergence
  ctrl_cv <- PLN_param(trace = 0, config_optim = list(maxit_em = 500))
  model1 <- PLN(Abundance ~ 1, data = trichoptera, control = ctrl_cv)

  ## initialization with the previous fit: should converge to the same point
  model2 <- PLN(Abundance ~ 1, data = trichoptera, control = PLN_param(inception = model1, trace = 0))

  expect_equal(model2$loglik   , model1$loglik   , tolerance = 0.1)
  tol <- 1e-2
  expect_lt(sum((model2$model_par$Theta - model1$model_par$Theta)^2), tol)
  expect_lt(sum((model2$model_par$Sigma - model1$model_par$Sigma)^2), tol)
  expect_lt(sum((model2$var_par$M - model1$var_par$M)^2), tol)
  expect_lt(sum((model2$var_par$S2 - model1$var_par$S2)^2), tol)

})

test_that("PLN: Check consistency of initialization - diagonal covariance",  {

  ## use default initialization, run enough EM iters to ensure convergence
  ctrl_cv <- PLN_param(trace = 0, covariance = "diagonal", config_optim = list(maxit_em = 500))
  model1 <- PLN(Abundance ~ 1, data = trichoptera, control = ctrl_cv)

  ## initialization with the previous fit: should converge to the same point
  model2 <- PLN(Abundance ~ 1, data = trichoptera, control = PLN_param(inception = model1, trace = 0, covariance = "diagonal"))

  expect_equal(model2$loglik   , model1$loglik   , tolerance = 0.1)
  tol <- 1e-2
  expect_lt(sum((model2$model_par$Theta - model1$model_par$Theta)^2), tol)
  expect_lt(sum((model2$model_par$Sigma - model1$model_par$Sigma)^2), tol)
  expect_lt(sum((model2$var_par$M - model1$var_par$M)^2), tol)
  expect_lt(sum((model2$var_par$S2 - model1$var_par$S2)^2), tol)
})

test_that("PLN: Check consistency of observation weights - fully parameterized covariance",  {
  tol <- 1e-2

  ## no weights
  model1 <- PLN(Abundance ~ 1, data = trichoptera, control = PLN_param(trace = 0))

  ## equivalent weigths
  expect_output(model2 <- PLN(Abundance ~ 1, data = trichoptera, weights = rep(1.0, nrow(trichoptera))),
                paste("\n Initialization...",
                      "Adjusting a full covariance PLN model with nlopt optimizer",
                      "Post-treatments...",
                      "DONE!", sep = "\n "), fixed = TRUE)

  expect_equal(model2$loglik, model1$loglik, tolerance = tol)


})

test_that("PLN: Check consistency of observation weights - diagonal covariance",  {
  tol <- 1e-2

  ## no weights
  model1 <- PLN(Abundance ~ 1, data = trichoptera, control = PLN_param(covariance = "diagonal", trace = 0))

  ## equivalent weights
  model2 <- PLN(Abundance ~ 1, data = trichoptera, weights = rep(1.0, nrow(trichoptera)), control = PLN_param(covariance = "diagonal", trace = 0))

  expect_equal(model2$loglik   , model1$loglik   , tolerance = tol)
})

test_that("PLN: Check consistency of observation weights - spherical covariance",  {
  tol <- 1e-2

  ## no weights
  model1 <- PLN(Abundance ~ 1, data = trichoptera, control = PLN_param(covariance = "spherical", trace = 0))

  ## equivalent weights
  model2 <- PLN(Abundance ~ 1, data = trichoptera, weights = rep(1.0, nrow(trichoptera)), control = PLN_param(covariance = "spherical", trace = 0))

  expect_equal(model2$loglik   , model1$loglik   , tolerance = tol)
})

test_that("PLN: Routine comparison between the different covariance models",  {
  model_full      <- PLN(Abundance ~ 1, data = trichoptera, control = PLN_param(covariance = "full"     , trace = 0))
  model_diagonal  <- PLN(Abundance ~ 1, data = trichoptera, control = PLN_param(covariance = "diagonal" , trace = 0))
  model_spherical <- PLN(Abundance ~ 1, data = trichoptera, control = PLN_param(covariance = "spherical", trace = 0))
  expect_gte(model_full$loglik  , model_diagonal$loglik)
  expect_gte(model_diagonal$loglik, model_spherical$loglik)
})

test_that("PLN is working with a single variable data matrix",  {
  Y <- matrix(rpois(10, exp(0.5)), ncol = 1)
  colnames(Y) <- "Y"
  expect_is(PLN(Y ~ 1), "PLNfit")

  Y <- matrix(rpois(10, exp(0.5)), ncol = 1)
  expect_is(PLN(Y ~ 1), "PLNfit")
})

test_that("PLN is working with unnamed data matrix",  {
  n = 10; d = 3; p = 10
  Y <- matrix(rpois(n*p, 1), n, p)
  X <- matrix(rnorm(n*d), n, d)
  expect_is(PLN(Y ~ X), "PLNfit")
})

 test_that("PLN is working with different optimization algorithm in NLopt",  {

    MMA    <- PLN(Abundance ~ 1, data = trichoptera, control = PLN_param(config_optim = list(algorithm = "MMA")))
    CCSAQ  <- PLN(Abundance ~ 1, data = trichoptera, control = PLN_param(config_optim = list(algorithm = "CCSAQ")))
    LBFGS  <- PLN(Abundance ~ 1, data = trichoptera, control = PLN_param(config_optim = list(algorithm = "LBFGS")))

    expect_equal(MMA$loglik, CCSAQ$loglik, tolerance = 1e-1) ## Almost equivalent, CCSAQ faster
    expect_equal(MMA$loglik, LBFGS$loglik, tolerance = 1e-1)

    expect_error(PLN(Abundance ~ 1, data = trichoptera, control = PLN_param(backend = "nlopt", config_optim = list(algorithm = "nawak"))))
 })

test_that("PLN: Check that univariate PLN models works, with matrix of numeric format",  {
  expect_no_error(uniPLN <- PLN(Abundance[,1,drop=FALSE] ~ 1, data = trichoptera))
  expect_no_error(uniPLN <- PLN(Abundance[,1] ~ 1, data = trichoptera))
   y <- trichoptera$Abundance[,1]
   expect_no_error(uniPLN <- PLN(y ~ 1))
})

test_that("PLN: Check that all univariate PLN models are equivalent with the multivariate diagonal case",  {

  p <- ncol(trichoptera$Abundance)
  Offset <- trichoptera$Offset

  univariate_full <- lapply(1:p, function(j) {
    Abundance <- trichoptera$Abundance[, j, drop = FALSE]
    PLN(Abundance ~ 1 + offset(log(Offset)), control = PLN_param(trace = 0, config_optim = list(algorithm = "LBFGS")))
  })

  univariate_diagonal <- lapply(1:p, function(j) {
    Abundance <- trichoptera$Abundance[, j, drop = FALSE]
    PLN(Abundance ~ 1 + offset(log(Offset)), control = PLN_param(covariance = "diagonal", trace = 0))
  })

  univariate_spherical <- lapply(1:p, function(j) {
    Abundance <- trichoptera$Abundance[, j, drop = FALSE]
    PLN(Abundance ~ 1 + offset(log(Offset)), control = PLN_param(covariance = "spherical", trace = 0))
  })

  multivariate_diagonal <-
    PLN(Abundance ~ 1 + offset(log(Offset)), data = trichoptera, control = PLN_param(covariance = "diagonal", trace = 0))

  expect_univariate_equivalence(univariate_full, univariate_diagonal, univariate_spherical,
                                multivariate_diagonal)

})

test_that("PLN: singular covariate model matrix is handled gracefully",  {
  n = 10; d = 1; p = 10
  Y <- matrix(rpois(n*p, 1), n, p)

  # Artifically build singular model matrix with continuous covariates
  x <- matrix(rnorm(n*d), n, d)
  X_singular <- cbind(x, x)

  # Artificially build correlated factors (discrete covariates)
  f1 <- gl(2, n/2, labels = c("1.1", "1.2"))
  f2 <- gl(2, n/2, labels = c("2.1", "2.2"))

  # singular.ok = TRUE: dropped coefficients are zeroed, PLN does not throw
  expect_no_error(PLN(Y ~ X_singular))
  expect_no_error(PLN(Y ~ f1 + f2))

})

test_that("PLN builtin: SQUAREM, Newton-CG VE step and tests on the bound with its constants", {
  tr <- trichoptera
  ctrl <- function(...) PLN_param(trace = 0, backend = "builtin", config_optim = list(...))
  ref <- PLN(Abundance ~ 1 + Wind, data = tr, control = ctrl())
  ## the options are off by default
  expect_equal(PLN(Abundance ~ 1 + Wind, data = tr,
                   control = ctrl(objective_constant = FALSE, accelerate = FALSE, cg_max = 1L))$loglik, ref$loglik)
  ## accelerated EM loop, loose VE steps, tests on the bound: at least as good as the default
  for (cg in c(1L, 5L)) {
    fit <- PLN(Abundance ~ 1 + Wind, data = tr,
               control = ctrl(objective_constant = TRUE, accelerate = TRUE, ftol_in = 1e-6,
                              ftol_em = 1e-8, maxit_em = 2000L, cg_max = cg))
    expect_gte(fit$loglik, ref$loglik - 1e-3)
  }
  for (cov in c("diagonal", "spherical")) {
    ref <- PLN(Abundance ~ 1, data = tr, control = PLN_param(trace = 0, backend = "builtin", covariance = cov))
    fit <- PLN(Abundance ~ 1, data = tr,
               control = PLN_param(trace = 0, backend = "builtin", covariance = cov,
                                   config_optim = list(objective_constant = TRUE, accelerate = TRUE, ftol_in = 1e-6,
                                                       ftol_em = 1e-8, maxit_em = 2000L, cg_max = 5L)))
    expect_gte(fit$loglik, ref$loglik - 1e-3)
  }

  ## Newton-CG VE step, B and Omega fixed: same optimum as the cell-wise step, fewer iterations
  Y <- as.matrix(tr$Abundance); n <- nrow(Y)
  data <- list(Y = Y, X = matrix(1, n, 1), O = matrix(0, n, ncol(Y)), w = rep(1, n))
  full <- PLN(Abundance ~ 1, data = tr, control = PLN_param(trace = 0, backend = "builtin"))
  params <- list(M = log(1 + Y), S2 = matrix(.1, n, ncol(Y)),
                 B = matrix(colMeans(log(1 + Y)), 1), Omega = as.matrix(full$model_par$Omega))
  ve <- lapply(c(cellwise = 1L, cg = 5L), function(k)
    PLNmodels:::builtin_optimize_vestep_full(data, params, list(cg_max = k, ftol_in = 1e-12, maxeval = 10000L)))
  expect_equal(sum(ve$cg$Ji), sum(ve$cellwise$Ji), tolerance = 1e-6)
  expect_lt(length(ve$cg$monitoring$objective), length(ve$cellwise$monitoring$objective))
})

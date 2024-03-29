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

  ## use default initialization (LM)
  model1 <- PLN(Abundance ~ 1, data = trichoptera, control = PLN_param(trace = 0))

  ## initialization with the previous fit
  model2 <- PLN(Abundance ~ 1, data = trichoptera, control = PLN_param(inception = model1, trace = 0))

  expect_equal(model2$loglik   , model1$loglik   , tolerance = 0.1)
  tol <- 1e-2
  expect_lt(sum((model2$model_par$Theta - model1$model_par$Theta)^2), tol)
  expect_lt(sum((model2$model_par$Sigma - model1$model_par$Sigma)^2), tol)
  expect_lt(sum((model2$var_par$M - model1$var_par$M)^2), tol)
  expect_lt(sum((model2$var_par$S2 - model1$var_par$S2)^2), tol)

})

test_that("PLN: Check consistency of initialization - diagonal covariance",  {

  ## use default initialization (GLM)
  model1 <- PLN(Abundance ~ 1, data = trichoptera, control = PLN_param(trace = 0, covariance = "diagonal"))

  ## initialization with the previous fit
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
  model1 <- PLN(Abundance ~ 1, data = trichoptera, control = PLN_param(covariance = "spherical", trace = 0))

  ## equivalent weigths
  model2 <- PLN(Abundance ~ 1, data = trichoptera, weights = rep(1.0, nrow(trichoptera)), control = PLN_param(covariance = "spherical", trace = 0))

  expect_equal(model2$loglik   , model1$loglik   , tolerance = tol)
})

test_that("PLN: Check consistency of observation weights - spherical covariance",  {
  tol <- 1e-2

  ## no weights
  model1 <- PLN(Abundance ~ 1, data = trichoptera, control = PLN_param(covariance = "spherical", trace = 0))

  ## equivalent weigths
  model2 <- PLN(Abundance ~ 1, data = trichoptera, weights = rep(1.0, nrow(trichoptera)), control = PLN_param(covariance = "spherical", trace = 0))
  model3 <- PLN(Abundance ~ 1, data = trichoptera, weights = runif(nrow(trichoptera)), control = PLN_param(covariance = "spherical", trace = 0))

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

    expect_error(PLN(Abundance ~ 1, data = trichoptera, control = PLN_param(config_optim = list(algorithm = "nawak"))))
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

  expect_true(all.equal(
    map_dbl(univariate_spherical, "nb_param"),
    map_dbl(univariate_full     , "nb_param")
  ))

  expect_true(all.equal(
    map_dbl(univariate_spherical, "nb_param"),
    map_dbl(univariate_diagonal , "nb_param")
  ))
  expect_true(all.equal(
    map_dbl(univariate_full     , "nb_param"),
    map_dbl(univariate_diagonal , "nb_param")
  ))

  expect_true(all.equal(
    map_dbl(univariate_full, "loglik") %>% sum(),
    multivariate_diagonal$loglik, tolerance = 1e-2)
  )

  expect_true(all.equal(
    map_dbl(univariate_diagonal, "loglik") %>% sum(),
    multivariate_diagonal$loglik, tolerance = 1e-2)
  )

   expect_true(all.equal(
    map_dbl(univariate_spherical, "loglik") %>% sum(),
    multivariate_diagonal$loglik, tolerance = 1e-2)
  )

  expect_true(all.equal(
    map(univariate_spherical, sigma) %>% map_dbl(as.double),
    map(univariate_diagonal , sigma) %>% map_dbl(as.double), tolerance = .25
  ))

  expect_true(all.equal(
    map(univariate_spherical, sigma) %>% map_dbl(as.double),
    map(univariate_full , sigma) %>% map_dbl(as.double), tolerance = .25
  ))

  expect_true(all.equal(
    map(univariate_diagonal, sigma) %>% map_dbl(as.double),
    map(univariate_full , sigma) %>% map_dbl(as.double), tolerance = .25
  ))

})

context("test-zipln")
require(purrr)

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance[1:20, 1:5], trichoptera$Covariate[1:20, ])

test_that("ZIPLN: Check that ZIPLN is running and robust",  {

  expect_is(zi_single <- ZIPLN(Abundance ~ 1, data = trichoptera), "ZIPLNfit")
  expect_is(zi_row <- ZIPLN(Abundance ~ 1, data = trichoptera, zi = "row"), "ZIPLNfit")
  expect_is(zi_col <- ZIPLN(Abundance ~ 1, data = trichoptera, zi = "col"), "ZIPLNfit")
  expect_is(zi_covar  <- ZIPLN(Abundance ~ 1 | Wind, data = trichoptera), "ZIPLNfit")
  expect_is(zi_covar  <- ZIPLN(Abundance ~ 0 + Wind | Wind, data = trichoptera), "ZIPLNfit")

  ## initialization could not work without any regressor...
  expect_error(ZIPLN(Abundance ~ 0, data = trichoptera))

  expect_is(ZIPLN(trichoptera$Abundance ~ 1), "ZIPLNfit")

  expect_equal(ZIPLN(trichoptera$Abundance ~ 1 + trichoptera$Wind)$fitted,
               ZIPLN(Abundance ~ Wind, data = trichoptera)$fitted, tolerance = 1e-3)

  expect_is(model_sparse <- ZIPLN(Abundance ~ 1, data = trichoptera, control = ZIPLN_param(penalty = 0.2, trace = 0)), "ZIPLNfit_sparse")
  expect_is(model_fixed <- ZIPLN(Abundance ~ 1, data = trichoptera, control = ZIPLN_param(Omega = zi_single$model_par$Omega, trace = 0)), "ZIPLNfit_fixed")

})

test_that("ZIPLN: Routine comparison between the different covariance models",  {
  model_full      <- ZIPLN(Abundance ~ 1, data = trichoptera, control = ZIPLN_param(covariance = "full"     , trace = 0))
  model_diagonal  <- ZIPLN(Abundance ~ 1, data = trichoptera, control = ZIPLN_param(covariance = "diagonal" , trace = 0))
  model_spherical <- ZIPLN(Abundance ~ 1, data = trichoptera, control = ZIPLN_param(covariance = "spherical", trace = 0))
  expect_gte(model_full$loglik  , model_diagonal$loglik)
  expect_gte(model_diagonal$loglik, model_spherical$loglik)
})

test_that("PLN is working with a single variable data matrix",  {
  Y <- matrix(rpois(10, exp(0.5)), ncol = 1)
  colnames(Y) <- "Y"
  expect_is(ZIPLN(Y ~ 1), "ZIPLNfit")

  Y <- matrix(rpois(10, exp(0.5)), ncol = 1)
  expect_is(ZIPLN(Y ~ 1), "ZIPLNfit")
})

test_that("PLN is working with unnamed data matrix",  {
  n = 10; d = 3; p = 10
  Y <- matrix(rpois(n*p, 1), n, p)
  X <- matrix(rnorm(n*d), n, d)
  expect_is(ZIPLN(Y ~ X), "ZIPLNfit")
})

 test_that("ZIPLN is working with different optimization algorithm in NLopt",  {

    MMA    <- ZIPLN(Abundance ~ 1, data = trichoptera, control = ZIPLN_param(config_optim = list(algorithm = "MMA")))
    CCSAQ  <- ZIPLN(Abundance ~ 1, data = trichoptera, control = ZIPLN_param(config_optim = list(algorithm = "CCSAQ")))
    LBFGS  <- ZIPLN(Abundance ~ 1, data = trichoptera, control = ZIPLN_param(config_optim = list(algorithm = "LBFGS")))

    expect_equal(MMA$loglik, CCSAQ$loglik, tolerance = 1e-1) ## Almost equivalent, CCSAQ faster

    expect_error(ZIPLN(Abundance ~ 1, data = trichoptera, control = ZIPLN_param(config_optim = list(algorithm = "nawak"))))
 })

test_that("ZIPLN is working with exact and variational inference for the conditional distribution of the ZI component",  {

   approx <- ZIPLN(Abundance ~ 1, data = trichoptera, control = ZIPLN_param(config_optim = list(approx_ZI = TRUE)))
   exact  <- ZIPLN(Abundance ~ 1, data = trichoptera, control = ZIPLN_param(config_optim = list(approx_ZI = FALSE)))

   expect_equal(approx$loglik, exact$loglik, tolerance = 1e-1) ## Almost equivalent
   expect_equal(approx$model_par$B, exact$model_par$B, tolerance = 1e-1) ## Almost equivalent
   expect_equal(approx$model_par$Sigma, exact$model_par$Sigma, tolerance = 1e-1) ## Almost equivalent

})

test_that("ZIPLN: Check that univariate ZIPLN models works, with matrix of numeric format",  {
  expect_no_error(uniZIPLN <- ZIPLN(Abundance[,1,drop=FALSE] ~ 1, data = trichoptera))
  expect_no_error(uniZIPLN <- ZIPLN(Abundance[,1] ~ 1, data = trichoptera))
   y <- trichoptera$Abundance[,1]
   expect_no_error(uniZIPLN <- ZIPLN(y ~ 1))
})

test_that("ZIPLN: Check that all univariate ZIPLN models are equivalent with the multivariate diagonal case",  {

  p <- ncol(trichoptera$Abundance)
  Offset <- trichoptera$Offset
  Wind <- trichoptera$Wind

  univariate_full <- lapply(1:p, function(j) {
    Abundance <- trichoptera$Abundance[, j, drop = FALSE]
    ZIPLN(Abundance ~ 1 + offset(log(Offset)) | Wind, control = ZIPLN_param(trace = 0))
  })

  univariate_diagonal <- lapply(1:p, function(j) {
    Abundance <- trichoptera$Abundance[, j, drop = FALSE]
    ZIPLN(Abundance ~ 1 + offset(log(Offset)) | Wind, control = ZIPLN_param(covariance = "diagonal", trace = 0))
  })

  univariate_spherical <- lapply(1:p, function(j) {
    Abundance <- trichoptera$Abundance[, j, drop = FALSE]
    ZIPLN(Abundance ~ 1 + offset(log(Offset)) | Wind, control = ZIPLN_param(covariance = "spherical", trace = 0))
  })

  multivariate_diagonal <-
    ZIPLN(Abundance ~ 1 + offset(log(Offset)) | Wind, data = trichoptera, control = ZIPLN_param(covariance = "diagonal", trace = 0))

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

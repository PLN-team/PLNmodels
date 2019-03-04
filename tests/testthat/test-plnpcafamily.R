context("test-plnpcafamily")

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

test_that("PLNPCAfamily: main function, field access and methods", {

  models <- PLNPCA(Abundance ~ 1, data = trichoptera)
  expect_is(models, "PLNPCAfamily")

  X <- model.matrix(Abundance ~ 1, data = trichoptera)
  Y <- as.matrix(trichoptera$Abundance)
  O <- matrix(0, nrow = nrow(Y), ncol = ncol(Y))
  w <- rep(1, nrow(Y))

  ## extract the data matrices and weights
  ctrl_init <- PLNmodels:::PLN_param(list(), nrow(Y), ncol(Y), ncol(X), weighted = FALSE)
  ctrl_main <- PLNmodels:::PLNPCA_param(list(), weighted = FALSE)

  ## instantiate
  myPLN <- PLNmodels:::PLNPCAfamily$new(1:5, Y, X, O, w, Abundance ~ 1, ctrl_init)

  ## optimize
  myPLN$optimize(ctrl_main)

  ## post-treatment
  myPLN$postTreatment()

  expect_equivalent(myPLN, models)

  ## S3 methods
  expect_true(PLNmodels:::isPLNPCAfamily(myPLN))
  expect_is(plot(myPLN), "ggplot")
  expect_is(plot(myPLN, map="individual"), "ggplot")
  expect_is(plot(myPLN, map="variable"), "ggplot")
  expect_is(getBestModel(myPLN), "PLNPCAfit")
  expect_is(getModel(myPLN, myPLN$ranks[1]), "PLNPCAfit")
})

test_that("PLNPCA is fast on low ranks", {

  n <- 100
  p <- 1000
  lambda <- exp(rnorm(n * p))
  Y <- matrix(rpois(n * p, lambda), n, p)

  models <- PLNPCA(Y ~ 1, ranks = 1:3)
  expect_is(models, "PLNPCAfamily")
})

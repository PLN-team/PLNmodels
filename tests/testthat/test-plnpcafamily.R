context("test-plnpcafamily")

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

test_that("PLNPCAfamily methods", {

  models <- PLNPCA(Abundance ~ 1, data = trichoptera)

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
  ## add some
})


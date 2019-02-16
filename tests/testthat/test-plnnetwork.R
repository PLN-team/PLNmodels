context("test-plnnetwork")

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

test_that("PLNnetwork runs", {

  models <- PLNnetwork(Abundance ~ 1, data = trichoptera)
  expect_is(models, "PLNnetworkfamily")

  ## this is more PLNPCAnetworkfit testing
  expect_equal(getBestModel(models), getBestModel(models, "BIC"))

})

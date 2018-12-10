context("test-plnnetwork")

test_that("PLNnetwork runs", {

  models <- PLNnetwork(Abundance ~ 1, data = trichoptera)
  expect_is(models, "PLNnetworkfamily")

  ## this is more PLNPCAnetworkfit testing
  expect_equal(getBestModel(models), getBestModel(models, "BIC"))

})

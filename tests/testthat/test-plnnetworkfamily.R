context("test-plnnetworkfamily")

data(trichoptera)

test_that("PLNnetwork methods", {

  models <- PLNnetwork(Abundance ~ 1, data = trichoptera)
  expect_equal(getBestModel(models), getBestModel(models, "BIC"))

  ## add some
})

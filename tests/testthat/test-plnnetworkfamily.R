context("test-plnnetworkfamily")

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

data(mollusk)
mollusc <- prepare_data(mollusk$Abundance, mollusk$Covariate)

test_that("PLNnetwork methods", {

  models <- PLNnetwork(Abundance ~ 1, data = trichoptera)
  expect_equal(getBestModel(models), getBestModel(models, "BIC"))

  ## add some
})

test_that("PLNnetwork computes the stability path only once.", {
  ## Compute network and stability selection once
  nets <- PLNnetwork(Abundance ~ 0 + site + offset(log(Offset)),
                     data = mollusc)
  set.seed(1077)
  subs <- replicate(2,
                    sample.int(nrow(mollusc), size = nrow(mollusc)/2),
                    simplify = FALSE)
  stability_selection(nets, subsamples = subs)
  ## Stability_path has correct dimensions
  p <- getModel(nets, index = 1)$p
  expect_equal(dim(nets$stability_path),
               c(length(nets$penalties) * p*(p-1)/2L, 5))
  ## try to compute it again
  expect_message(stability_selection(nets),
                 "Previous stability selection detected. Use \"force = TRUE\" to recompute it.")

})

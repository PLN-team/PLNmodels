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

test_that("PLNnetwork computes a the stability path only once.", {

  ## Compute network and stability selection once
  nets <- PLNnetwork(Abundance ~ 0 + site + offset(log(Offset)),
                     data = mollusc)
  set.seed(1077)
  subs <- replicate(2,
                    sample.int(nrow(mollusc), size = round(nrow(mollusc)/2)),
                    simplify = FALSE)
  stability_selection(nets, subsamples = subs)
  ## try to compute it again
  expect_message(stability_selection(nets),
                 "Previous stability selection detected. Use \"force = TRUE\" to recompute it.")

})

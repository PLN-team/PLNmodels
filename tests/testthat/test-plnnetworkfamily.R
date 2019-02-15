context("test-plnnetworkfamily")

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

data(mollusk)
mollusc <- suppressWarnings(
  prepare_data(mollusk$Abundance, mollusk$Covariate)
)

test_that("PLNnetwork methods", {

  models <- PLNnetwork(Abundance ~ 1, data = trichoptera)
  expect_equal(getBestModel(models), getBestModel(models, "BIC"))

  ## add some
})

test_that("PLNnetwork computes the stability path only once.", {
  ## Compute network and stability selection once
  nets <- PLNnetwork(Abundance ~ 0 + site + offset(log(Offset)),
                     data = mollusc)
  ## extract_probs fails if stability selection has not been performed.
  expect_error(extract_probs(nets),
               "Please perform stability selection using stability_selection(Robject) first", fixed = TRUE)
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
  ## extracts the inclusion frequencies
  expect_equal(dim(extract_probs(nets, index = 1, format = "matrix")),
               c(p, p))
  expect_length(extract_probs(nets, index = 1, format = "vector"),
               p*(p-1)/2)
})

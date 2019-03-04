context("test-plnnetworkfamily")

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

data(mollusk)
mollusc <- suppressWarnings(
  prepare_data(mollusk$Abundance, mollusk$Covariate)
)

test_that("PLNnetwork: main function, fields access and methods", {

  models <- PLNnetwork(Abundance ~ 1, data = trichoptera)
  expect_equal(getBestModel(models), getBestModel(models, "BIC"))

  X <- model.matrix(Abundance ~ 1, data = trichoptera)
  Y <- as.matrix(trichoptera$Abundance)
  O <- matrix(0, nrow(Y),ncol(Y))
  w <- rep(1, nrow(Y))

  ## extract the data matrices and weights
  ctrl_init <- PLNmodels:::PLN_param(list(), nrow(Y), ncol(Y), ncol(X), weighted = FALSE)
  ctrl_init$trace <- 0; ctrl_init$nPenalties <- 30; ctrl_init$min.ratio   <- .1
  ctrl_main <- PLNmodels:::PLNnetwork_param(list(), nrow(Y), ncol(Y),ncol(X), weighted = FALSE)

  ## instantiate
  myPLN <- PLNmodels:::PLNnetworkfamily$new(NULL, Y, X, O, w, Abundance ~ 1, ctrl_init)

  ## optimize
  myPLN$optimize(ctrl_main)

  ## post-treatment
  myPLN$postTreatment()

  expect_equivalent(myPLN, models)

  ## S3 methods
  expect_true(PLNmodels:::isPLNnetworkfamily(myPLN))
  expect_is(plot(myPLN), "ggplot")
  expect_is(getBestModel(myPLN), "PLNnetworkfit")
  expect_is(getModel(myPLN, myPLN$penalties[1]), "PLNnetworkfit")

  ## Field access
  expect_true(all(myPLN$penalties > 0))
  expect_null(myPLN$stability_path)
  expect_true(anyNA(myPLN$stability))

  ## Other R6 methods
  expect_true(is.data.frame(myPLN$coefficient_path()))
  subs <- replicate(2,
                    sample.int(nrow(trichoptera), size = nrow(trichoptera)/2),
                    simplify = FALSE)
  myPLN$stability_selection(subsamples = subs)
  expect_true(!is.null(myPLN$stability_path))
  expect_true(inherits(myPLN$plot(), "ggplot"))
  expect_true(inherits(myPLN$plot_objective(), "ggplot"))
  expect_true(inherits(myPLN$plot_stars(), "ggplot"))
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

context("test-plnnetworkfamily")

data(trichoptera)
## use a subset t save some time
trichoptera <- prepare_data(trichoptera$Abundance[1:20, 1:5], trichoptera$Covariate[1:20, ])

models <- PLNnetwork(Abundance ~ 1, data = trichoptera)

test_that("PLNnetwork: main function, fields access and methods", {

  expect_equal(getBestModel(models), getBestModel(models, "BIC"))

  X <- model.matrix(Abundance ~ 1, data = trichoptera)
  Y <- as.matrix(trichoptera$Abundance)
  O <- matrix(0, nrow(Y),ncol(Y))
  w <- rep(1, nrow(Y))

  ## extract the data matrices and weights
  ctrl <- PLNmodels:::PLNnetwork_param(trace = 0)

  ## instantiate
  myPLN <- PLNmodels:::PLNnetworkfamily$new(NULL, Y, X, O, w, Abundance ~ 1, ctrl)

  ## optimize
  myPLN$optimize(ctrl$config_optim)

  ## post-treatment
  config_post <- PLNmodels:::config_post_default_PLNnetwork
  config_post$trace <- 0
  myPLN$postTreatment(config_post)

  expect_equal(myPLN$criteria$BIC, models$criteria$BIC)

  ## S3 methods
  expect_true(PLNmodels:::isPLNnetworkfamily(myPLN))
  expect_is(plot(myPLN), "ggplot")
  expect_is(plot(myPLN, reverse = TRUE), "ggplot")
  expect_is(plot(myPLN, type = "diagnostic"), "ggplot")
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
  expect_is(plot(myPLN, type = "stability"), "ggplot")
  expect_true(!is.null(myPLN$stability_path))
  expect_true(inherits(myPLN$plot(), "ggplot"))
  expect_true(inherits(myPLN$plot_objective(), "ggplot"))
  expect_true(inherits(myPLN$plot_stars(), "ggplot"))
})

test_that("PLNnetwork computes the stability path only once.", {

  ## extract_probs fails if stability selection has not been performed.
  expect_error(extract_probs(models),
               "Please perform stability selection using stability_selection(Robject) first", fixed = TRUE)
  set.seed(1077)
  subs <- replicate(2,
                    sample.int(nrow(trichoptera), size = nrow(trichoptera)/2),
                    simplify = FALSE)
  stability_selection(models, subsamples = subs)
  ## Stability_path has correct dimensions
  p <- getModel(models, index = 1)$p
  expect_equal(dim(models$stability_path),
               c(length(models$penalties) * p*(p-1)/2L, 5))
  ## try to compute it again
  expect_message(stability_selection(models),
                 "Previous stability selection detected. Use \"force = TRUE\" to recompute it.")
  ## extracts the inclusion frequencies
  expect_equal(dim(extract_probs(models, index = 1, format = "matrix")),
               c(p, p))
  expect_length(extract_probs(models, index = 1, format = "vector"),
               p*(p-1)/2)
})

test_that("PLNnetwork: matrix of penalties work", {

  p <- ncol(trichoptera$Abundance)
  W <- diag(1, p, p)
  W[upper.tri(W)] <- runif(p*(p-1)/2, min = 1, max = 5)
  W[lower.tri(W)] <- t(W)[lower.tri(W)]
  myPLN <- PLNnetwork(Abundance ~ 1, data = trichoptera, control = PLNnetwork_param(penalty_weights = W))

  ## S3 methods
  expect_true(PLNmodels:::isPLNnetworkfamily(myPLN))
  expect_is(plot(myPLN), "ggplot")
  expect_is(plot(myPLN, reverse = TRUE), "ggplot")
  expect_is(plot(myPLN, type = "diagnostic"), "ggplot")
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
  myPLN$stability_selection(subsamples = subs, control = PLNnetwork_param(penalty_weights = W))
  expect_is(plot(myPLN, type = "stability"), "ggplot")
  expect_true(!is.null(myPLN$stability_path))
  expect_true(inherits(myPLN$plot(), "ggplot"))
  expect_true(inherits(myPLN$plot_objective(), "ggplot"))
  expect_true(inherits(myPLN$plot_stars(), "ggplot"))

  ## missspecification of penlaty weights should induce errors
  ## not symmetric
  W <- diag(1, p, p)
  W[upper.tri(W)] <- runif(p*(p-1)/2, min = 1, max = 5)
  expect_error(PLNnetwork(Abundance ~ 1, data = trichoptera, control = PLNnetwork_param(penalty_weights = W)))

  ## not square
  W <- matrix(1, p + 1, p)
  expect_error(PLNnetwork(Abundance ~ 1, data = trichoptera, control = PLNnetwork_param(penalty_weights = W)))

  ## not-positive entries
  W <- matrix(0, p, p)
  expect_error(PLNnetwork(Abundance ~ 1, data = trichoptera, control = PLNnetwork_param(penalty_weights = W)))

})

test_that("PLNnetwork: list of matrices of penalties work", {

  p <- ncol(trichoptera$Abundance)
  W <- diag(1, p, p)
  W[upper.tri(W)] <- runif(p*(p-1)/2, min = 1, max = 5)
  W[lower.tri(W)] <- t(W)[lower.tri(W)]
  list_W <- lapply(seq(1, 1e-2, len = 30), function(rho) rho * W)

  myPLN <- PLNnetwork(Abundance ~ 1, data = trichoptera, control = PLNnetwork_param(penalty_weights = list_W))

  ## S3 methods
  expect_true(PLNmodels:::isPLNnetworkfamily(myPLN))
  expect_is(plot(myPLN), "ggplot")
  expect_is(plot(myPLN, reverse = TRUE), "ggplot")
  expect_is(plot(myPLN, type = "diagnostic"), "ggplot")
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
  myPLN$stability_selection(subsamples = subs, control = PLNnetwork_param(penalty_weights = W))
  expect_is(plot(myPLN, type = "stability"), "ggplot")
  expect_true(!is.null(myPLN$stability_path))
  expect_true(inherits(myPLN$plot(), "ggplot"))
  expect_true(inherits(myPLN$plot_objective(), "ggplot"))
  expect_true(inherits(myPLN$plot_stars(), "ggplot"))

})

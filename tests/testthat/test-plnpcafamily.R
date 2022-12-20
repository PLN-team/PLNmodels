context("test-plnpcafamily")

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

test_that("PLNPCAfamily: main function, field access and methods", {

  ##  does not work because of future evaluating in random order the models
  # output <- "\n Initialization...\n\n Adjusting 5 PLN models for PCA analysis.\n Rank approximation = 1\n\t conservative convex separable approximation for gradient descent Rank approximation = 2\n\t conservative convex separable approximation for gradient descent Rank approximation = 3\n\t conservative convex separable approximation for gradient descent Rank approximation = 4\n\t conservative convex separable approximation for gradient descent Rank approximation = 5\n\t conservative convex separable approximation for gradient descent\n Post-treatments\n DONE!"
  #
  models <- PLNPCA(Abundance ~ 1, data = trichoptera,
                                  ranks = 1:5, control = PLNPCA_param(trace = 0))

  expect_is(models, "PLNPCAfamily")
  expect_is(plot(models), "ggplot")
  expect_is(plot(models, reverse = TRUE), "ggplot")
  expect_is(plot(models, map="individual"), "ggplot")
  expect_is(plot(models, map="variable"), "ggplot")
  expect_is(getBestModel(models), "PLNPCAfit")
  expect_is(getModel(models, models$ranks[1]), "PLNPCAfit")

  X <- model.matrix(Abundance ~ 1, data = trichoptera)
  Y <- as.matrix(trichoptera$Abundance)
  O <- matrix(0, nrow = nrow(Y), ncol = ncol(Y))
  w <- rep(1, nrow(Y))

  ## extract the data matrices and weights
  ctrl <-PLNPCA_param()

  ## instantiate
  myPLN <- PLNmodels:::PLNPCAfamily$new(1:5, Y, X, O, w, Abundance ~ 1, ctrl)

  ## optimize
  myPLN$optimize(ctrl$config_optim)

  ## post-treatment
  myPLN$postTreatment(ctrl$config_post)

  ## S3 methods
  expect_true(PLNmodels:::isPLNPCAfamily(myPLN))
  expect_is(plot(myPLN), "ggplot")
  expect_is(plot(myPLN, map="individual"), "ggplot")
  expect_is(plot(myPLN, map="variable"), "ggplot")
  expect_is(getBestModel(myPLN), "PLNPCAfit")
  expect_is(getModel(myPLN, myPLN$ranks[1]), "PLNPCAfit")

  ## test fail on Mac OS : fit is slightly different between OS X and Ubuntu
  ## probably due to the version of optimization libraries...
  ## Show method
#   expect_output(models$show(),
# "--------------------------------------------------------
# COLLECTION OF 5 POISSON LOGNORMAL MODELS
# --------------------------------------------------------
#  Task: Principal Component Analysis
# ========================================================
#  - Ranks considered: from 1 to 5
#  - Best model (greater BIC): rank = 5
#  - Best model (greater ICL): rank = 4",
#   fixed = TRUE)

})

# test_that("PLNPCA is fast on low ranks", {
#
#   n <- 100
#   p <- 1000
#   lambda <- exp(rnorm(n * p))
#   Y <- matrix(rpois(n * p, lambda), n, p)
#
#   models <- PLNPCA(Y ~ 1, ranks = 1:3)
#   expect_is(models, "PLNPCAfamily")
# })

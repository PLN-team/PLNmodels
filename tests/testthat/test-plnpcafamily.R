context("test-plnpcafamily")

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

test_that("PLNPCAfamily: main function, field access and methods", {

  output <- "\n Initialization...\n\n Adjusting 5 PLN models for PCA analysis.\n Rank approximation = 1\n\t conservative convex separable approximation for gradient descent Rank approximation = 2\n\t conservative convex separable approximation for gradient descent Rank approximation = 3\n\t conservative convex separable approximation for gradient descent Rank approximation = 4\n\t conservative convex separable approximation for gradient descent Rank approximation = 5\n\t conservative convex separable approximation for gradient descent\n Post-treatments\n DONE!"

  expect_output(models <- PLNPCA(Abundance ~ 1, data = trichoptera,
                                 ranks = 1:5, control_main = list(trace = 2)),
              output, fixed = TRUE)

  expect_is(models, "PLNPCAfamily")

  X <- model.matrix(Abundance ~ 1, data = trichoptera)
  Y <- as.matrix(trichoptera$Abundance)
  O <- matrix(0, nrow = nrow(Y), ncol = ncol(Y))
  w <- rep(1, nrow(Y))
  xlevels <- NULL

  ## extract the data matrices and weights
  ctrl_init <- PLNmodels:::PLN_param(list(), nrow(Y), ncol(Y), ncol(X))
  ctrl_main <- PLNmodels:::PLNPCA_param(list())

  ## instantiate
  myPLN <- PLNmodels:::PLNPCAfamily$new(1:5, Y, X, O, w, Abundance ~ 1, xlevels, ctrl_init)

  ## optimize
  myPLN$optimize(ctrl_main)

  ## post-treatment
  myPLN$postTreatment()

  expect_equivalent(myPLN, models)

  ## S3 methods
  expect_true(PLNmodels:::isPLNPCAfamily(myPLN))
  expect_is(plot(myPLN), "ggplot")
  expect_is(plot(myPLN, map="individual"), "ggplot")
  expect_is(plot(myPLN, map="variable"), "ggplot")
  expect_is(getBestModel(myPLN), "PLNPCAfit")
  expect_is(getModel(myPLN, myPLN$ranks[1]), "PLNPCAfit")

  ## test fail on Mac OS : fit is slightly different between OS X and Ubuntu
  ## probably due to the version of optimization libraries...
#   ## Show method
#   expect_output(models$show(),
# "--------------------------------------------------------
# COLLECTION OF 5 POISSON LOGNORMAL MODELS
# --------------------------------------------------------
#  Task: Principal Component Analysis
# ========================================================
#  - Ranks considered: from 1 to 5
#  - Best model (greater BIC): rank = 4 - R2 = 0.98
#  - Best model (greater ICL): rank = 4 - R2 = 0.98",
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

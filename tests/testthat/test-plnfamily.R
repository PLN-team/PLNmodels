context("test-plnfamily")

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

test_that("PLNfamily: main function, field access and methods", {

  X <- model.matrix(Abundance ~ Group, data = trichoptera)
  Y <- as.matrix(trichoptera$Abundance)
  O <- matrix(0, nrow = nrow(Y), ncol = ncol(Y))
  n <- nrow(Y); p <- ncol(Y); d <- ncol(X)

  ## extract the data matrices and weights
  ctrl <- PLNmodels:::PLN_param()

  ## Simple PLN models
  model1 <- PLN(Abundance ~ 1, data = trichoptera)
  model2 <- PLN(Abundance ~ Group, data = trichoptera)

  ## instantiate with unnamed matrices
  models <- PLNmodels:::PLNfamily$new(responses = unname(Y),
                                     covariates = unname(X),
                                     offsets = NULL,
                                     weights = NULL,
                                     control = ctrl)
  ## Set params (hacky)
  models$.__enclos_env__$private$params <- c(0.1, 3)
  models$models <- list(model1, model2)

  ## Check default names
  expect_equal(dimnames(models$responses),
               list(as.character(1:n), as.character(1:p)))
  expect_equal(dimnames(models$covariates),
               list(as.character(1:n), as.character(1:d)))

  ## Check print methods
  expect_equal(capture_output(models$show()),
               capture_output(models$print()))
  expect_output(models$show(),
"--------------------------------------------------------
COLLECTION OF 2 POISSON LOGNORMAL MODELS
--------------------------------------------------------", fixed = TRUE)

  ## Check extractor via parameter value
  expect_equal(models$getModel(0.1), model1)
  expect_equal(models$getModel(3), model2)
  expect_warning(res <- models$getModel(1),
                 paste("No such a model in the collection. Acceptable parameter values can be found via",
                       "$ranks (for PCA)",
                       "$clusters (for mixture models)",
                       "$penalties (for network)",
                       "Returning model with closest value. Requested: 1 , returned: 0.1",
                       sep = "\n"), fixed = TRUE)
  expect_equal(res, model1)

  ## Check convergence
# Commented to suppred purr warnings
#  expect_equal(class(models$convergence), "data.frame")
#  expect_named(models$convergence, as.character(c("param", "nb_param", names(model1$optim_par))))
})

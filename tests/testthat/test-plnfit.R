context("test-plnfit")

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

test_that("PLN fit: check classes, getters and field access",  {

  model <- PLN(Abundance ~ 1, data = trichoptera)

  expect_is(model, "PLNfit")

  expect_equal(model$n, nrow(trichoptera$Abundance))
  expect_equal(model$p, ncol(trichoptera$Abundance))
  expect_equal(model$d, 1)

  ## S3 methods: values
  expect_equal(coef(model), model$model_par$Theta)
  expect_equal(coef(model, type = "covariance"), sigma(model))
  expect_equal(sigma(model), model$model_par$Sigma)
  expect_equal(vcov(model), model$fisher$mat)

  ## S3 methods: class
  expect_equal(class(coef(model)), "matrix")
  expect_equal(class(sigma(model)), "matrix")
  expect_true(class(vcov(model)) == "dgCMatrix")

  ## S3 methods: dimensions
  expect_equal(dim(vcov(model)), c(model$d * model$p, model$d * model$p))

  ## S3 methods errors
  expect_error(standard_error(model, type = "louis"),
               "Standard errors were not computed using the louis approximation. Try another approximation scheme.")

})

test_that("PLN fit: Check prediction",  {

  model1 <- PLN(Abundance ~ 1, data = trichoptera, subset = 1:30)
  model2 <- PLN(Abundance ~ Pressure + Humidity, data = trichoptera, subset = 1:30)

  newdata <- trichoptera[31:49, ]
  newdata$Abundance <- NULL

  expect_gt(
    mean((trichoptera$Abundance[31:49, ] - predict(model1, newdata = newdata, type = "response"))^2),
    mean((trichoptera$Abundance[31:49, ] - predict(model2, newdata = newdata, type = "response"))^2)
  )
})

test_that("PLN fit: Check number of parameters",  {

  p <- ncol(trichoptera$Abundance)

  model <- PLN(Abundance ~ 1, data = trichoptera)
  expect_equal(model$nb_param, p*(p+1)/2 + p * 1)

  model <- PLN(Abundance ~ 1 + Wind, data = trichoptera)
  expect_equal(model$nb_param, p*(p+1)/2 + p * 2)

  model <- PLN(Abundance ~ Group + 0 , data = trichoptera)
  expect_equal(model$nb_param, p*(p+1)/2 + p * nlevels(trichoptera$Group))

  model <- PLN(Abundance ~ 1, data = trichoptera, control = list(covariance = "diagonal"))
  expect_equal(model$nb_param, p + p * 1)

  model <- PLN(Abundance ~ 1, data = trichoptera, control = list(covariance = "spherical"))
  expect_equal(model$nb_param, 1 + p * 1)

})

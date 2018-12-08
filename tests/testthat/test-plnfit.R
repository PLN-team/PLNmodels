context("test-plnfit")

data("trichoptera")

test_that("Check classes, getters and field access",  {

  model <- PLN(Abundance ~ 1, data = trichoptera)

  expect_is(model, "PLNfit")

  expect_equal(model$n, nrow(trichoptera$Abundance))
  expect_equal(model$p, ncol(trichoptera$Abundance))
  expect_equal(model$d, 1)

  expect_equal(coef(model), model$model_par$Theta)
  expect_equal(vcov(model), model$model_par$Sigma)

  expect_equal(class(coef(model)), "matrix")
  expect_equal(class(vcov(model)), "matrix")

})

test_that("Check prediction",  {

  model1 <- PLN(Abundance ~ 1, data = trichoptera, subset = 1:30)
  model2 <- PLN(Abundance ~ Pressure + Humidity, data = trichoptera, subset = 1:30)

  expect_gt(
    mean((trichoptera$Abundance[31:49, ] - predict(model1, trichoptera[31:49, ], type = "response"))^2),
    mean((trichoptera$Abundance[31:49, ] - predict(model2, trichoptera[31:49, ], type = "response"))^2)
  )
})

test_that("Check degrees of freedom",  {

  p <- ncol(trichoptera$Abundance)

  model <- PLN(Abundance ~ 1, data = trichoptera)
  expect_equal(model$degrees_freedom, p*(p+1)/2 + p * 1)

  model <- PLN(Abundance ~ 1 + Wind, data = trichoptera)
  expect_equal(model$degrees_freedom, p*(p+1)/2 + p * 2)

  model <- PLN(Abundance ~ Group + 0 , data = trichoptera)
  expect_equal(model$degrees_freedom, p*(p+1)/2 + p * nlevels(trichoptera$Group))

  model <- PLN(Abundance ~ 1, data = trichoptera, control = list(covariance = "diagonal"))
  expect_equal(model$degrees_freedom, p + p * 1)

  model <- PLN(Abundance ~ 1, data = trichoptera, control = list(covariance = "spherical"))
  expect_equal(model$degrees_freedom, 1 + p * 1)

})

context("test-ziplnfit")

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

test_that("ZIPLN fit: check classes, getters and field access",  {

  model <- ZIPLN(Abundance ~ 1, data = trichoptera, control = ZIPLN_param(trace = 1))
  expect_is(model, "ZIPLNfit")

  expect_equal(model$n, nrow(trichoptera$Abundance))
  expect_equal(model$p, ncol(trichoptera$Abundance))
  expect_equal(model$d, 1)
  expect_equal(model$d0, 0)

  ## S3 methods: values
  expect_equal(coef(model), model$model_par$B)
  expect_equal(coef(model, type = "count"), model$model_par$B)
  expect_equal(coef(model, type = "zero"), model$model_par$B0)
  expect_equal(coef(model, type = "precision"), model$model_par$Omega)
  expect_equal(coef(model, type = "covariance"), model$model_par$Sigma)
  expect_equal(sigma(model), model$model_par$Sigma)

  ## S3 methods: class
  expect_true(inherits(coef(model), "matrix"))
  expect_true(inherits(sigma(model), "matrix"))

  ## R6 bindings
  expect_is(model$latent, "matrix")
  expect_true(is.numeric(model$latent))
  expect_equal(dim(model$latent), c(model$n, model$p))

})

test_that("ZIPLN fit: check print message",  {

  expect_output(model <- ZIPLN(Abundance ~ 1, data = trichoptera))

  ## show and print are equivalent
  expect_equal(capture_output(model$show()),
               capture_output(model$print()))
})

test_that("PLN fit: Check prediction",  {

  model1     <- ZIPLN(Abundance ~ 1, data = trichoptera, subset = 1:30)
  model1_off <- ZIPLN(Abundance ~ 1 + offset(log(Offset)), data = trichoptera, subset = 1:30)
  model2     <- ZIPLN(Abundance ~ 1 + Pressure , data = trichoptera, subset = 1:30)
  model3     <- ZIPLN(Abundance ~ 1 + Pressure + offset(log(Offset)) | Wind , data = trichoptera, subset = 1:30)

  newdata <- trichoptera[31:49, ]
  # newdata$Abundance <- NULL

  pred1     <- predict(model1, newdata = newdata, type = "response")
  pred1_off <- predict(model1_off, newdata = newdata, type = "response")
  pred2     <- predict(model2, newdata = newdata, type = "response")
  pred2_def <- predict(model2, newdata = newdata, type = "deflated")
  pred2_ve  <- predict(model2, newdata = newdata, type = "response",
                      responses = newdata$Abundance)

  ## predict returns fitted values if no data is provided
  expect_equal(model2$predict(), model2$fitted)

  ## Adding covariates improves fit
  expect_gt(
    mean((newdata$Abundance - pred1)^2),
    mean((newdata$Abundance - pred2)^2)
  )

  ## Doing one VE step improves fit
  expect_gt(
    mean((newdata$Abundance - pred2)^2),
    mean((newdata$Abundance - pred2_ve)^2)
  )

  ## Removing zero-inflation leads to higher predicted values
  expect_gt(min(pred2_ve, pred2_def), 0)

  ## R6 methods
  ## without offsets, predictions should be the same for all samples
  expect_equal(unname(apply(pred1, 2, sd)), rep(0, ncol(pred1)))

  ## Unequal factor levels in train and prediction datasets
  suppressWarnings(
    toy_data <- prepare_data(
    counts     = matrix(c(1, 3, 1, 1), ncol = 1),
    covariates = data.frame(Cov = c("A", "B", "A", "A")),
    offset     = rep(1, 4))
  )
  model <- ZIPLN(Abundance ~ Cov + offset(log(Offset)), data = toy_data[1:2,])
  expect_length(predict(model, newdata = toy_data[3:4, ], type = "r"), 2L)
})

test_that("ZIPLN fit: Check number of parameters",  {

  p <- ncol(trichoptera$Abundance)

  model <- ZIPLN(Abundance ~ 1, data = trichoptera)
  expect_equal(model$nb_param, p*(p+1)/2 + p * 1 + 1)

  model <- ZIPLN(Abundance ~ 1 + Wind, data = trichoptera)
  expect_equal(model$nb_param, p*(p+1)/2 + p * 2 + 1)

  model <- ZIPLN(Abundance ~ Group + 0 , data = trichoptera)
  expect_equal(model$nb_param, p*(p+1)/2 + p * nlevels(trichoptera$Group) + 1)

  modelS <- ZIPLN(Abundance ~ 1, data = trichoptera, control = ZIPLN_param(covariance = "spherical"))
  expect_equal(modelS$nb_param, 1 + p * 1 + 1)
  expect_equal(modelS$vcov_model, "spherical")

  modelD <- ZIPLN(Abundance ~ 1, data = trichoptera, control = ZIPLN_param(covariance = "diagonal"))
  expect_equal(modelD$nb_param, p + p * 1 + 1)
  expect_equal(modelD$vcov_model, "diagonal")

  model <- ZIPLN(Abundance ~ 1, data = trichoptera, control = ZIPLN_param(covariance = "fixed", Omega = as.matrix(modelD$model_par$Omega)))
  expect_equal(model$nb_param, 0 + p * 1 + 1)
  expect_equal(model$model_par$Omega, modelD$model_par$Omega)
  expect_equal(model$vcov_model, "fixed")

})

test_that("ZIPLN fit: check sparse output and plot",  {

  myPLNfit <- ZIPLN(Abundance ~ 1, data = trichoptera, control = ZIPLN_param(trace = 1, penalty = 0.1))

  expect_is(myPLNfit, "ZIPLNfit_sparse")
  expect_is(myPLNfit, "ZIPLNfit")

  expect_equal(myPLNfit$vcov_model, "sparse")
  expect_true(igraph::is.igraph(myPLNfit$plot_network(output = "igraph", plot = FALSE)))
  expect_true(inherits(myPLNfit$plot_network(output = "corrplot", plot = FALSE), "Matrix"))

})

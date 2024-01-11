context("test-plnfit")

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

test_that("PLN fit: check classes, getters and field access",  {

  expect_output(model <- PLN(Abundance ~ 1, data = trichoptera,
                             control = PLN_param(trace = 1)),
"
 Initialization...
 Adjusting a full covariance PLN model with nlopt optimizer
 Post-treatments...
 DONE!"
  )

  expect_output(model <- PLN(Abundance ~ 1, data = trichoptera,
                             control = PLN_param(trace = 1, inception = model)),
"
 Initialization...
 Adjusting a full covariance PLN model with nlopt optimizer
 Post-treatments...
 DONE!"
  )


  expect_is(model, "PLNfit")

  expect_equal(model$n, nrow(trichoptera$Abundance))
  expect_equal(model$p, ncol(trichoptera$Abundance))
  expect_equal(model$d, 1)

  ## S3 methods: values
  expect_equal(coef(model), model$model_par$B)
  expect_equal(coef(model, type = "covariance"), sigma(model))
  expect_equal(sigma(model), model$model_par$Sigma)
  # expect_equal(vcov(model), model$vcov_coef)

  ## S3 methods: class
  expect_true(inherits(coef(model), "matrix"))
  expect_true(inherits(sigma(model), "matrix"))
  # expect_true(inherits(vcov(model), "dsCMatrix"))

  ## S3 methods: dimensions
  ## expect_equal(dim(vcov(model)), c(model$d * model$p, model$d * model$p))

  ## R6 bindings
  expect_is(model$latent, "matrix")
  expect_true(is.numeric(model$latent))
  expect_equal(dim(model$latent), c(model$n, model$p))

})

test_that("PLN fit: check print message",  {

  expect_output(model <- PLN(Abundance ~ 1, data = trichoptera))

  output <- paste(
"A multivariate Poisson Lognormal fit with full covariance model.
==================================================================",
capture_output(print(as.data.frame(round(model$criteria, digits = 3), row.names = ""))),
"==================================================================
* Useful fields
    $model_par, $latent, $latent_pos, $var_par, $optim_par
    $loglik, $BIC, $ICL, $loglik_vec, $nb_param, $criteria
* Useful S3 methods
    print(), coef(), sigma(), vcov(), fitted()
    predict(), predict_cond(), standard_error()",
    sep = "\n")

  expect_output(model$show(),
                output,
                fixed = TRUE)
  ## show and print are equivalent
  expect_equal(capture_output(model$show()),
               capture_output(model$print()))
})

test_that("PLN fit: Check prediction",  {

  model1     <- PLN(Abundance ~ 1, data = trichoptera, subset = 1:30)
  model1_off <- PLN(Abundance ~ 1 + offset(log(Offset)), data = trichoptera, subset = 1:30)
  model2     <- PLN(Abundance ~ Pressure, data = trichoptera, subset = 1:30)

  newdata <- trichoptera[31:49, ]
  # newdata$Abundance <- NULL

  pred1     <- predict(model1, newdata = newdata, type = "response")
  pred1_off <- predict(model1_off, newdata = newdata, type = "response")
  pred2     <- predict(model2, newdata = newdata, type = "response")
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

  ## R6 methods
  ## with offset, predictions should vary across samples
  expect_gte(min(apply(pred1_off, 2, sd)), .Machine$double.eps)
  newdata$Offset <- NULL
  ## without offsets, predictions should be the same for all samples
  expect_equal(unname(apply(pred1, 2, sd)), rep(0, ncol(pred1)))

  ## Unequal factor levels in train and prediction datasets
  suppressWarnings(
    toy_data <- prepare_data(
    counts     = matrix(c(1, 3, 1, 1), ncol = 1),
    covariates = data.frame(Cov = c("A", "B", "A", "A")),
    offset     = rep(1, 4))
  )
  model <- PLN(Abundance ~ Cov + offset(log(Offset)), data = toy_data[1:2,])
  expect_length(predict(model, newdata = toy_data[3:4, ], type = "r"), 2L)
})



test_that("PLN fit: Check conditional prediction",  {

  n_cond = 10
  p_cond = 2
  p <- ncol(trichoptera$Abundance)

  myPLN <- PLN(Abundance ~ Temperature, trichoptera)
  Yc <- trichoptera$Abundance[1:n_cond, 1:p_cond, drop=FALSE]

  newX <- data.frame(1, Temperature = trichoptera$Temperature[1:n_cond])

  pred <- predict_cond(myPLN, newX, Yc, type = "response")

  # check dimensions of the predictions (#TODO: modify pred$pred if we decide not to return M,S)
  expect_equal(dim(pred), c(n_cond,p-p_cond))

  # check if the RMSE of conditional predictions are greater than the marginal ones
  expect_gt(
    mean((trichoptera$Abundance[1:n_cond, (p_cond+1):p] -
            predict(myPLN, newdata = newX, type = "response")[1:n_cond, (p_cond+1):p])^2),
    mean((trichoptera$Abundance[1:n_cond, (p_cond+1):p] - pred)^2)
  )

  # check the dimension of the variational parameters when sent back
  pred <- predict_cond(myPLN, newX, Yc, type = "response", var_par = TRUE)
  expect_equal(dim(attr(pred, "M")), dim(pred))
  expect_equal(dim(attr(pred, "S")), c(p-p_cond, p-p_cond, n_cond))

})

test_that("PLN fit: Check number of parameters",  {

  p <- ncol(trichoptera$Abundance)

  model <- PLN(Abundance ~ 1, data = trichoptera)
  expect_equal(model$nb_param, p*(p+1)/2 + p * 1)

  model <- PLN(Abundance ~ 1 + Wind, data = trichoptera)
  expect_equal(model$nb_param, p*(p+1)/2 + p * 2)

  model <- PLN(Abundance ~ Group + 0 , data = trichoptera)
  expect_equal(model$nb_param, p*(p+1)/2 + p * nlevels(trichoptera$Group))

  modelS <- PLN(Abundance ~ 1, data = trichoptera, control = PLN_param(covariance = "spherical"))
  expect_equal(modelS$nb_param, 1 + p * 1)
  expect_equal(modelS$vcov_model, "spherical")

  modelD <- PLN(Abundance ~ 1, data = trichoptera, control = PLN_param(covariance = "diagonal"))
  expect_equal(modelD$nb_param, p + p * 1)
  expect_equal(modelD$vcov_model, "diagonal")

  model <- PLN(Abundance ~ 1, data = trichoptera, control = PLN_param(covariance = "fixed", Omega = as.matrix(modelD$model_par$Omega)))
  expect_equal(model$nb_param, 0 + p * 1)
  expect_equal(model$vcov_model, "fixed")

})

context("test-plnfit")

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

test_that("PLN fit: check classes, getters and field access",  {

  expect_output(model <- PLN(Abundance ~ 1, data = trichoptera,
                             control = list(trace = 2)),
"
 Initialization...
 Use LM after log transformation to define the inceptive model
 Adjusting a PLN model with full covariance model
 Post-treatments...
 DONE!"
  )

  expect_output(model <- PLN(Abundance ~ 1, data = trichoptera,
                             control = list(trace = 2, inception = model)),
"
 Initialization...
 User defined inceptive PLN model
 Adjusting a PLN model with full covariance model
 Post-treatments...
 DONE!"
  )


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
  expect_true(inherits(coef(model), "matrix"))
  expect_true(inherits(sigma(model), "matrix"))
  expect_true(class(vcov(model)) == "dgCMatrix")

  ## S3 methods: dimensions
  expect_equal(dim(vcov(model)), c(model$d * model$p, model$d * model$p))

  ## S3 methods errors
  expect_error(standard_error(model, type = "louis"),
               "Standard errors were not computed using the louis approximation. Try another approximation scheme.")

  ## R6 bindings
  expect_is(model$latent, "matrix")
  expect_true(is.numeric(model$latent))
  expect_equal(dim(model$latent), c(model$n, model$p))

  ## Fisher
  X <- model.matrix(Abundance ~ 1, data = trichoptera)
  fisher_louis <- model$compute_fisher(type = "louis", X = X)
  fisher_wald  <- model$compute_fisher(type = "wald", X = X)
  ## Louis fisher matrix is (component-wise) larger than its wald counterpart
  expect_gte(min(fisher_louis - fisher_wald), 0)

})

test_that("PLN fit: check print message",  {

  expect_output(model <- PLN(Abundance ~ 1, data = trichoptera))

  output <- paste(
"A multivariate Poisson Lognormal fit with full covariance model.
==================================================================",
capture_output(print(as.data.frame(round(model$criteria, digits = 3), row.names = ""))),
"==================================================================
* Useful fields
    $model_par, $latent, $var_par, $optim_par
    $loglik, $BIC, $ICL, $loglik_vec, $nb_param, $criteria
* Useful S3 methods
    print(), coef(), sigma(), vcov(), fitted(), predict(), standard_error()",
    sep = "\n")

  expect_output(model$show(),
                output,
                fixed = TRUE)
  ## show and print are equivalent
  expect_equal(capture_output(model$show()),
               capture_output(model$print()))
})

test_that("standard error fails for degenerate models", {
  trichoptera$X1 <- ifelse(trichoptera$Cloudiness <= 50, 0, 1)
  trichoptera$X2 <- 1 - trichoptera$X1
  # expect_warning(model <- PLN(Abundance ~ 1 + X1 + X2, data = trichoptera),
  #                "Something went wrong during model fitting!!\nMatrix A has missing values.")
  model <- PLN(Abundance ~ 1 + X1, data = trichoptera)
  ## Force a degenerate matrix in the FIM slot
  model$.__enclos_env__$private$FIM <- diag(0, nrow = model$p * model$d)
  expect_warning(std_err <- model$compute_standard_error(),
                 "Inversion of the Fisher information matrix failed with following error message:")
  expect_equal(std_err,
               matrix(NA, nrow = model$p, ncol = model$d,
                      dimnames = dimnames(coef(model))))
})

test_that("PLN fit: Check prediction",  {

  model1 <- PLN(Abundance ~ 1, data = trichoptera, subset = 1:30)
  model2 <- PLN(Abundance ~ Pressure, data = trichoptera, subset = 1:30)

  newdata <- trichoptera[31:49, ]
  newdata$Abundance <- NULL

  expect_gt(
    mean((trichoptera$Abundance[31:49, ] - predict(model1, newdata = newdata, type = "response"))^2),
    mean((trichoptera$Abundance[31:49, ] - predict(model2, newdata = newdata, type = "response"))^2)
  )

  ## R6 methods
  predictions <- model1$predict(newdata = newdata, type = "response")
  ## with offset, predictions should vary across samples
  expect_gte(min(apply(predictions, 2, sd)), 0)
  newdata$Offset <- NULL
  ## without offsets, predictions should be the same for all samples
  expect_equal(unname(apply(predictions, 2, sd)), rep(0, ncol(predictions)))

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

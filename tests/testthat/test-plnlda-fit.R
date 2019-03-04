context("test-plnfit")

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

test_that("Check classes, getters and field access",  {

  myPLNfit <- PLNLDA(Abundance ~ 1, data = trichoptera, grouping = Group)
  n <- nrow(trichoptera$Abundance)
  p <- ncol(trichoptera$Abundance)
  g <- length(unique(trichoptera$Group))

  expect_is(myPLNfit, "PLNLDAfit")
  expect_equal(myPLNfit$n, n)
  expect_equal(myPLNfit$p, p)
  expect_equal(myPLNfit$d, g)
  expect_equal(myPLNfit$rank, g - 1)

  expect_null(coef(myPLNfit))
  expect_lt(sum((myPLNfit$group_means - myPLNfit$myPLNfit_par$Theta)^2), .Machine$double.eps)
  expect_equal(sigma(myPLNfit), myPLNfit$model_par$Sigma)
  expect_equal(vcov(myPLNfit), myPLNfit$fisher$mat)

  expect_equal(class(myPLNfit$group_means), "matrix")
  expect_equal(class(sigma(myPLNfit)), "matrix")
  expect_true(class(vcov(myPLNfit)) == "dgCMatrix")

  expect_equal(dim(vcov(myPLNfit)), c(myPLNfit$d * myPLNfit$p, myPLNfit$d * myPLNfit$p))

  ## fields and active bindings
  expect_equal(dim(myPLNfit$model_par$B), c(p, p))
  expect_equal(dim(myPLNfit$model_par$Sigma), c(p, p))
  expect_equal(dim(myPLNfit$var_par$M), c(n, p))
  expect_equal(dim(myPLNfit$var_par$S), c(n, p))
  expect_equal(sum(myPLNfit$loglik_vec), myPLNfit$loglik)
  expect_lt(myPLNfit$BIC, myPLNfit$loglik)
  expect_lt(myPLNfit$ICL, myPLNfit$loglik)
  expect_lt(myPLNfit$ICL, myPLNfit$BIC)
  expect_gt(myPLNfit$R_squared, 0)
  expect_equal(myPLNfit$nb_param, p * (2 *g - 1))
  expect_equal(dim(myPLNfit$group_means), c(p, g))
  expect_equal(dim(myPLNfit$scores), c(n, myPLNfit$rank))
  expect_lt(sum(myPLNfit$percent_var) - 1, 1e-3)
  expect_equal(dim(myPLNfit$corr_map), c(p, myPLNfit$rank))

  ## S3 methods
  expect_equal(dim(fitted(myPLNfit)), c(n, p))
  expect_equal(sigma(myPLNfit), myPLNfit$model_par$Sigma)
  expect_equal(vcov(myPLNfit, "main"), myPLNfit$fisher$mat)
  expect_equal(vcov(myPLNfit, "covariance"), myPLNfit$model_par$Sigma)
  expect_equal(vcov(myPLNfit, "covariance"), sigma(myPLNfit))

  expect_true(inherits(plot(myPLNfit, map = "variable", plot = FALSE), "ggplot"))
  expect_true(inherits(plot(myPLNfit, map = "individual", plot = FALSE), "ggplot"))
  expect_true(inherits(plot(myPLNfit, map = "both", plot = FALSE), "grob"))

  ## R6 methods
  expect_true(inherits(myPLNfit$plot_correlation_map(plot = FALSE), "ggplot"))
  expect_true(inherits(myPLNfit$plot_individual_map(plot = FALSE), "ggplot"))
  expect_true(inherits(myPLNfit$plot_LDA(plot = FALSE), "grob"))

})

test_that("Check number of parameters",  {

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

## add tests for predictions, tests for fit
test_that("Predictions have the right dimensions.", {
  model1 <- PLNLDA(Abundance ~ 0 + offset(log(Offset)),
                   grouping = Group, data = trichoptera)

  expect_length(predict(model1, newdata = trichoptera, type = "response"),
                nrow(trichoptera))
  expect_is(predict(model1, newdata = trichoptera, type = "response"),
            "factor")
  expect_equal(dim(predict(model1, newdata = trichoptera)),
               c(nrow(trichoptera), length(levels(trichoptera$Group))))
  expect_equal(dim(predict(model1, newdata = trichoptera, type = "scores")),
               c(nrow(trichoptera), model1$rank))
  ## log-posterior probabilities are nonpositive
  expect_lt(max(predict(model1, newdata = trichoptera)), 0)
  ## Posterior probabilities are between 0 and 1
  expect_lte(max(predict(model1, newdata = trichoptera, scale = "prob")), 1)
  expect_gte(min(predict(model1, newdata = trichoptera, scale = "prob")), 0)

})



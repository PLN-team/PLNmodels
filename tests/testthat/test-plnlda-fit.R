context("test-plnldafit")

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

test_that("PLNLDA fit: check classes, getters and field access",  {

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
  expect_true(all(myPLNfit$percent_var >= 0))
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

test_that("PLNLDA fit: check print message",  {
  model <- PLNLDA(Abundance ~ 1, data = trichoptera, grouping = Group)
  expect_equal(capture_output(model$show()),
"Linear Discriminant Analysis for Poisson Lognormal distribution
==================================================================
 nb_param    loglik       BIC       ICL R_squared
      391 -895.3299 -1656.181 -1884.964 0.9549509
==================================================================
* Useful fields
    $model_par, $latent, $var_par, $optim_par
    $loglik, $BIC, $ICL, $loglik_vec, $nb_param, $criteria
* Useful S3 methods
    print(), coef(), vcov(), sigma(), fitted(), predict(), standard_error()
* Additional fields for LDA
    $percent_var, $corr_map, $scores, $group_means
* Additional S3 methods for LDA
    plot.PLNLDAfit(), predict.PLNLDAfit()"
)
})

test_that("plot_LDA works for binary groups:", {
  trichoptera$custom_group <- ifelse(trichoptera$Cloudiness <= 50, "Sunny", "Cloudy")
  model <- PLNLDA(Abundance ~ 1, data = trichoptera, grouping = custom_group)
  ## One axis only
  expect_true(inherits(model$plot_LDA(plot = FALSE), "grob"))
})

test_that("plot_LDA works for 4 or more axes:", {
  model <- PLNLDA(Abundance ~ 1, data = trichoptera, grouping = Group)
  expect_true(inherits(model$plot_LDA(nb_axes = 4, plot = FALSE), "grob"))
})

test_that("PLNLDA fit: Check number of parameters",  {

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

## add tests for predictions, tests for fit --------------------------------------------
test_that("Predictions have the right dimensions.", {
  model <- PLNLDA(Abundance ~ 0 + offset(log(Offset)),
                   grouping = Group, data = trichoptera)
  expect_length(predict(model, newdata = trichoptera, type = "response"),
                nrow(trichoptera))
  expect_is(predict(model, newdata = trichoptera, type = "response"),
            "factor")
  expect_equal(dim(predict(model, newdata = trichoptera)),
               c(nrow(trichoptera), length(levels(trichoptera$Group))))
  expect_equal(dim(predict(model, newdata = trichoptera, type = "scores")),
               c(nrow(trichoptera), model$rank))
  ## log-posterior probabilities are nonpositive
  expect_lt(max(predict(model, newdata = trichoptera)), 0)
  ## Posterior probabilities are between 0 and 1
  expect_lte(max(predict(model, newdata = trichoptera, scale = "prob")), 1)
  expect_gte(min(predict(model, newdata = trichoptera, scale = "prob")), 0)
})

test_that("Predictions are not affected by inclusion of an intercept.", {
  model0 <- PLNLDA(Abundance ~ 0 + offset(log(Offset)),
                  grouping = Group, data = trichoptera)
  model1 <- PLNLDA(Abundance ~ 1 + offset(log(Offset)),
                   grouping = Group, data = trichoptera)
  expect_equal(predict(model0, newdata = trichoptera),
               predict(model1, newdata = trichoptera))
})

test_that("Prediction fails for non positive prior probabilities.", {
  model <- PLNLDA(Abundance ~ 0 + offset(log(Offset)),
                  grouping = Group, data = trichoptera)
  nb_groups <- length(levels(trichoptera$Group))
  expect_error(predict(model, newdata = trichoptera, type = "response", prior = rep(0, nb_groups)),
               "Prior group proportions should be positive.")
  expect_error(predict(model, newdata = trichoptera, type = "response", prior = rep(NA, nb_groups)),
               "Prior group proportions should be positive.")
})

test_that("Predictions succeeds for positive prior probabilities.", {
  model <- PLNLDA(Abundance ~ 0 + offset(log(Offset)),
                  grouping = Group, data = trichoptera)
  nb_groups <- model$rank + 1
  expect_length(predict(model, newdata = trichoptera, type = "response", prior = rep(1, nb_groups)),
                nrow(trichoptera))
  ## Predictions can be set to 6 (sixth class) by putting enough prior weight on it
  prior <- rep(1, nb_groups)
  ## the posterior difference between class 6 and best class in log scale is at most 248.5,
  ## which corresponds to 1e108
  prior[6] <- 1e108
  expect_equal(predict(model, newdata = trichoptera, type = "response", prior = prior),
               setNames(factor(rep(6, nrow(trichoptera)), levels = levels(trichoptera$Group)),
                        rownames(trichoptera)))
})


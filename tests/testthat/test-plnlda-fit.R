context("test-plnldafit")

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
model0 <- PLNLDA(Abundance ~ 0 + offset(log(Offset)),
                 grouping = Group, data = trichoptera)
model  <- PLNLDA(Abundance ~ 1 + offset(log(Offset)),
                 grouping = Group, data = trichoptera)

test_that("PLNLDA fit: check classes, getters and field access",  {

  n <- nrow(trichoptera$Abundance)
  p <- ncol(trichoptera$Abundance)
  g <- length(unique(trichoptera$Group))

  expect_is(model, "PLNLDAfit")
  expect_equal(model$n, n)
  expect_equal(model$p, p)
  expect_equal(model$d, g)
  expect_equal(model$rank, g - 1)

  expect_null(coef(model))
  expect_lt(sum((model$group_means - model$model_par$Theta)^2), .Machine$double.eps)
  expect_equal(sigma(model), model$model_par$Sigma)
  expect_equal(vcov(model), model$fisher$mat)

  expect_true(inherits(model$group_means, "matrix"))
  expect_true(inherits(sigma(model), "matrix"))
  expect_true(class(vcov(model)) == "dgCMatrix")

  expect_equal(dim(vcov(model)), c(model$d * model$p, model$d * model$p))

  ## fields and active bindings
  expect_equal(dim(model$model_par$B), c(p, p))
  expect_equal(dim(model$model_par$Sigma), c(p, p))
  expect_equal(dim(model$var_par$M), c(n, p))
  expect_equal(dim(model$var_par$S), c(n, p))
  expect_equal(sum(model$loglik_vec), model$loglik)
  expect_lt(model$BIC, model$loglik)
  expect_lt(model$ICL, model$loglik)
  expect_gt(model$R_squared, 0)
  expect_equal(model$nb_param, p * (2 *g - 1))
  expect_equal(dim(model$group_means), c(p, g))
  expect_equal(dim(model$scores), c(n, model$rank))
  expect_true(all(model$percent_var >= 0))
  expect_equal(dim(model$corr_map), c(p, model$rank))

  ## S3 methods
  expect_equal(dim(fitted(model)), c(n, p))
  expect_equal(sigma(model), model$model_par$Sigma)
  expect_equal(vcov(model, "main"), model$fisher$mat)
  expect_equal(vcov(model, "covariance"), model$model_par$Sigma)
  expect_equal(vcov(model, "covariance"), sigma(model))

  expect_true(inherits(plot(model, map = "variable", plot = FALSE), "ggplot"))
  expect_true(inherits(plot(model, map = "individual", plot = FALSE), "ggplot"))
  expect_true(inherits(plot(model, map = "both", plot = FALSE), "grob"))

  ## R6 methods
  expect_true(inherits(model$plot_correlation_map(plot = FALSE), "ggplot"))
  expect_true(inherits(model$plot_individual_map(plot = FALSE), "ggplot"))
  expect_true(inherits(model$plot_LDA(plot = FALSE), "grob"))
})

test_that("PLNLDA fit: check print message",  {

  output <- paste(
"Linear Discriminant Analysis for Poisson Lognormal distribution
==================================================================",
capture_output(print(as.data.frame(round(model$criteria, digits = 3), row.names = ""))),
"==================================================================
* Useful fields
    $model_par, $latent, $var_par, $optim_par
    $loglik, $BIC, $ICL, $loglik_vec, $nb_param, $criteria
* Useful S3 methods
    print(), coef(), sigma(), vcov(), fitted(), predict(), standard_error()
* Additional fields for LDA
    $percent_var, $corr_map, $scores, $group_means
* Additional S3 methods for LDA
    plot.PLNLDAfit(), predict.PLNLDAfit()",
sep = "\n"
)

  expect_output(model$show(),
                output,
                fixed = TRUE)
})

test_that("plot_LDA works for binary groups:", {
  trichoptera$custom_group <- ifelse(trichoptera$Cloudiness <= 50, "Sunny", "Cloudy")

  ## One axis only
  expect_true(inherits(model$plot_LDA(plot = FALSE), "grob"))
})

test_that("plot_LDA works for 4 or more axes:", {
  expect_true(inherits(model$plot_LDA(nb_axes = 4, plot = FALSE), "grob"))
})

# test_that("PLNLDA fit: Check number of parameters",  {
#
#   p <- ncol(trichoptera$Abundance)
#
#   mdl <- PLN(Abundance ~ 1, data = trichoptera)
#   expect_equal(mdl$nb_param, p*(p+1)/2 + p * 1)
#
#   mdl <- PLN(Abundance ~ 1 + Wind, data = trichoptera)
#   expect_equal(mdl$nb_param, p*(p+1)/2 + p * 2)
#
#   mdl <- PLN(Abundance ~ Group + 0 , data = trichoptera)
#   expect_equal(mdl$nb_param, p*(p+1)/2 + p * nlevels(trichoptera$Group))
#
#   mdl <- PLN(Abundance ~ 1, data = trichoptera, control = list(covariance = "diagonal"))
#   expect_equal(mdl$nb_param, p + p * 1)
#
#   mdl <- PLN(Abundance ~ 1, data = trichoptera, control = list(covariance = "spherical"))
#   expect_equal(mdl$nb_param, 1 + p * 1)
#
# })

## add tests for predictions, tests for fit --------------------------------------------
test_that("Predictions have the right dimensions.", {
  predictions_response <- predict(model, newdata = trichoptera, type = "response")
  predictions_post <- predict(model, newdata = trichoptera)
  predictions_prob <- predict(model, newdata = trichoptera, scale = "prob")
  predictions_score <- predict(model, newdata = trichoptera)
  ## Train = Test
  expect_length(predictions_response, nrow(trichoptera))
  expect_is(predictions_response, "factor")
  expect_equal(dim(predictions_post),
               c(nrow(trichoptera), length(levels(trichoptera$Group))))
  expect_equal(dim(predictions_score),
               c(nrow(trichoptera), length(levels(trichoptera$Group))))
  ## log-posterior probabilities are nonpositive
  expect_lt(max(predictions_post), 0)
  ## Posterior probabilities are between 0 and 1
  expect_lte(max(predictions_prob), 1)
  expect_gte(min(predictions_prob), 0)
  ## Train != Test

  ## test failing due to core dump
  # test <- 1:nrow(trichoptera) < (nrow(trichoptera)/2)
  # expect_equal(dim(predict(model, newdata = trichoptera[test, ], type = "scores")),
  #              c(sum(test), model$rank))


})

test_that("Predictions are not affected by inclusion of an intercept.", {
  expect_equal(predict(model0, newdata = trichoptera),
               predict(model , newdata = trichoptera))
})

# test_that("Predictions work when train and test data have different factor levels.", {
#   suppressWarnings(
#     toy_data <- prepare_data(
#       counts     = matrix(c(1, 4, 2, 1,
#                             1, 8, 2, 1),
#                           ncol = 2),
#       covariates = data.frame(Cov   = c("A", "B", "B", "A"),
#                               Group = c("a", "b", "a", "a")),
#       offset     = matrix(rep(1, 8), ncol = 2)
#   ))
#   suppressWarnings(
#     model <- PLNLDA(Abundance ~ Cov + offset(log(Offset)),
#                     grouping = Group,
#                     data     = toy_data[1:3,])
#   )
#   # expect_length(predict(model, newdata = toy_data[c(1,4), ], type = "r"),
#   #               2L)
#   # expect_identical(predict(model, newdata = toy_data[c(1,4), ], type = "r"),
#   #                  factor(c(`1` = "a", `4` = "a"), levels = c("a", "b")))
# })


test_that("Prediction fails for non positive prior probabilities.", {
  nb_groups <- length(levels(trichoptera$Group))
  expect_error(predict(model, newdata = trichoptera, type = "response", prior = rep(0, nb_groups)),
               "Prior group proportions should be positive.")
  expect_error(predict(model, newdata = trichoptera, type = "response", prior = rep(NA, nb_groups)),
               "Prior group proportions should be positive.")
})

# Commenting due to failure on CRAN
# test_that("Predictions succeeds for positive prior probabilities.", {
#   nb_groups <- model$rank + 1
#   expect_length(predict(model, newdata = trichoptera, type = "response", prior = rep(1, nb_groups)),
#                 nrow(trichoptera))
#   ## Predictions can be set to 6 (sixth class) by putting enough prior weight on it
#   prior <- rep(1, nb_groups)
#   ## the posterior difference between class 6 and best class in log scale is at most 248.5,
#   ## which corresponds to 1e108
#   prior[6] <- 1e108
#   expect_equal(predict(model, newdata = trichoptera, type = "response", prior = prior),
#                setNames(factor(rep(6, nrow(trichoptera)), levels = levels(trichoptera$Group)),
#                         rownames(trichoptera)))
# })
#

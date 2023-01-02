context("test-plnpcafit")

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
models <- PLNPCA(Abundance ~ 1, data = trichoptera)
X <- model.matrix(Abundance ~ 1, data = trichoptera)
myPLNfit <- getModel(models, 3)

test_that("PLNPCA fit: check classes, getters and field access", {

  Y <- as.matrix(trichoptera$Abundance)
  n <- nrow(Y); p <- ncol(Y)
  O <- matrix(0, nrow = n, ncol = p)
  w <- rep(1, n)

  ## fields and active bindings
  expect_equal(dim(myPLNfit$latent), dim(Y))
  expect_equal(dim(myPLNfit$model_par$Theta), c(ncol(Y), ncol(X)))
  expect_equal(dim(myPLNfit$model_par$C), c(ncol(Y), myPLNfit$rank))
  expect_equal(dim(myPLNfit$model_par$Sigma), c(ncol(Y), ncol(Y)))
  expect_equal(dim(myPLNfit$var_par$M), c(nrow(Y), myPLNfit$rank))
  expect_equal(dim(myPLNfit$var_par$S), c(nrow(Y), myPLNfit$rank))
  expect_equal(sum(myPLNfit$loglik_vec), myPLNfit$loglik)
  expect_lt(myPLNfit$BIC, myPLNfit$loglik)
  expect_lt(myPLNfit$ICL, myPLNfit$loglik)
#  expect_lt(myPLNfit$ICL, myPLNfit$BIC) ## entropy could be positive
  expect_gt(myPLNfit$R_squared, 0)
  expect_equal(myPLNfit$nb_param, p + p * myPLNfit$rank - myPLNfit$rank * (myPLNfit$rank - 1) / 2 )
  expect_equal(dim(myPLNfit$rotation), c(p, myPLNfit$rank))
  expect_equal(dim(myPLNfit$scores), c(n, myPLNfit$rank))
  expect_true(all(myPLNfit$percent_var >= 0))
  expect_equal(dim(myPLNfit$corr_circle), c(p, myPLNfit$rank))
  ## Eigenvalues, informations about individuals and variables
  expect_equal(dim(myPLNfit$eig), c(myPLNfit$rank, 3))
  ## $var
  expect_equal(dim(myPLNfit$var$coord), c(p, myPLNfit$rank))
  expect_equal(dim(myPLNfit$var$cor), c(p, myPLNfit$rank))
  expect_equal(dim(myPLNfit$var$cos2), c(p, myPLNfit$rank))
  expect_equal(dim(myPLNfit$var$contrib), c(p, myPLNfit$rank))
  ## $ind
  expect_equal(dim(myPLNfit$ind$coord), c(n, myPLNfit$rank))
  expect_equal(dim(myPLNfit$ind$cos2), c(n, myPLNfit$rank))
  expect_equal(dim(myPLNfit$ind$contrib), c(n, myPLNfit$rank))
  expect_equal(length(myPLNfit$ind$dist), n)

  ## S3 methods
  expect_equal(coefficients(myPLNfit), myPLNfit$model_par$B)
  expect_equal(dim(fitted(myPLNfit)), dim(Y))
  expect_equal(sigma(myPLNfit), myPLNfit$model_par$Sigma)
  expect_equal(vcov(myPLNfit, "main"), myPLNfit$vcov_coef)
  expect_equal(vcov(myPLNfit, "covariance"), myPLNfit$model_par$Sigma)
  expect_equal(vcov(myPLNfit, "covariance"), sigma(myPLNfit))
  expect_warning(sd_plnfit <- standard_error(myPLNfit), "Standard error of B is not implemented yet for PLNPCA models")
  expect_equal(dim(sd_plnfit), dim(coefficients(myPLNfit)))
  expect_error(standard_error(myPLNfit, parameter = "Omega"), "Omega is not estimated as such in PLNPCA models")

  expect_true(inherits(plot(myPLNfit, map = "variable", plot = FALSE), "ggplot"))
  expect_true(inherits(plot(myPLNfit, map = "individual", plot = FALSE), "ggplot"))
  expect_true(inherits(plot(myPLNfit, map = "both", plot = FALSE), "grob"))

  ## R6 methods: Plots
  expect_true(inherits(myPLNfit$plot_correlation_circle(plot = FALSE), "ggplot"))
  expect_true(inherits(myPLNfit$plot_individual_map(plot = FALSE), "ggplot"))
  expect_true(inherits(myPLNfit$plot_PCA(plot = FALSE), "grob"))

  ## R6 methods: VEstep
  ve_results <- myPLNfit$optimize_vestep(covariates = X, offsets = O, responses = Y)
  expect_equal(dim(ve_results$M), c(n, myPLNfit$rank))
  expect_equal(dim(ve_results$S), c(n, myPLNfit$rank))
  expect_length(ve_results$Ji, n)
  expect_equal(ve_results$M, unname(myPLNfit$var_par$M), tolerance = 1e-1)
  ## R6 methods: project()
  scores <- myPLNfit$project(newdata = trichoptera)
  expect_equal(dim(scores), dim(myPLNfit$scores))
  expect_equal(dimnames(scores), dimnames(myPLNfit$scores))
  expect_equal(scores, myPLNfit$scores, tolerance = 1e-1)

  ## Class
  expect_true(inherits(myPLNfit, "PCA"))
})

test_that("Bindings for factoextra return sensible values", {
  ## $eig
  expect_gte(min(myPLNfit$eig[, "eigenvalue"]), 0)
  expect_gte(min(myPLNfit$eig[, "percentage of variance"]), 0)
  expect_lte(max(myPLNfit$eig[, "percentage of variance"]), 100 * myPLNfit$R_squared)
  expect_equal(tail(myPLNfit$eig[, "cumulative percentage of variance"], n = 1), 100 * myPLNfit$R_squared, check.attributes = FALSE)
  ## $var
  .var <- myPLNfit$var
  cor_range <- range(.var$cor)
  expect_gte(cor_range[1], -1)
  expect_lte(cor_range[2], 1)
  cos2_range <- range(.var$cos2)
  expect_gte(cos2_range[1], 0)
  expect_lte(cos2_range[2], 1)
  expect_equal(rowSums(.var$cos2), rep(1, myPLNfit$p), check.attributes = FALSE)
  expect_equal(colSums(.var$contrib), rep(100, myPLNfit$rank), check.attributes = FALSE)
  ## $ind
  .ind <- myPLNfit$ind
  cos2_range <- range(.ind$cos2)
  expect_gte(cos2_range[1], 0)
  expect_lte(cos2_range[2], 1)
  expect_equal(rowSums(.ind$cos2), rep(1, myPLNfit$n), check.attributes = FALSE)
  expect_equal(colSums(.ind$contrib), rep(100, myPLNfit$rank), check.attributes = FALSE)
  expect_equal(colSums(.ind$coord), rep(0, myPLNfit$rank), check.attributes = FALSE)
})


test_that("Louis-type Fisher matrices are not available for PLNPCA", {
  expect_error(myPLNfit$compute_fisher(type = "louis", X = X))
})

test_that("plot_PCA works one axis only:", {
  model1 <- getModel(models, 1)
  ## One axis only
  expect_true(inherits(model1$plot_PCA(plot = FALSE), "grob"))
})

test_that("plot_PCA works for 4 or more axes:", {
  model4 <- getModel(models, 4)
  expect_true(inherits(model4$plot_PCA(nb_axes = 4, plot = FALSE), "grob"))
})

test_that("PLNPCA fit: check print message",  {
  output <- paste(
"Poisson Lognormal with rank constrained for PCA (rank = 3)
==================================================================",
capture_output(print(as.data.frame(round(myPLNfit$criteria, digits = 3), row.names = ""))),
"==================================================================
* Useful fields
    $model_par, $latent, $latent_pos, $var_par, $optim_par
    $loglik, $BIC, $ICL, $loglik_vec, $nb_param, $criteria
* Useful S3 methods
    print(), coef(), sigma(), vcov(), fitted()
    predict(), predict_cond(), standard_error()
* Additional fields for PCA
    $percent_var, $corr_circle, $scores, $rotation, $eig, $var, $ind
* Additional S3 methods for PCA
    plot.PLNPCAfit()",
sep = "\n")

  expect_output(myPLNfit$show(),
                output,
                fixed = TRUE)
})

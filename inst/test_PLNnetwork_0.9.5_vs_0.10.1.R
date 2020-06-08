## Need Shift + Ctrl + F10 first
rm(list=ls())
## get oaks data set
load("inst/case_studies/oaks_mildew/oaks_alphitoides.RData")
class(oaks$Abundancies) <- class(oaks$Abundancies)[-match("AsIs", class(oaks$Abundancies))]

model <- model.frame(Abundancies ~ 0 + treeStatus + offset(log(sequencingEffort)), data = oaks)
Y <- model.response(model)
X <- model.matrix(terms(model), model)

loglik_param <- function(myPLN_fit) {
  Theta <- myPLN_fit$model_par$Theta
  Omega <- solve(myPLN_fit$model_par$Sigma)
  mu <- X %*% t(Theta)
  M <- myPLN_fit$var_par$M
  S <- myPLN_fit$var_par$S
  if (ncol(S) == 1) S <- matrix(c(S), nrow(Y), ncol(Y))
  Z <- myPLN_fit$latent
  A <- myPLN_fit$fitted

  Ki <- - rowSums(lfactorial(Y)) + .5 * (1+(1-ncol(Y))* log(2*pi))
  loglik <- rowSums(Y * Z - A + .5*log(S) - .5*( (M %*% Omega) * M + S %*% diag(diag(Omega))) ) + .5 * determinant(Omega,logarithm = TRUE)$modulus + Ki
  loglik
}

library(PLNmodels, lib.loc = "~/R/x86_64-pc-linux-gnu-library/4.0/")
cat("\nVersion", as.character(packageVersion("PLNmodels")), "\n")

## simple PLN
time_dev <- system.time(myPLN_dev <- PLNnetwork(Abundancies ~ 0 + treeStatus + offset(log(sequencingEffort)), data = oaks, control_main = list(trace = 2)))

detach("package:PLNmodels")

library(PLNmodels, lib.loc = "~/R/x86_64-pc-linux-gnu-library/PLNmodels_9.5/")
cat("\nVersion", as.character(packageVersion("PLNmodels")), "\n")

## simple PLN
time_CRAN <- system.time(myPLN_CRAN <- PLNnetwork(Abundancies ~ 0 + treeStatus + offset(log(sequencingEffort)), data = oaks, control_main = list(trace = 2)))


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
time_full      <- system.time(myPLN <- PLN(Abundancies ~ 0 + treeStatus + offset(log(sequencingEffort)), data = oaks, control = list(trace = 0)))
time_diagonal  <- system.time(myPLN_diagonal <- PLN(Abundancies ~ 0 + treeStatus + offset(log(sequencingEffort)), data = oaks, control = list(trace = 0, covariance = "diagonal")))
time_spherical <- system.time(myPLN_spherical <- PLN(Abundancies ~ 0 + treeStatus + offset(log(sequencingEffort)), data = oaks, control = list(trace = 0, covariance = "spherical")))

data.frame(nb_param = c(myPLN$nb_param, myPLN_diagonal$nb_param, myPLN_spherical$nb_param),
            time = c(time_full[3], time_diagonal[3], time_spherical[3]),
            loglik = c(sum(loglik_param(myPLN)), sum(loglik_param(myPLN_diagonal)), sum(loglik_param(myPLN_spherical))),
           row.names = c("full", "diagonal", "spherical")) %>%
  knitr::kable() %>% print()


detach("package:PLNmodels")

library(PLNmodels, lib.loc = "~/R/x86_64-pc-linux-gnu-library/PLNmodels_9.5/")
cat("\nVersion", as.character(packageVersion("PLNmodels")), "\n")

## simple PLN
time_full      <- system.time(myPLN <- PLN(Abundancies ~ 0 + treeStatus + offset(log(sequencingEffort)), data = oaks, control = list(trace = 0)))
time_diagonal  <- system.time(myPLN_diagonal <- PLN(Abundancies ~ 0 + treeStatus + offset(log(sequencingEffort)), data = oaks, control = list(trace = 0, covariance = "diagonal")))
time_spherical <- system.time(myPLN_spherical <- PLN(Abundancies ~ 0 + treeStatus + offset(log(sequencingEffort)), data = oaks, control = list(trace = 0, covariance = "spherical")))

data.frame(nb_param = c(myPLN$nb_param, myPLN_diagonal$nb_param, myPLN_spherical$nb_param),
            time = c(time_full[3], time_diagonal[3], time_spherical[3]),
            loglik = c(sum(loglik_param(myPLN)), sum(loglik_param(myPLN_diagonal)), sum(loglik_param(myPLN_spherical))),
           row.names = c("full", "diagonal", "spherical")) %>%
  knitr::kable() %>% print()

detach("package:PLNmodels", unload = TRUE)

#  loglik = sum(data.Y % Z - A + .5*log(S) - .5*( ((M - mu) * Omega) % (M - mu) + S * diagmat(Omega)), 1) + .5 * real(log_det(Omega)) + data.Ki ;

# loglik = sum(data.Y % Z - A + .5*log(S) - .5*( (M * Omega) % M + S * diagmat(Omega)), 1) + .5 * real(log_det(Omega))




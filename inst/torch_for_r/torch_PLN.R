data("oaks")
system.time(myPLN_nlopt <-
              PLN(Abundance ~ 1  + offset(log(Offset)),
                  data = oaks, control = PLN_param(backend = "nlopt"))
)
system.time(myPLN_torch <-
              PLN(Abundance ~ 1  + offset(log(Offset)),
                  data = oaks, control = PLN_param(backend = "torch",
                                                   inception = myPLN_nlopt))
)

plot(-myPLN_torch$optim_par$objective[-(1:10)], type = "l", log= 'y')
x11()
par(mfrow = c(2,2))
plot(myPLN_torch$model_par$Theta,
     myPLN_nlopt$model_par$Theta); abline(0, 1)
plot(myPLN_torch$model_par$Sigma,
     myPLN_nlopt$model_par$Sigma); abline(0, 1)
plot(myPLN_torch$var_par$M,
     myPLN_nlopt$var_par$M); abline(0, 1)
plot(myPLN_torch$var_par$S2,
     myPLN_nlopt$var_par$S2); abline(0, 1)
par()
myPLN_torch$loglik
myPLN_nlopt$loglik

system.time(myPLN_torch <-
              PLN(Abundance ~ 1  + offset(log(Offset)),
                  data = oaks, control = PLN_param(backend = "torch", covariance = "spherical"))
)
system.time(myPLN_nlopt <-
              PLN(Abundance ~ 1  + offset(log(Offset)),
                  data = oaks, control = PLN_param(backend = "nlopt", covariance = "spherical"))
)
plot(myPLN_torch$model_par$Theta,
     myPLN_nlopt$model_par$Theta); abline(0, 1)
plot(myPLN_torch$model_par$Sigma,
     myPLN_nlopt$model_par$Sigma); abline(0, 1)
myPLN_torch$loglik
myPLN_nlopt$loglik

system.time(myPLN_torch <-
              PLN(Abundance ~ 1  + offset(log(Offset)),
                  data = oaks, control = PLN_param(backend = "torch", covariance = "diagonal"))
)
system.time(myPLN_nlopt <-
              PLN(Abundance ~ 1  + offset(log(Offset)),
                  data = oaks, control = PLN_param(backend = "nlopt", covariance = "diagonal"))
)
plot(myPLN_torch$model_par$Theta,
     myPLN_nlopt$model_par$Theta); abline(0, 1)
plot(myPLN_torch$model_par$Sigma,
     myPLN_nlopt$model_par$Sigma); abline(0, 1)
myPLN_torch$loglik
myPLN_nlopt$loglik

data("trichoptera")
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
system.time(myPLN_torch <-
              PLN(Abundance ~ 1  + offset(log(Offset)),
                  data = trichoptera, control = PLN_param(backend = "torch"))
)
system.time(myPLN_nlopt <-
              PLN(Abundance ~ 1  + offset(log(Offset)),
                  data = trichoptera, control = PLN_param(backend = "nlopt"))
)

x11()
par(mfrow = c(2,2))
plot(myPLN_torch$model_par$Theta,
     myPLN_nlopt$model_par$Theta); abline(0, 1)
plot(myPLN_torch$model_par$Sigma,
     myPLN_nlopt$model_par$Sigma); abline(0, 1)
plot(myPLN_torch$var_par$M,
     myPLN_nlopt$var_par$M); abline(0, 1)
plot(myPLN_torch$var_par$S2,
     myPLN_nlopt$var_par$S2); abline(0, 1)
myPLN_torch$loglik
myPLN_nlopt$loglik

system.time(myPLN_torch <-
              PLN(Abundance ~ 1  + offset(log(Offset)),
                  data = trichoptera, control = PLN_param(backend = "torch", covariance = "spherical"))
)
system.time(myPLN_nlopt <-
              PLN(Abundance ~ 1  + offset(log(Offset)),
                  data = trichoptera, control = PLN_param(backend = "nlopt", covariance = "spherical"))
)
plot(myPLN_torch$model_par$Theta,
     myPLN_nlopt$model_par$Theta); abline(0, 1)
plot(myPLN_torch$model_par$Sigma,
     myPLN_nlopt$model_par$Sigma); abline(0, 1)
myPLN_torch$loglik
myPLN_nlopt$loglik

system.time(myPLN_torch <-
              PLN(Abundance ~ 1  + offset(log(Offset)),
                  data = trichoptera, control = PLN_param(backend = "torch", covariance = "diagonal"))
)
system.time(myPLN_nlopt <-
              PLN(Abundance ~ 1  + offset(log(Offset)),
                  data = trichoptera, control = PLN_param(backend = "nlopt", covariance = "diagonal"))
)

par(mfrow = c(2,2))
plot(myPLN_torch$model_par$Theta,
     myPLN_nlopt$model_par$Theta); abline(0, 1)
plot(myPLN_torch$model_par$Sigma,
     myPLN_nlopt$model_par$Sigma); abline(0, 1)
plot(myPLN_torch$var_par$M,
     myPLN_nlopt$var_par$M); abline(0, 1)
plot(myPLN_torch$var_par$S2,
     myPLN_nlopt$var_par$S2); abline(0, 1)
myPLN_torch$loglik
myPLN_nlopt$loglik

data("mollusk")
mollusk <- prepare_data(mollusk$Abundance, mollusk$Covariate)
system.time(myPLN_torch <-
              PLN(Abundance ~ 1  + offset(log(duration)),
                  data = mollusk, control = PLN_param(backend = "torch"))
)
system.time(myPLN_nlopt <-
              PLN(Abundance ~ 1  + offset(log(duration)),
                  data = mollusk, control = PLN_param(backend = "nlopt"))
)

par(mfrow = c(2,2))
plot(myPLN_torch$model_par$Theta,
     myPLN_nlopt$model_par$Theta); abline(0, 1)
plot(myPLN_torch$model_par$Sigma,
     myPLN_nlopt$model_par$Sigma); abline(0, 1)
plot(myPLN_torch$var_par$M,
     myPLN_nlopt$var_par$M); abline(0, 1)
plot(myPLN_torch$var_par$S2,
     myPLN_nlopt$var_par$S2); abline(0, 1)
myPLN_torch$loglik
myPLN_nlopt$loglik

data("barents")
system.time(myPLN_torch <-
              PLN(Abundance ~ 1 + Depth + Temperature + offset(log(Offset)),
                  data = barents, control = PLN_param(backend = "torch"))
)
system.time(myPLN_nlopt <-
              PLN(Abundance ~ 1 + Depth + Temperature + offset(log(Offset)),
                  data = barents, control = PLN_param(backend = "nlopt",
                                                      config_optim = list(xtol_rel=1e-9)))
)

par(mfrow = c(2,2))
plot(myPLN_torch$model_par$Theta,
     myPLN_nlopt$model_par$Theta); abline(0, 1)
plot(myPLN_torch$model_par$Sigma,
     myPLN_nlopt$model_par$Sigma); abline(0, 1)
plot(myPLN_torch$var_par$M,
     myPLN_nlopt$var_par$M); abline(0, 1)
plot(myPLN_torch$var_par$S2,
     myPLN_nlopt$var_par$S2); abline(0, 1)
myPLN_torch$loglik
myPLN_nlopt$loglik


params <- PLNmodels:::create_parameters()
Theta <- params$Theta
## Extract X
X <- params$X
## Extract Y
Y <- rPLN(n = nrow(X), mu = tcrossprod(X, Theta), Sigma = params$Sigma, depths = params$depths)
data <- prepare_data(Y, X, offset = "none")
O <- rowSums(Y)
myPLN_nlopt <- PLN(Abundance ~ 0 + . + offset(log(O)), data = data,
             control = PLN_param(backend = "nlopt", covariance = "fixed", Omega = solve(params$Sigma)))
myPLN_torch <- PLN(Abundance ~ 0 + . + offset(log(O)), data = data,
             control = PLN_param(backend = "torch", covariance = "fixed", Omega = solve(params$Sigma)))
par(mfrow = c(2,2))
plot(myPLN_torch$model_par$Theta,
     myPLN_nlopt$model_par$Theta); abline(0, 1)
plot(myPLN_torch$model_par$Sigma,
     myPLN_nlopt$model_par$Sigma); abline(0, 1)
plot(myPLN_torch$var_par$M,
     myPLN_nlopt$var_par$M); abline(0, 1)
plot(myPLN_torch$var_par$S2,
     myPLN_nlopt$var_par$S2); abline(0, 1)
myPLN_torch$loglik
myPLN_nlopt$loglik

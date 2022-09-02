data("oaks")

system.time(myPLN_torch <-
              PLN(Abundance ~ 1  + offset(log(Offset)),
                  data = oaks, control = list(backend = "torch"))
)
system.time(myPLN_nlopt <-
              PLN(Abundance ~ 1  + offset(log(Offset)),
                  data = oaks, control = list(backend = "nlopt"))
)

plot(myPLN_torch$model_par$Theta,
     myPLN_nlopt$model_par$Theta); abline(0, 1)
plot(myPLN_torch$model_par$Sigma,
     myPLN_nlopt$model_par$Sigma); abline(0, 1)
myPLN_torch$loglik
myPLN_nlopt$loglik

system.time(myPLN_torch <-
              PLN(Abundance ~ 1  + offset(log(Offset)),
                  data = oaks, control = list(backend = "torch", covariance = "spherical"))
)
system.time(myPLN_nlopt <-
              PLN(Abundance ~ 1  + offset(log(Offset)),
                  data = oaks, control = list(backend = "nlopt", covariance = "spherical"))
)
plot(myPLN_torch$model_par$Theta,
     myPLN_nlopt$model_par$Theta); abline(0, 1)
plot(myPLN_torch$model_par$Sigma,
     myPLN_nlopt$model_par$Sigma); abline(0, 1)
myPLN_torch$loglik
myPLN_nlopt$loglik

system.time(myPLN_torch <-
              PLN(Abundance ~ 1  + offset(log(Offset)),
                  data = oaks, control = list(backend = "torch", covariance = "diagonal"))
)
system.time(myPLN_nlopt <-
              PLN(Abundance ~ 1  + offset(log(Offset)),
                  data = oaks, control = list(backend = "nlopt", covariance = "diagonal"))
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
                  data = trichoptera, control = list(backend = "torch"))
)
system.time(myPLN_nlopt <-
              PLN(Abundance ~ 1  + offset(log(Offset)),
                  data = trichoptera, control = list(backend = "nlopt"))
)

plot(myPLN_torch$model_par$Theta,
     myPLN_nlopt$model_par$Theta); abline(0, 1)
plot(myPLN_torch$model_par$Sigma,
     myPLN_nlopt$model_par$Sigma); abline(0, 1)
myPLN_torch$loglik
myPLN_nlopt$loglik

system.time(myPLN_torch <-
              PLN(Abundance ~ 1  + offset(log(Offset)),
                  data = trichoptera, control = list(backend = "torch", covariance = "spherical"))
)
system.time(myPLN_nlopt <-
              PLN(Abundance ~ 1  + offset(log(Offset)),
                  data = trichoptera, control = list(backend = "nlopt", covariance = "spherical"))
)
plot(myPLN_torch$model_par$Theta,
     myPLN_nlopt$model_par$Theta); abline(0, 1)
plot(myPLN_torch$model_par$Sigma,
     myPLN_nlopt$model_par$Sigma); abline(0, 1)
myPLN_torch$loglik
myPLN_nlopt$loglik

system.time(myPLN_torch <-
              PLN(Abundance ~ 1  + offset(log(Offset)),
                  data = trichoptera, control = list(backend = "torch", covariance = "diagonal"))
)
system.time(myPLN_nlopt <-
              PLN(Abundance ~ 1  + offset(log(Offset)),
                  data = trichoptera, control = list(backend = "nlopt", covariance = "diagonal"))
)
plot(myPLN_torch$model_par$Theta,
     myPLN_nlopt$model_par$Theta); abline(0, 1)
plot(myPLN_torch$model_par$Sigma,
     myPLN_nlopt$model_par$Sigma); abline(0, 1)
myPLN_torch$loglik
myPLN_nlopt$loglik

data("mollusk")
mollusk <- prepare_data(mollusk$Abundance, mollusk$Covariate)
system.time(myPLN_torch <-
              PLN(Abundance ~ 1  + offset(log(duration)),
                  data = mollusk, control = list(backend = "torch"))
)
system.time(myPLN_nlopt <-
              PLN(Abundance ~ 1  + offset(log(duration)),
                  data = mollusk, control = list(backend = "nlopt"))
)

plot(myPLN_torch$model_par$Theta,
     myPLN_nlopt$model_par$Theta); abline(0, 1)
plot(myPLN_torch$model_par$Sigma,
     myPLN_nlopt$model_par$Sigma); abline(0, 1)
myPLN_torch$loglik
myPLN_nlopt$loglik




load("../BarentsFish.Rdata")
myPLN_nlopt <- PLN(count ~ -1 + covariates, data = Data)
myPLN_torch <- PLN(count ~ -1 + covariates, data = Data, control = list(backend = "torch"))

plot(myPLN_torch$model_par$Theta,
     myPLN_nlopt$model_par$Theta); abline(0, 1)
plot(myPLN_torch$model_par$Sigma,
     myPLN_nlopt$model_par$Sigma); abline(0, 1)
myPLN_torch$loglik
myPLN_nlopt$loglik

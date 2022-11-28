data("oaks")
system.time(myPLN_torch <-
              PLN(Abundance ~ 1  + offset(log(Offset)),
                  data = oaks, control = list(backend = "torch"))
)
system.time(myPLN_nlopt <-
              PLN(Abundance ~ 1  + offset(log(Offset)),
                  data = oaks, control = list(backend = "nlopt"))
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
                  data = mollusk, control = list(backend = "torch"))
)
system.time(myPLN_nlopt <-
              PLN(Abundance ~ 1  + offset(log(duration)),
                  data = mollusk, control = list(backend = "nlopt"))
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
                  data = barents, control = list(backend = "torch"))
)
system.time(myPLN_nlopt <-
              PLN(Abundance ~ 1 + Depth + Temperature + offset(log(Offset)),
                  data = barents, control = list(backend = "nlopt", xtol_rel = 1e-8))
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



library(PLNmodels)
data("trichoptera")

# Fit PLN on trichoptera
trichoptera = prepare_data(trichoptera$Abundance, trichoptera$Covariate)
model = PLN(Abundance ~ 1 + offset(log(Offset)), data = trichoptera)

# take training data
newdata = trichoptera[1,]

# try to do 1 VEstep
VE = model$optimize_vestep(covariates = matrix(1), offsets = matrix(rep(log(newdata$Offset), dim(newdata$Abundance)[2]), nrow = 1), responses = newdata$Abundance, weights = 1)
print(VE$M) #  0 (which corresponds to initialization)

# Try forcing (B, Omega) and tweek optimizer parameters
VE = model$optimize_vestep(covariates = matrix(1), offsets = matrix(rep(log(newdata$Offset), dim(newdata$Abundance)[2]), nrow = 1), responses = newdata$Abundance, weights = 1, 
                           B = model$model_par$B,
                           Omega = model$model_par$Omega,
                           control = PLN_param(backend = "nlopt",
                                               ### Comment/ uncomment to try to give Omega directly to the optimizer
                                               covariance = model$vcov_model,
                                               # covariance = "fixed",
                                               # Omega = matrix(c(1,2)), # wrong shape but still works, why ?
                                               # Omega = model$model_par$Omega,
                                               ### Try lowering optimizer threshold to avoid local minima
                                               config_optim = list(maxeval=1e6, ftol_rel=1e-32, xtol_rel=1e-32)))
print(VE$M) # still 0

# try with random integers
newdata = list( Abundance = matrix(c(3,0,1,0,0,0,0,0,0,1,0,1,2,0,2,14,1), nrow = 1), Offset = 1)
VE = model$optimize_vestep(covariates = matrix(1), offsets = matrix(rep(log(newdata$Offset), dim(newdata$Abundance)[2]), nrow = 1), responses = newdata$Abundance, weights = 1)
print(VE$M)

# try with higher order of magnitude for observations
newdata = list( Abundance = matrix(1:17, nrow = 1), Offset = 1)
VE = model$optimize_vestep(covariates = matrix(1), offsets = matrix(rep(log(newdata$Offset), dim(newdata$Abundance)[2]), nrow = 1), responses = newdata$Abundance, weights = 1,
                           control = PLN_param(config_optim = list(maxeval=1e2)))
print(VE$M) # not 0, and tweeking the config_optim params seems to have an effect on the result


library(PLNmodels)

data("trichoptera")
trichoptera <- prepare_data(trichoptera$Abundance, covariates = trichoptera$Covariate)
p <- ncol(trichoptera$Abundance)

PLN_full <- PLN(Abundance ~ 1 + offset(log(Offset)), data = trichoptera)
PLN_spherical <- PLN(Abundance ~ 1 + offset(log(Offset)), data = trichoptera, control = PLN_param(covariance = "spherical"))

PLN_full_block <- PLNblock(Abundance ~ 1 + offset(log(Offset)), nb_blocks = p, data = trichoptera)$models[[1]]
PLN_one_block  <- PLNblock(Abundance ~ 1 + offset(log(Offset)), nb_blocks = 1, data = trichoptera)$models[[1]]

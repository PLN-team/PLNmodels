library(PLNmodels)

## get oaks data set
load("inst/oaks_study/oaks_Ealphitoides.RData")

myPLNs <- PLN(Abundancies ~ 1 + log(Offsets), data = oaks)

system.time(covariates_tree <- PLN(Y ~ 1 + covariates$tree + covariates$orientation + offset(log(O))))


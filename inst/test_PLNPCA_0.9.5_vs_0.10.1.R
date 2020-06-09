## Need Shift + Ctrl + F10 first
rm(list=ls())
## get oaks data set
load("inst/case_studies/oaks_mildew/oaks_alphitoides.RData")
class(oaks$Abundancies) <- class(oaks$Abundancies)[-match("AsIs", class(oaks$Abundancies))]

library(PLNmodels, lib.loc = "~/R/x86_64-pc-linux-gnu-library/4.0/")
cat("\nVersion", as.character(packageVersion("PLNmodels")), "\n")

## simple PLN
# time_dev <- system.time(myPLN_dev <- PLNPCA(Abundancies ~ 1 + offset(log(sequencingEffort)), data = oaks, ranks = 1:25, control_main = list(cores = 10, trace = 0)))
time_dev <- system.time(myPLN_dev <- PLNPCA(Abundancies ~ 0 + treeStatus + offset(log(sequencingEffort)), data = oaks, control_main = list(trace = 0)))

detach("package:PLNmodels")

library(PLNmodels, lib.loc = "~/R/x86_64-pc-linux-gnu-library/PLNmodels_9.5/")
cat("\nVersion", as.character(packageVersion("PLNmodels")), "\n")

## simple PLN
time_CRAN <- system.time(myPLN_CRAN <- PLNPCA(Abundancies ~ 0 + treeStatus + offset(log(sequencingEffort)), data = oaks, control_main = list(trace = 0)))


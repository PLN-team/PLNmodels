library(PLNmodels)

## get oaks data set
load("inst/case_studies/oaks_mildew/oaks_alphitoides.RData")
class(oaks$Abundancies) <- class(oaks$Abundancies)[-match("AsIs", class(oaks$Abundancies))]

## simple PLN
system.time(myPLN <- PLN(Abundancies ~ 0 + treeStatus + offset(log(sequencingEffort)), data = oaks))
system.time(myPLN_diagonal <- PLN(Abundancies ~ 0 + treeStatus + offset(log(sequencingEffort)), data = oaks, control = list(covariance = "diagonal")))
system.time(myPLN_spherical <- PLN(Abundancies ~ 0 + treeStatus + offset(log(sequencingEffort)), data = oaks, control = list(covariance = "spherical")))

## Discriminant Analysis with LDA
myLDA_tree <- PLNLDA(Abundancies ~ 1 + offset(log(sequencingEffort)), grouping = oaks$treeStatus, data = oaks)
plot(myLDA_tree)
plot(myLDA_tree, "individual")

myLDA_tree_diagonal <- PLNLDA(Abundancies ~ 1 + offset(log(sequencingEffort)), grouping = oaks$treeStatus, data = oaks, control = list(covariance = "diagonal"))
plot(myLDA_tree_diagonal)
otu.family <- factor(rep(c("fungi", "E. aphiltoides", "bacteria"), c(47, 1, 66)))
plot(myLDA_tree, "variable", var_cols = otu.family) ## TODO: add color for arrows to check

## One dimensional check of plot
myLDA_orientation <- PLNLDA(Abundancies ~ 1 + offset(log(sequencingEffort)), grouping = oaks$orientation, data = oaks)
plot(myLDA_orientation)

## Dimension reduction with PCA
system.time(myPLNPCAs <- PLNPCA(Abundancies ~ 1 + offset(log(sequencingEffort)), data = oaks, ranks = 1:30, control_main = list(cores = 10))) # about 55 sec.
myPLNPCA <-
plot(getBestModel(myPLNPCAs), ind_cols = oaks$treeStatus)

## Network inference with sparce covariance estimation
system.time(myPLNnets <- PLNnetwork(Abundancies ~ 0 + treeStatus + offset(log(sequencingEffort)), data = oaks, control_main = list(trace = 2)))
stability_selection(myPLNnets)
plot(getBestModel(myPLNnets, "StARS", stability = .985))


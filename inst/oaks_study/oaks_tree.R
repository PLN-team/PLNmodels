library(PLNmodels)

## get oaks data set
load("inst/oaks_study/oaks_alphitoides.RData")

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
plot(myLDA_tree, "variable") ## TODO: add color for arrows to check

## One dimensional check of plot
myLDA_orientation <- PLNLDA(Abundancies ~ 1 + offset(log(sequencingEffort)), grouping = oaks$orientation, data = oaks)
plot(myLDA_orientation)

## Dimension reduction with PCA
system.time(myPLNPCAs <- PLNPCA(Abundancies ~ 1 + offset(log(sequencingEffort)), data = oaks, ranks = 1:30)) # about 250 sec.
myPLNPCA <-
plot(getBestModel(myPLNPCAs), ind_cols = oaks$treeStatus)

## Network inference with sparce covariance estimation
myPLNnets <- PLNnetwork(Abundancies ~ 1 + treeStatus + offset(log(sequencingEffort)), data = oaks)
stability_selection(myPLNnets)
plot(getBestModel(myPLNnets, "StARS", stability = .985))


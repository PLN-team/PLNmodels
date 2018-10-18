library(PLNmodels)

## get oaks data set
load("inst/oaks_study/oaks_alphitoides.RData")

## simple PLN
system.time(myPLN <- PLN(Abundancies ~ 1 + offset(log(sequencingEffort)), data = oaks, control = list(maxtime = 1)))
system.time(myPLN_diagonal <- PLN(Abundancies ~ 1 + offset(log(sequencingEffort)), data = oaks, control = list(covariance = "diagonal")))
system.time(myPLN_spherical <- PLN(Abundancies ~ 1 + offset(log(sequencingEffort)), data = oaks, control = list(covariance = "spherical")))

## Discriminant Analysis with LDA
myLDA_tree <- PLNLDA(Abundancies ~ 1 + offset(log(sequencingEffort)), grouping = oaks$treeStatus, data = oaks)
myLDA_tree$plot_LDA()

myLDA_tree_spherical <- PLNLDA(Abundancies ~ 1 + offset(log(sequencingEffort)), grouping = oaks$treeStatus, data = oaks, control = list(covariance = "diagonal"))
myLDA_tree_spherical $plot_LDA()

myLDA_orientation <- PLNLDA(Abundancies ~ 1 + offset(log(sequencingEffort)), grouping = oaks$orientation, data = oaks)
myLDA_orientation$plot_LDA()

## Dimension reduction with PCA
system.time(myPLNPCAs <- PLNPCA(Abundancies ~ 1 + offset(log(sequencingEffort)), data = oaks, ranks = 1:30)) # about 50 sec.
system.time(myPLNPCAs_spherical <- PLNPCA(Abundancies ~ 1 + offset(log(sequencingEffort)), data = oaks, ranks = 1:30, control_init = list(covariance = "spherical"))) # about 50 sec.

myPLNPCA <- myPLNPCAs$getBestModel('ICL')
myPLNPCA$plot_PCA(cols.ind = oaks$treeStatus)

## Network inference with sparce covariance estimation
myPLNnets <- PLNnetwork(Abundancies ~ 1 + treeStatus + offset(log(sequencingEffort)), data = oaks)
myPLNnets$stability_selection()
myPLNnet <- myPLNnets$getBestModel("StARS", .9825)
myPLNnet$plot_network()

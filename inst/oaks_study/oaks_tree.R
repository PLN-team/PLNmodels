library(PLNmodels)

## get oaks data set
load("inst/oaks_study/oaks_alphitoides.RData")

## simple PLN
myPLN <- PLN(Abundancies ~ 1 + offset(log(sequencingEffort)), data = oaks)

## Discriminant Analysis with LDA
myLDA_tree <- PLNLDA(Abundancies ~ 1 + offset(log(sequencingEffort)), grouping = oaks$treeStatus, data = oaks)
myLDA_tree$plot_LDA()

myLDA_orientation <- PLNLDA(Abundancies ~ 1 + offset(log(sequencingEffort)), grouping = oaks$orientation, data = oaks)
myLDA_orientation$plot_LDA()

## Dimension reduction with PCA
myPLNPCAs <- PLNPCA(Abundancies ~ 1 + offset(log(sequencingEffort)), data = oaks, ranks = 1:30, control.main = list(cores = 10)) # about 50 sec.
myPLNPCA <- myPLNPCAs$getBestModel('ICL')
myPLNPCA$plot_PCA(cols.ind = oaks$treeStatus)

## Network inference with sparce covariance estimation
myPLNnets <- PLNnetwork(Abundancies ~ 1 + treeStatus + offset(log(sequencingEffort)), data = oaks)
myPLNnets$stability_selection()
myPLNnet <- myPLNnets$getBestModel("StARS", .9925)
myPLNnet$plot_network()

## Clustering with mixture model
myPLNMMs <- PLNMM(Abundancies ~ 1 + offset(log(sequencingEffort)), data = oaks, clusters = 1:5)
plot(purrr::map_dbl(myPLNMMs$models, ~aricode::NID(.$memberships, oaks$treeStatus)),
     xlab = "number of components", ylab = "NID",  type = "l")


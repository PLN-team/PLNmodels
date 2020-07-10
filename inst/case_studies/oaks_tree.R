library(PLNmodels)

nb_cores <- 10

## get oaks data set
data(oaks)

## simple PLN
system.time(myPLN <- PLN(Abundance ~ 0 + tree + offset(log(Offset)), data = oaks))
system.time(myPLN_diagonal <- PLN(Abundance ~ 0 + tree + offset(log(Offset)), data = oaks, control = list(covariance = "diagonal")))
system.time(myPLN_spherical <- PLN(Abundance ~ 0 + tree + offset(log(Offset)), data = oaks, control = list(covariance = "spherical")))

rbind(
  myPLN$criteria,
  myPLN_diagonal$criteria,
  myPLN_spherical$criteria
) %>%
  as.data.frame(row.names = c("full", "diagonal", "spherical")) %>%
  knitr::kable()

## Discriminant Analysis with LDA
myLDA_tree <- PLNLDA(Abundance ~ 1 + offset(log(Offset)), grouping = tree, data = oaks)
plot(myLDA_tree)
plot(myLDA_tree, "individual")

myLDA_tree_diagonal <- PLNLDA(Abundance ~ 1 + offset(log(Offset)), grouping = tree, data = oaks, control = list(covariance = "diagonal"))
plot(myLDA_tree_diagonal)
otu.family <- factor(rep(c("fungi", "E. aphiltoides", "bacteria"), c(47, 1, 66)))
plot(myLDA_tree, "variable", var_cols = otu.family) ## TODO: add color for arrows to check

## One dimensional check of plot
myLDA_orientation <- PLNLDA(Abundance ~ 1 + offset(log(Offset)), grouping = orientation, data = oaks)
plot(myLDA_orientation)

## Dimension reduction with PCA
system.time(myPLNPCAs <- PLNPCA(Abundance ~ 1 + offset(log(Offset)), data = oaks, ranks = 1:30, control_main = list(cores = nb_cores))) # about 40 sec.
plot(myPLNPCAs)
myPLNPCA <- getBestModel(myPLNPCAs)
plot(myPLNPCA, ind_cols = oaks$tree)

## Network inference with sparce covariance estimation
system.time(myPLNnets <- PLNnetwork(Abundance ~ 0 + tree + offset(log(Offset)), data = oaks, control_main = list(trace = 2)))
stability_selection(myPLNnets, mc.cores = nb_cores)
plot(getBestModel(myPLNnets, "StARS", stability = .985))


library(PLNmodels)

## get oaks data set
data(mollusk)
mollusc <- prepare_data(mollusk$Abundance, mollusk$Covariate)#> Warning: Sample(s) 134 and 137 and 145 and 146 dropped for lack of positive counts.

## simple PLN
system.time(myPLN_M0 <- PLN(Abundance ~ 1 + offset(log(Offset)), data = mollusc))
system.time(myPLN <- PLN(Abundance ~ 0 + site + offset(log(Offset)), data = mollusc))
system.time(myPLN_diagonal <- PLN(Abundance ~ 0 + site + offset(log(Offset)), data = mollusc, control = list(covariance = "diagonal")))
system.time(myPLN_spherical <- PLN(Abundance ~ 0 + site + offset(log(Offset)), data = mollusc, control = list(covariance = "spherical")))

rbind(
  myPLN_M0$criteria,
  myPLN$criteria,
  myPLN_diagonal$criteria,
  myPLN_spherical$criteria
) %>%
  as.data.frame(row.names = c("full", "full_site", "diagonal_site", "spherical_site")) %>%
  knitr::kable()

## Discriminant Analysis with LDA
myLDA_site <- PLNLDA(Abundance ~ season + offset(log(Offset)), grouping = site, data = mollusc)
plot(myLDA_site)

myLDA_season <- PLNLDA(Abundance ~ site + offset(log(Offset)), grouping = season, data = mollusc)
plot(myLDA_season)


## Dimension reduction with PCA
system.time(myPLNPCAs <- PLNPCA(Abundance ~ 1 + offset(log(Offset)), data = mollusc, ranks = 1:10, control_main = list(cores = 10))) # about 40 sec.
plot(myPLNPCAs)
myPLNPCA <- getBestModel(myPLNPCAs)
plot(myPLNPCA, ind_cols = mollusc$site)

## Dimension reduction with PCA
system.time(myPLNPCAs_site <- PLNPCA(Abundance ~ 0 + site + offset(log(Offset)), data = mollusc, ranks = 1:6, control_main = list(cores = 6))) # about 40 sec.
plot(myPLNPCAs_site)
myPLNPCA_site <- getBestModel(myPLNPCAs_site)
plot(myPLNPCA_site, ind_cols = mollusc$method)


## Network inference with sparce covariance estimation
myPLNnets <- PLNnetwork(Abundance ~ 0 + treeStatus + offset(log(sequencingEffort)), data = oaks)
stability_selection(myPLNnets, mc.cores = 10)
plot(getBestModel(myPLNnets, "StARS", stability = .985))


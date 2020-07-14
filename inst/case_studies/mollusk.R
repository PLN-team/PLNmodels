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
myPLNPCA_ICL <- getBestModel(myPLNPCAs, crit = 'ICL')
myPLNPCA_BIC <- getBestModel(myPLNPCAs, crit = 'BIC')
plot(myPLNPCA_ICL, ind_cols = mollusc$site)
plot(myPLNPCA_BIC, ind_cols = mollusc$site)

## Dimension reduction with PCA
system.time(myPLNPCAs_site <- PLNPCA(Abundance ~ 0 + site + offset(log(Offset)), data = mollusc, ranks = 1:6, control_main = list(cores = 6))) # about 40 sec.
plot(myPLNPCAs_site)
myPLNPCA_site <- getBestModel(myPLNPCAs_site, "ICL")
plot(myPLNPCA_site, ind_cols = mollusc$season)




## clustering
my_mixtures <-  PLNmixture(Abundance ~ 1 + offset(log(Offset)), clusters = 1:5, data = mollusc, control_main = list(covariance = "diagonal", core = nb_cores))

plot(my_mixtures)
my_mixtures$plot_objective()

myPLN <- getBestModel(my_mixtures)
myPLN$plot_clustering_pca()
myPLN$plot_clustering_data()

aricode::ARI(myPLN$memberships, mollusc$site)

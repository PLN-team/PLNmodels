library(PLNmodels)
library(factoextra)

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

# fancy graph with factoextra
factoextra::fviz_pca_biplot(
  myPLNPCA, select.var = list(contrib = 10), addEllipses = TRUE, habillage = oaks$tree,
  title = "Biplot (10 most contributing species, samples colored by susceptibility)"
  ) + labs(col = "distance (cm)") + scale_color_viridis_d()

## Dimension reduction with PCA
system.time(myPLNPCAs_tree <- PLNPCA(Abundance ~ 0 + tree + offset(log(Offset)), data = oaks, ranks = 1:30, control_main = list(cores = nb_cores))) # about 40 sec.
plot(myPLNPCAs_tree)
myPLNPCA_tree <- getBestModel(myPLNPCAs_tree)

# fancy graph with factoextra
factoextra::fviz_pca_biplot(
  myPLNPCA_tree, select.var = list(contrib = 10), col.ind  = oaks$distTOground,
  title = "Biplot after correction (10 most contributing species, samples colored by distance to ground)") +
  labs(col = "distance (cm)") + scale_color_viridis_c()

factoextra::fviz_pca_ind(myPLNPCA_tree, axes = c(1,2), col.ind = oaks$distTOground)
factoextra::fviz_pca_var(myPLNPCA_tree, axes = c(1,2), select.var = list(contrib = 10))

## Network inference with sparce covariance estimation
system.time(myPLNnets <- PLNnetwork(Abundance ~ 0 + tree + offset(log(Offset)), data = oaks, control_main = list(trace = 2)))
stability_selection(myPLNnets, mc.cores = nb_cores)
plot(getBestModel(myPLNnets, "StARS", stability = .985))

## Mixture model to recover tree structure
system.time(my_mixtures <- PLNmixture(Abundance ~ 1 + offset(log(Offset)), data = oaks, clusters = 1:7, control_main = list(cores = nb_cores)))

plot(my_mixtures, criteria = c("loglik", "ICL", "BIC"))

myPLN <- my_mixtures %>% getModel(3)

plot(myPLN, "pca", main = 'clustering memberships in individual factor map')
plot(myPLN, "matrix")

aricode::ARI(myPLN$memberships, oaks$tree)
table(myPLN$memberships, oaks$tree)

data.frame(
  nb_components  = sapply(my_mixtures$models, function(model) model$k),
  ARI = sapply(lapply(my_mixtures$models, function(model) model$memberships), aricode::ARI, oaks$tree),
  AMI = sapply(lapply(my_mixtures$models, function(model) model$memberships), aricode::AMI, oaks$tree),
  NID = sapply(lapply(my_mixtures$models, function(model) model$memberships), aricode::NID, oaks$tree)
) %>%
  tidyr::pivot_longer(-nb_components,names_to = "score") %>%
  dplyr::group_by(score) %>%
  ggplot(aes(x = nb_components, y = value, colour = score)) + geom_line() + theme_bw() + labs(y = "clustering similarity", x = "number of components")

## Mixture model to recover tree structure - with covariates
system.time(my_mixtures <- PLNmixture(Abundance ~ 0 + tree + distTOground + offset(log(Offset)), data = oaks, clusters = 1:5, control_main = list(cores = nb_cores)))

plot(my_mixtures, criteria = c("loglik", "ICL", "BIC"))

myPLN <- my_mixtures %>% getBestModel()

myPLN$plot_clustering_pca(main = 'clustering memberships in individual factor map')
p <- myPLN$plot_clustering_data()

aricode::ARI(myPLN$memberships, oaks$tree)


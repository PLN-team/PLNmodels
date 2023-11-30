library(PLNmodels)
library(factoextra)

## setting up future for parallelism
nb_cores <- 10
options(future.fork.enable = TRUE)
future::plan("multicore", workers = nb_cores)

## get oaks data set
data(stoc)

## simple PLN
system.time(myPLN_diag_0 <- PLN(Abundance ~ 1 + offset(log(Offset)), data = stoc, control = PLN_param(covariance = "diagonal")))
system.time(myPLN_diag_temp <- PLN(Abundance ~ 1 + temp + offset(log(Offset)), data = stoc, control = PLN_param(covariance = "diagonal")))
system.time(myPLN_diag_precip <- PLN(Abundance ~ 1 + precip + offset(log(Offset)), data = stoc, control = PLN_param(covariance = "diagonal", config_optim = list(xtol_rel=1e-8))))
system.time(myPLN_diag_div <- PLN(Abundance ~ 1 + div + offset(log(Offset)), data = stoc, control = PLN_param(covariance = "diagonal")))
system.time(myPLN_diag_zonebio <- PLN(Abundance ~ 0 + zonebio + offset(log(Offset)), data = stoc, control = PLN_param(covariance = "diagonal")))
system.time(myPLN_diag_Agricultural <- PLN(Abundance ~ 1 + cover_Agricultural + offset(log(Offset)), data = stoc, control = PLN_param(covariance = "diagonal")))
system.time(myPLN_diag_Forest <- PLN(Abundance ~ 1 + cover_Forest + offset(log(Offset)), data = stoc, control = PLN_param(covariance = "diagonal")))
system.time(myPLN_diag_Artificial <- PLN(Abundance ~ 1 + cover_Artificial + offset(log(Offset)), data = stoc, control = PLN_param(covariance = "diagonal")))

cover <- dplyr::select(stoc, starts_with('cover')) %>%  scale(TRUE, FALSE)
## Loadings/rotation matrix
U <- eigen(cov(cover))$vectors
## Function for projection
stoc$cover_PC1 <- cover %*% U[, 1, drop = FALSE]
system.time(myPLN_diag_cover <- PLN(Abundance ~ 1 + cover_PC1 + offset(log(Offset)), data = stoc, control = PLN_param(covariance = "diagonal", config_optim = list(xtol_rel = 1e-8))))

system.time(myPLN_diag_temp_Forest <- PLN(Abundance ~ 1 + cover_Forest + temp + offset(log(Offset)), data = stoc, control = PLN_param(covariance = "diagonal")))

rbind(
  myPLN_diag_0$criteria,
  myPLN_diag_temp$criteria,
  myPLN_diag_precip$criteria,
  myPLN_diag_div$criteria,
  myPLN_diag_zonebio$criteria,
  myPLN_diag_Agricultural$criteria,
  myPLN_diag_Forest$criteria,
  myPLN_diag_Artificial$criteria,
  myPLN_diag_cover$criteria,
  myPLN_diag_temp_Forest$criteria
) %>%
  as.data.frame(row.names = c("1", "temp", "precip", "div", "zonebio", "Agricultural", "Forest", "Artificial", "cover", "temp-Forest")) %>%
  dplyr::arrange(loglik) %>%
  knitr::kable()

## --> best model is temp + cover_Forest

myPLN_temp_Forest <- PLN(Abundance ~ 1 + cover_Forest + temp + offset(log(Offset)), data = stoc)

nb_blocks <- seq(20, 80, by=4)
system.time(myPLN_blocks <- PLNblock(Abundance ~ 1 + cover_Forest + temp + offset(log(Offset)), nb_blocks = nb_blocks, data = stoc, control = PLNblock_param(inception = myPLN_temp_Forest)))
myPLN_block <- getBestModel(myPLN_blocks)

data.frame(
  fitted   = c(as.vector(fitted(myPLN)), as.vector(fitted(myPLN_diagonal)),
               as.vector(fitted(myPLN_spherical)), as.vector(fitted(myPLN_block))),
  observed = rep(as.vector(stoc$Abundance), 2),
  method   = factor(rep(c("full", "diagonal", "spherical", "block"), each = length(stoc$Abundance)))
) %>%
  ggplot(aes(x = observed, y = fitted)) +
  geom_point(size = .5, alpha =.25 ) +
  facet_wrap( ~ method) +
  scale_x_log10() +
  scale_y_log10(limits  =c(1e-2, 1e4)) +
  theme_bw() + annotation_logticks()

rbind(
  myPLN$criteria,
  myPLN_diagonal$criteria,
  myPLN_spherical$criteria,
  myPLN_block$criteria
) %>%
  as.data.frame(row.names = c("full", "diagonal", "spherical", "block")) %>%
  knitr::kable()

## Discriminant Analysis with LDA
myLDA_zonebio <- PLNLDA(Abundance ~ 1 + offset(log(Offset)), grouping = zonebio, data = stoc)
plot(myLDA_zonebio)
plot(myLDA_zonebio, "individual")

## Dimension reduction with PCA
system.time(myPLNPCAs <- PLNPCA(Abundance ~ cover_Agriculture + offset(log(Offset)), data = stoc, ranks = 1:30)) # about 40 secs.
plot(myPLNPCAs)
myPLNPCA <- getBestModel(myPLNPCAs)

# fancy graph with factoextra
factoextra::fviz_pca_biplot(
  myPLNPCA, select.var = list(contrib = 10), col.ind = stoc$cover_Forest,
  title = "Biplot (10 most contributing species)"
) + labs(col = "cover") + scale_color_viridis_c()

## Dimension reduction with PCA
system.time(myPLNPCAs_zonebio <- PLNPCA(Abundance ~ 0 + zonebio + offset(log(Offset)), data = stoc, ranks = 1:30)) # about 40 sec.
plot(myPLNPCAs_tree)
myPLNPCA_tree <- getBestModel(myPLNPCAs_tree)

## Network inference with sparce covariance estimation
system.time(myPLNnets <- PLNnetwork(Abundance ~ 0 + tree + offset(log(Offset)), data = oaks, control = PLNnetwork_param(min_ratio = 0.1, penalize_diagonal = FALSE)))
plot(myPLNnets)
plot(getBestModel(myPLNnets, "EBIC"))
stability_selection(myPLNnets)
plot(getBestModel(myPLNnets, "StARS", stability = .975))

future::plan("sequential")

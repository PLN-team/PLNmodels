library(PLNmodels)
library(profvis)
data("trichoptera")

TotalCount <- matrix(
  rowSums(trichoptera$Abundance),
  nrow = nrow(trichoptera$Abundance),
  ncol = ncol(trichoptera$Abundance)
)

system.time(model1 <- PLNPCA(Abundance ~ 1 + offset(log(TotalCount)), data = trichoptera, ranks = 1:8))
system.time(model2 <- PLNPCA(Abundance ~ 1 + offset(log(TotalCount)), data = trichoptera, ranks = 1:8, control_init = list(covariance = "spherical")))

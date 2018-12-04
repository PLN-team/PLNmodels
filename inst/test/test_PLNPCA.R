library(PLNmodels)
library(profvis)
data("trichoptera")

TotalCount <- matrix(
  rowSums(trichoptera$Abundance),
  nrow = nrow(trichoptera$Abundance),
  ncol = ncol(trichoptera$Abundance)
)

system.time(models <- PLNPCA(Abundance ~ 1 + offset(log(TotalCount)), data = trichoptera, ranks = 1:8))

plot(getBestModel(models), ind_cols = trichoptera$Group)
plot(getBestModel(models), "variable")

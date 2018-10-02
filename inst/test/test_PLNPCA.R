library(PLNmodels)
library(profvis)
data("trichoptera")

TotalCount <- matrix(
  rowSums(trichoptera$Abundance),
  nrow = nrow(trichoptera$Abundance),
  ncol = ncol(trichoptera$Abundance)
)

profvis(model <- PLNPCA(Abundance ~ 1 + offset(log(TotalCount)), data = trichoptera, ranks = 1:8))

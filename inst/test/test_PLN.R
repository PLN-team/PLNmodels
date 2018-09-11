library(PLNmodels)
library(profvis)

data("trichoptera")

profvis::profvis(model <- PLN(abundance ~ 1, data = trichopetra))

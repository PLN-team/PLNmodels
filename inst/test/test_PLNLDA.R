library(PLNmodels)
library(profvis)
data("trichometeo")

profvis(model <- PLNLDA(Abundance ~ 0, grouping = trichoptera$Group, data = trichoptera))


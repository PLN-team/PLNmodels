library(PLNmodels)
library(ggplot2)
library(gridExtra)
library(profvis)
data("trichoptera")

## warm / cold start are very similar in terms of timings
# timings <- microbenchmark::microbenchmark(
#   warm  = PLNnetwork(Abundance ~ 1, data = trichoptera, control.main = list(warm = TRUE , trace = 0)),
#   cold  = PLNnetwork(Abundance ~ 1, data = trichoptera, control.main = list(warm = FALSE, trace = 0)),
#   times = 20)

profvis(fits <- PLNnetwork(Abundance ~ 1, data = trichoptera))

fits$stability_selection()

p0 <- fits$plot()
p1 <- fits$plot_objective()
p2 <- ggplot(fits$coefficient_path()) +
  aes(x = Penalty, y = Coeff, group = Edge, color = Edge) +
  geom_line(show.legend = FALSE) + theme_bw() +
  coord_trans(x = "log10")
p3 <- ggplot(fits$stability_path) +
  aes(x = Penalty, y = Prob, group = Edge, color = Edge) +
  geom_line(show.legend = FALSE) + theme_bw() +
  coord_trans(x = "log10")
p4 <- fits$plot_stars(.95)

grid.arrange(p0, p1, p2, p4, ncol = 2)

network_BIC    <- fits$getBestModel("BIC")
network_EBIC   <- fits$getBestModel("BIC")
network_StARS  <- fits$getBestModel("StARS", stability = 0.95)

par(mfrow = c(1,2))
network_BIC$plot_network(type = "partial_cor")
network_StARS$plot_network()

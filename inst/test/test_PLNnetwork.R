library(PLNmodels)
library(ggplot2)
library(gridExtra)
library(testthat)
library(profr)
library(ade4)
data("trichometeo")

abundance <- as.matrix(trichometeo$fau)

profiling1 <- profr::profr(fits     <- PLNnetwork(abundance))
profiling2 <- profr::profr(fits_alt <- PLNnetwork(abundance ~ 1))

expect_equivalent(fits, fits_alt)

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
network_BIC$plot_network()
network_StARS$plot_network()

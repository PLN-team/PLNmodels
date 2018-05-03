library(PLNmodels)
library(ggplot2)
library(gridExtra)
library(testthat)
library(profr)
library(ade4)
data("trichometeo")

abundance <- as.matrix(trichometeo$fau)

profiling1 <- profr::profr(model1 <- PLNnetwork(abundance))
profiling2 <- profr::profr(model2 <- PLNnetwork(abundance ~ 1))

expect_equivalent(model1, model2)

p0 <- model1$plot()
p1 <- model1$plot_objective()
p2 <- ggplot(model1$coefficient_path()) +
  aes(x = Penalty, y = Coeff, group = Edge, color = Edge) +
  geom_line(show.legend = FALSE) +
  coord_trans(x = "log10")

grid.arrange(p0, p1, p2, ncol = 1)
model1$getBestModel()$plot_network("partial_cor")

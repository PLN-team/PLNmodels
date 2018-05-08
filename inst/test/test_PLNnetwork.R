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

model1$stability_selection()

p0 <- model1$plot()
p1 <- model1$plot_objective()
p2 <- ggplot(model1$coefficient_path()) +
  aes(x = Penalty, y = Coeff, group = Edge, color = Edge) +
  geom_line(show.legend = FALSE) + theme_bw() +
  coord_trans(x = "log10")
p22 <- ggplot(model1$stability_path) +
  aes(x = Penalty, y = Prob, group = Edge, color = Edge) +
  geom_line(show.legend = FALSE) + theme_bw() +
  coord_trans(x = "log10")
p3 <- model1$plot_stars()

grid.arrange(p0, p1, p2, p3, ncol = 2)

model1$getBestModel()$plot_network("partial_cor")

library(PLNmodels)
library(ggplot2)
library(gridExtra)
library(profr)
library(ade4)
data("trichometeo")
abundance <- as.matrix(trichometeo$fau)

profiling <- profr::profr(model <- PLNnetwork(abundance ~ 1))
plot(profiling)

p0 <- model$plot()
p1 <- model$plot_objective()
p2 <- ggplot(model$coefficient_path(), aes(x=Penalty, y=Coeff, group = Edge, color = Edge)) + geom_line(show.legend = FALSE) + coord_trans(x = "log10")

grid.arrange(p0, p1, p2, ncol=1)
model$getBestModel()$plot_network()

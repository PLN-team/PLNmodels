library(PLNmodels)
library(ggplot2)
library(gridExtra)
library(profr)
library(ade4)
data("trichometeo")
abundance <- as.matrix(trichometeo$fau)

profiling_1 <- profr::profr(model_1 <- PLNnetwork(abundance ~ 1, control.init = list(nPenalties = 50, min.ratio = 1e-2), control.main = list(nloptr = TRUE)))
profiling_2 <- profr::profr(model_2 <- PLNnetwork(abundance ~ 1, control.init = list(nPenalties = 50, min.ratio = 1e-2), control.main = list()))

p01 <- model_1$plot()
p02 <- model_2$plot()

p11 <- model_1$plot_objective()
p12 <- model_2$plot_objective()

p21 <- ggplot(model_1$coefficient_path(), aes(x=Penalty, y=Coeff, group = Edge, color = Edge)) + geom_line(show.legend = FALSE) + coord_trans(x = "log10")
p22 <- ggplot(model_2$coefficient_path(), aes(x=Penalty, y=Coeff, group = Edge, color = Edge)) + geom_line(show.legend = FALSE) + coord_trans(x = "log10")

grid.arrange(p01, p02, p11, p12, p21, p22, ncol=2)

par(mfrow=c(2,2))
plot(profiling_1)
plot(profiling_2)
model_1$getBestModel()$plot_network("partial_cor")
model_2$getBestModel()$plot_network("partial_cor")

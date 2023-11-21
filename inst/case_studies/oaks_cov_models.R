library(PLNmodels)
library(tidyverse)

## setting up future for parallelism
nb_cores <- 20
options(future.fork.enable = TRUE)
future::plan("multicore", workers = nb_cores)

## get oaks data set
data(oaks)

## simple PLN
system.time(myPLN <- PLN(Abundance ~ 0 + tree + offset(log(Offset)), data = oaks))
system.time(myPLN_diagonal <- PLN(Abundance ~ 0 + tree + offset(log(Offset)), data = oaks, control = PLN_param(covariance = "diagonal")))
system.time(myPLN_spherical <- PLN(Abundance ~ 0 + tree + offset(log(Offset)), data = oaks, control = PLN_param(covariance = "spherical")))

## Blockwise covariance
system.time(myPLN_blocks <- PLNblock(Abundance ~ 0 + tree + offset(log(Offset)), nb_blocks = 1:10, data = oaks, control = PLNblock_param(inception = myPLN)))
myPLN_block <- getBestModel(myPLN_blocks)

data.frame(
  fitted   = c(as.vector(fitted(myPLN)), as.vector(fitted(myPLN_diagonal)),
               as.vector(fitted(myPLN_spherical)), as.vector(fitted(myPLN_block))),
  observed = rep(as.vector(oaks$Abundance), 2),
  method   = factor(rep(c("full", "diagonal", "spherical", "block"), each = length(oaks$Abundance)))
) %>%
  ggplot(aes(x = observed, y = fitted)) +
  geom_point(size = .5, alpha =.25 ) +
  facet_wrap( ~ method) +
  scale_x_log10() +
  scale_y_log10(limits  =c(1e-2, 1e4)) +
  theme_bw() + annotation_logticks()

rbind(
  myPLN$criteria,
  myPLN_diagonal$criteria,
  myPLN_spherical$criteria,
  myPLN_block$criteria
) %>%
  as.data.frame(row.names = c("full", "diagonal", "spherical", "block")) %>%
  knitr::kable()

library(PLNmodels)
library(tidyverse)

## setting up future for parallelism
nb_cores <- 20
options(future.fork.enable = TRUE)
future::plan("multicore", workers = nb_cores)

## get oaks data set
data(oaks)

## simple PLN
system.time(myPLN <- PLN(Abundance ~ 1 + offset(log(Offset)), data = oaks))

## Blockwise covariance
system.time(myPLN_blocks <- PLNblock(Abundance ~ 1 + offset(log(Offset)),
                                     nb_blocks = 1:114, data = oaks,
                                     control = PLNblock_param(inception = myPLN)))
plot(myPLN_blocks)
myPLN_blocks$plot_objective()

best_block <- getBestModel(myPLN_blocks)
list_Tau <- best_block$optim_par$posteriorProb

clusterings <- map(list_Tau, ~apply(., 2, which.max))

future::plan("sequential")

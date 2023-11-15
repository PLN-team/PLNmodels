library(PLNmodels)

## setting up future for parallelism
nb_cores <- 10
options(future.fork.enable = TRUE)
future::plan("multicore", workers = nb_cores)

data("trichoptera")
trichoptera <- prepare_data(trichoptera$Abundance, covariates = trichoptera$Covariate)
p <- ncol(trichoptera$Abundance)

PLN_full <- PLN(Abundance ~ 1 + offset(log(Offset)), data = trichoptera)
PLN_spherical <- PLN(Abundance ~ 1 + offset(log(Offset)), data = trichoptera, control = PLN_param(covariance = "spherical"))

PLN_full_block <- PLNblock(Abundance ~ 1 + offset(log(Offset)), nb_blocks = p, data = trichoptera)$models[[1]]
PLN_one_block  <- PLNblock(Abundance ~ 1 + offset(log(Offset)), nb_blocks = 1, data = trichoptera)$models[[1]]

PLN_blocks <- PLNblock(Abundance ~ 1 + offset(log(Offset)), nb_blocks = 1:p, data = trichoptera)
plot(PLN_blocks)


rbind(
  PLN_full$criteria,
  PLN_spherical$criteria,
  PLN_full_block$criteria,
  PLN_one_block$criteria
) %>%
  as.data.frame(row.names = c("full", "spherical", "p block", "1 block")) %>%
  knitr::kable()

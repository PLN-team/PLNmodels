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

## Check PLN full and PLN p block are similar (up to the term sum(tau * log(alpha)) )
vll_PLN_block <- PLN_full_block$loglik_vec - sum(PLN_full_block$posteriorProb * log(PLN_full_block$groupProportion))
vll_PLN <- PLN_full$loglik_vec
plot(vll_PLN, vll_PLN_block)
abline(0,1)
## YES
print(PLN_full$loglik - (PLN_full_block$loglik - sum(PLN_full_block$posteriorProb * log(PLN_full_block$groupProportion))))


## Check PLN full and PLN p block are similar (up to the term sum(tau * log(alpha)) )
vll_PLN_one_block <- PLN_one_block$loglik_vec - sum(PLN_one_block$posteriorProb * log(PLN_one_block$groupProportion))
vll_PLN_spherical <- PLN_spherical$loglik_vec
plot(vll_PLN_spherical, vll_PLN_one_block)
abline(0,1)
## YES
print(PLN_full$loglik - (PLN_full_block$loglik - sum(PLN_full_block$posteriorProb * log(PLN_full_block$groupProportion))))



library(PLNmodels); library(mvtnorm)


data("barents")
n <- nrow(barents$Abundance); p <- ncol(barents$Abundance)
fit <- PLN(Abundance ~ Latitude + Longitude + Depth + Temperature, data = barents)
X <- model.matrix(Abundance ~ Latitude + Longitude + Depth + Temperature, data = barents)

# Simulate data
n <- 1e2
X <- X[sample(1:nrow(X), n, replace=TRUE), ]
B0 <- fit$model_par$B; Sigma0 <- fit$model_par$Sigma
Y <- matrix(rpois(n*p, exp(X %*% B0 + mvtnorm::rmvnorm(n, sigma= Sigma0))), n, p)
d <- ncol(X)

# Jacknife sd estimate
tictoc::tic()
pln <- PLN(Y ~ -1 + X, control=PLN_param(config_post = list(jackknife=TRUE)))
tictoc::toc()

betaSdPackage <- as.vector(standard_error(pln, type='jackknife'))
jackknife_estimates <- attr(pln$model_par$B, "estimates_jaccknife") |> purrr::map("B")
betaJKPackage <- matrix(NA, n, p*d)
for(i in 1:n){
  betaJKPackage[i, ]  <- as.vector(jackknife_estimates[[i]])
}


# Jackknife sd estimate by hand
tictoc::tic()
betaJK <- matrix(NA, n, p*d)
for(i in 1:n){
  plnTmp <- PLN(Y[-i, ] ~ -1 + X[-i, ], control=PLN_param(trace=0))
  betaJK[i, ]  <- as.vector(plnTmp$model_par$B)
}
tictoc::toc()
betaSdByHand <- sqrt((n-1)*apply(betaJK, 2, function(x) { mean(x^2) - mean(x)^2 }))

## Check that "automatic" and "manual" jackknife estimates are identical
plot(as.vector(betaJK), as.vector(betaJKPackage),
     ylab = "Computed from package", xlab = "Computed by hand",
     main = "Jackknife estimates (all coefficients, all samples)")
abline(0, 1, col = "red")
max(abs(betaJK - betaJKPackage))

# Empirical sd estimate
B <- n
betaEmp <- matrix(NA, B, p*d)
for(b in 1:B){
  if(b %% 10==0){cat(b, '')}
  Ytmp <- matrix(rpois(n*p, exp(X %*% B0 + rmvnorm(n, sigma= Sigma0))), n, p)
  plnTmp <- PLN(Ytmp ~ -1 + X, control=PLN_param(trace=0))
  betaEmp[b, ]  <- as.vector(plnTmp$model_par$B)
}
cat("\n")
betaSdEmpirical <- apply(betaEmp, 2, sd)

# Comparison
plot(as.data.frame(cbind(betaSdPackage, betaSdByHand, betaSdEmpirical)))

library(lme4)
library(Matrix)
library(PLNmodels)

FitPLNGLMM <- function(X, Y, init=NULL){
  n <- nrow(Y); p <- ncol(Y); d <- ncol(X); site <- 1:n
  # Fit glmmm
  fitList <- list()
  Beta <- matrix(NA, d, p)
  diagSigma <- rep(NA, p)
  for(j in 1:p){
    if(!is.null(init)){
      start <- list(fixef=init$model_par$B[, j], theta=sqrt(init$model_par$Sigma[j, j]))
    } else {
      start <- NULL
    }
    fitList[[j]] <- glmer(Y[, j] ~ -1 + X + (1|site), family='poisson', start=start)
    Beta[, j] <- fitList[[j]]@beta
    diagSigma[j] <- fitList[[j]]@theta^2
  }
  # Estimate Sigma
  A <- exp(X%*%Beta + rep(1, n)%o%diagSigma/2)
  Sigma <- matrix(0, p, p)
  for(i in 1:n){
    Sigma <- Sigma + (Y[i, ] - A[i, ])%o%(Y[i, ] - A[i, ]) /
      (A[i, ] %o% A[i, ])
  }
  Sigma <- Sigma / n
  Sigma <- 1 + Sigma
  return(list(Beta=Beta, Sigma=Sigma, diagSigma=diagSigma, A=A, fitList=fitList))
}


data("oaks")
n <- nrow(oaks$Abundance)
p <- ncol(oaks$Abundance)
x <- model.matrix( ~ 1 + tree, oaks)
y <- oaks$Abundance
X <- as.matrix(bdiag(rep(list(x), p)))
Y <- c(y)

pln <- PLN(Abundance ~ 1 + tree, data = oaks, control = PLN_param(covariance = "diagonal"))
glm <- glm.fit(X, Y, family = poisson())
plot(coef(pln), coef(glm))
abline(0,1)

glmm <- FitPLNGLMM(x, y, init = pln)

plot(coef(pln), glmm$Beta)
abline(0,1)


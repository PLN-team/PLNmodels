library(torch)
library(R6)

log_stirling <- function(n_){
  n_ <- n_+ (n_==0)
  torch_log(torch_sqrt(2*pi*n_)) + n_*log(n_/exp(1))
}

PLN <-
  R6Class("PLN",
          public = list(
            Y = NULL,
            O = NULL,
            X = NULL,
            n = NULL,
            p = NULL,
            d = NULL,
            M = NULL,
            S = NULL,
            A = NULL,
            B = NULL,
            Sigma = NULL,
            Omega = NULL,
            ELBO_list = NULL,

            ## Constructor
            initialize = function(Y, O, X){
              self$Y <- torch_tensor(Y)
              self$O <- torch_tensor(O)
              self$X <- torch_tensor(X)
              self$n <- nrow(Y)
              self$p <- ncol(Y)
              self$d <- ncol(X)
              ## Variational parameters
              self$M <- torch_zeros(self$n, self$p, requires_grad = TRUE)
              self$S <- torch_ones(self$n , self$p, requires_grad = TRUE)
              ## Model parameters
              self$B <- torch_zeros(self$d, self$p, requires_grad = TRUE)
              self$Sigma <- torch_eye(self$p)
              self$Omega <- torch_eye(self$p)
              ## Monitoring
              self$ELBO_list <- c()
            },

            get_Sigma = function(M, S){
              1/self$n * (torch_matmul(torch_transpose(M,2,1),M) + torch_diag(torch_sum(torch_multiply(S,S), dim = 1)))
            },

            get_ELBO = function(B, M, S, Omega){
              S2 <- torch_multiply(S, S)
              XB <- torch_matmul(self$X, B)
              A  <- torch_exp(self$O + M + XB + S2/2)
              elbo <- self$n/2 * torch_logdet(Omega) +
                torch_sum(- A + torch_multiply(self$Y, self$O + M + XB) + .5 * torch_log(S2)) -
                .5 * torch_trace(torch_matmul(torch_matmul(torch_transpose(M, 2, 1), M) + torch_diag(torch_sum(S2, dim = 1)), Omega)) +
                .5 * self$n * self$p - torch_sum(log_stirling(self$Y))
              elbo
            },

            fit = function(N_iter, lr, tol = 1e-8, verbose = FALSE){
              self$ELBO_list <- double(length = N_iter)
              optimizer <- optim_rprop(c(self$B, self$M, self$S), lr = lr)
              objective0 <- Inf
              for (i in 1:N_iter){
                ## reinitialize gradients
                optimizer$zero_grad()

                ## compute current ELBO
                loss <- - self$get_ELBO(self$B, self$M, self$S, self$Omega)

                ## backward propagation and optimization
                loss$backward()
                optimizer$step()

                ## update parameters with close form
                self$Sigma <- self$get_Sigma(self$M, self$S)
                self$Omega <- torch_inverse(self$Sigma)

                objective <- -loss$item()
                if(verbose && (i %% 50 == 0)){
                  pr('i : ', i )
                  pr('ELBO', objective)
                }
                self$ELBO_list[i] <- objective
                if (abs(objective0 - objective)/abs(objective) < tol) {
                  self$ELBO_list <- self$ELBO_list[1:i]
                  break
                } else {
                  objective0 <- objective
                }
              }
            },

            plotLogNegElbo = function(from = 10){
              plot(log(-self$ELBO_list[from:length(self$ELBO_list) ]), type = "l")
            }
          )
  )

myPLN <- PLN$new(Y = oaks$Abundance, O = log(oaks$Offset), X = cbind(1, oaks$tree))
myPLN$fit(30, 0.1, tol = 1e-6)
myPLN$plotLogNegElbo()
